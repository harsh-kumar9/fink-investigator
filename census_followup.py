"""Check open October objects across four VizieR tables and archive multi-epoch images."""

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import census
import cutouts
import motion
import pilot
import vsx

TABLES = [
    "I/355/gaiadr3",
    "I/358/vclassre",
    "J/ApJS/249/18/table2",
    "J/ApJS/249/18/table3",
]
FIELDS = "Source,RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag,Class,ClassSc,ID,RAJ2000,DEJ2000,Per,Per-g,Per-r,Type"


def catalogs(case, fetch):
    url = (
        vsx.SERVICE
        + "?"
        + urlencode(
            {
                "-source": ",".join(TABLES),
                "-c": f"{case['ra']} {case['dec']:+}",
                "-c.rs": "5",
                "-out": FIELDS,
                "-out.add": "_r",
                "-out.max": "50",
            }
        )
    )
    folder = census.CENSUS / "sources/extended"
    if fetch:
        folder.mkdir(parents=True, exist_ok=True)
    receipt = folder / (case["id"] + ".json")
    path = folder / (case["id"] + ".xml.gz")
    if not receipt.exists() and fetch:
        raw = pilot.curl(url, [])
    else:
        raw = gzip.decompress(path.read_bytes())
    sha = hashlib.sha256(raw).hexdigest()
    root = ET.fromstring(raw)
    if any(
        i.get("value", "").upper() in {"ERROR", "OVERFLOW"}
        or i.get("name", "").lower() == "error"
        for i in root.findall(".//{*}INFO")
    ):
        raise ValueError("VizieR lookup failed")
    results = {}
    for table in root.findall(".//{*}TABLE"):
        name = table.get("name")
        fields = [f.get("name") for f in table.findall("{*}FIELD")]
        rows = [
            dict(
                zip(
                    fields,
                    [(td.text or "").strip() for td in r.findall("{*}TD")],
                    strict=True,
                )
            )
            for r in table.findall(".//{*}TR")
        ]
        if len(rows) >= 50 or any(not 0 <= float(r["_r"]) <= 5.01 for r in rows):
            raise ValueError("VizieR lookup exceeded bounds")
        if name == TABLES[0]:
            for r in rows:
                x, y = pilot.offsets(
                    case["ra"], case["dec"], float(r["RA_ICRS"]), float(r["DE_ICRS"])
                )
                elapsed = 2000 + (case["jd"] - 2451545.0) / 365.25 - 2016
                if r["pmRA"] and r["pmDE"]:
                    x -= elapsed * float(r["pmRA"]) / 1000
                    y -= elapsed * float(r["pmDE"]) / 1000
                    r["at_alert_epoch_arcsec"] = round(math.hypot(x, y), 3)
                else:
                    r["at_alert_epoch_arcsec"] = None
                r["j2016_separation_arcsec"] = round(
                    math.hypot(
                        *pilot.offsets(
                            case["ra"],
                            case["dec"],
                            float(r["RA_ICRS"]),
                            float(r["DE_ICRS"]),
                        )
                    ),
                    3,
                )
        results[name] = sorted(rows, key=lambda r: float(r["_r"]))
    if set(results) != set(TABLES):
        raise ValueError("Incomplete extended-catalog response")
    if not receipt.exists() and fetch:
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256_uncompressed": sha,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    m = json.loads(receipt.read_text())
    if m["url"] != url or m["sha256_uncompressed"] != sha:
        raise ValueError("Extended catalog checksum or request mismatch")
    return {"id": case["id"], "query": url, "tables": results}


def run(fetch):
    _, cases = census.assess(False)
    _, anchors, _ = census.locked_cohort()
    history = pilot.histories(census.queries("fink", anchors, False))
    gaia = census.gaia(False)
    results = []
    for r in cases:
        need_catalogs = r.get("adjudication") and not vsx.closed_alias(r)
        test_motion = (
            r["status"] == "moving-star alternative" and len(r["gaia_matches"]) == 1
        )
        if not (need_catalogs or test_motion):
            continue
        if need_catalogs:
            results.append(catalogs(r, fetch))
        folder = census.CENSUS / "review" / r["id"]
        if test_motion:
            star = gaia[r["gaia_matches"][0]["source_id"]]
            motion.inspect_match(r, star, folder / "reference", fetch)
        gallery = folder / "images"
        if not gallery.exists() and fetch:
            cutouts.capture_images(r, history[r["id"]], gallery)
        archived = json.loads((gallery / "manifest.json").read_text())
        expected = cutouts.epochs(history[r["id"]], r)
        if [str(row["i:candid"]) for row in expected] != [
            e["candid"] for e in archived["epochs"]
        ]:
            raise ValueError("image selection differs from frozen epochs")
        for e in archived["epochs"]:
            if len(e["images"]) != 3:
                raise ValueError("incomplete image triplet")
            for im in e["images"]:
                if (
                    hashlib.sha256((gallery / im["file"]).read_bytes()).hexdigest()
                    != im["sha256"]
                ):
                    raise ValueError("image bytes differ from capture")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    extra = run(args.fetch)
    pilot.json_write(census.CENSUS / "results/extended.json", extra)
    for r in extra:
        print(r["id"], [(name, len(v)) for name, v in r["tables"].items()])
