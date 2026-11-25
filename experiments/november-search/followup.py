"""Repeat the October catalog, imaging and direct photometry checks on audit survivors."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import census
import census_followup
import census_photometry
import cutouts
import pilot
import vsx


def histories(folder):
    _, raw = pilot.snapshot(folder)
    return pilot.histories(
        [r for name, pack in raw.items() if name.startswith("history-") for r in pack]
    )


def check_images(case, history, fetch):
    folder = HERE / "review" / case["id"]
    rows = history[case["id"]]
    if fetch and not folder.exists():
        cutouts.capture_images(case, rows, folder)
    images = json.loads((folder / "manifest.json").read_text())
    expected = [str(r["i:candid"]) for r in cutouts.epochs(rows, case)]
    if [e["candid"] for e in images["epochs"]] != expected:
        raise ValueError("Unexpected image epochs")
    for e in images["epochs"]:
        if [im["kind"] for im in e["images"]] != ["Science", "Template", "Difference"]:
            raise ValueError("Incomplete image triplet")
        for im in e["images"]:
            raw = (folder / im["file"]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != im["sha256"]:
                raise ValueError("Image checksum mismatch")


def variability_catalogs(case, fetch):
    # An exploratory prior-art check in two independent optical surveys.
    tables = ["J/AJ/156/241/table4", "II/366/catv2021"]
    url = (
        vsx.SERVICE
        + "?"
        + urlencode(
            {
                "-source": ",".join(tables),
                "-c": f"{case['ra']} {case['dec']:+}",
                "-c.rs": "5",
                "-out.add": "_r",
                "-out.max": "50",
                "-out": "ATOID,ID,ASASSN-V,RAJ2000,DEJ2000,Class,Type,ddc,prox,ls-Pday,fp-period,Per",
            }
        )
    )
    folder = HERE / "sources/prior-variables"
    receipt = folder / (case["id"] + ".json")
    path = folder / (case["id"] + ".xml.gz")
    new = not receipt.exists() and fetch
    if new:
        folder.mkdir(exist_ok=True, parents=True)
        raw = pilot.curl(url, [])
    else:
        raw = gzip.decompress(path.read_bytes())
    tree = ET.fromstring(raw)
    if any(
        i.get("value", "").upper() in {"ERROR", "OVERFLOW"}
        or i.get("name", "").lower() == "error"
        for i in tree.findall(".//{*}INFO")
    ):
        raise ValueError("Prior variability catalog error")
    matches = {}
    for table in tree.findall(".//{*}TABLE"):
        cols = [f.get("name") for f in table.findall("{*}FIELD")]
        rows = [
            dict(
                zip(
                    cols,
                    [(s.text or "").strip() for s in tr.findall("{*}TD")],
                    strict=True,
                )
            )
            for tr in table.findall(".//{*}TR")
        ]
        if len(rows) >= 50:
            raise ValueError("Prior variability catalog overflow")
        for r in rows:
            sep = math.hypot(
                *pilot.offsets(
                    case["ra"], case["dec"], float(r["RAJ2000"]), float(r["DEJ2000"])
                )
            )
            if float(r["_r"]) > 5.01 or abs(sep - float(r["_r"])) > 0.1:
                raise ValueError("Catalog position inconsistent")
        matches[table.get("name")] = rows
    if set(matches) != set(tables):
        raise ValueError("Missing prior variability table")
    digest = hashlib.sha256(raw).hexdigest()
    if new:
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256_uncompressed": digest,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    stored = json.loads(receipt.read_text())
    if stored["url"] != url or stored["sha256_uncompressed"] != digest:
        raise ValueError("Prior variability checksums or query changed")
    return {"id": case["id"], "query": url, "tables": matches}


def run(fetch):
    # Re-establish the direct photometry calibration before redirecting the shared
    # archive helpers to this experiment's own source directory.
    control_snapshot = census.BASE / "snapshot"
    prior = pilot.evaluate(control_snapshot)
    control = next(r for r in prior["objects"] if r["id"] == "ZTF18achfhai")
    calibration = census_photometry.compare(control, histories(control_snapshot), False)
    if not prior["gate_passed"] or not calibration["threshold_supported"]:
        raise ValueError("October control gate failed")
    snapshot = HERE / "snapshot"
    current = pilot.reconcile_snapshot(snapshot)
    if not current["gate_passed"]:
        raise ValueError("November control gate failed")
    ids = json.loads((HERE / "results/audit.json").read_text())[
        "audit_quality_survivors"
    ][:3]
    cases = [next(r for r in current["objects"] if r["id"] == oid) for oid in ids]
    history = histories(snapshot)
    archive = census.CENSUS
    census.CENSUS = HERE
    try:
        # Each service receives its small list sequentially. The three services are independent.
        with ThreadPoolExecutor(max_workers=3) as pool:
            cats = pool.submit(
                lambda: [census_followup.catalogs(c, fetch) for c in cases]
            )
            images = pool.submit(
                lambda: [check_images(c, history, fetch) for c in cases]
            )
            phot = pool.submit(
                lambda: [census_photometry.compare(c, history, fetch) for c in cases]
            )
            pilot.json_write(HERE / "results/extended.json", cats.result())
            pilot.json_write(
                HERE / "results/prior-variables.json",
                [variability_catalogs(c, fetch) for c in cases],
            )
            images.result()
            pilot.json_write(
                HERE / "results/direct.json",
                {"control": calibration, "objects": phot.result()},
            )
    finally:
        census.CENSUS = archive


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    run(parser.parse_args().fetch)
