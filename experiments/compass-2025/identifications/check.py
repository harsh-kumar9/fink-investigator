"""Check dated BTS and SIMBAD, VSX, ATLAS, ASAS-SN identities against the 2025 anchors."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
from xml.etree import ElementTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import pilot
import vsx

CATALOGS = ("B/vsx/vsx", "J/AJ/156/241/table4", "II/366/catv2021")
BTS = "https://sites.astro.caltech.edu/ztf/bts/explorer.php"


def source(name, url, fetch):
    archive = HERE / (name + ".gz")
    meta = HERE / (name + ".json")
    if fetch and not meta.exists():
        raw = pilot.curl(url, [])
        archive.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            meta,
            {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            },
        )
    raw = gzip.decompress(archive.read_bytes())
    receipt = json.loads(meta.read_text())
    if (
        receipt["url"] != url
        or receipt["sha256_uncompressed"] != hashlib.sha256(raw).hexdigest()
    ):
        raise ValueError(f"Unexpected source or checksum for {name}")
    return raw


def target_distance(r, ra, dec):
    return round(math.hypot(*pilot.offsets(ra, dec, r["ra"], r["dec"])), 3)


def csv_rows(raw, columns):
    reader = csv.DictReader(io.StringIO(raw.decode().lstrip()))
    if not set(columns).issubset(reader.fieldnames or ()):
        raise ValueError("Missing CSV columns")
    rows = list(reader)
    if len(rows) >= 100:
        raise ValueError("Too many cone rows")
    return rows


def bts(case, fetch):
    url = (
        BTS
        + "?"
        + urlencode(
            {
                "f": "s",
                "startra": "136.157",
                "endra": "136.163",
                "startdec": "18.464",
                "enddec": "18.469",
                "showdiscdate": "y",
                "showpeakdate": "y",
                "format": "csv",
            }
        )
    )
    rows = csv_rows(
        source("bts.csv", url, fetch),
        ["ZTFID", "IAUID", "RA", "Dec", "discdate", "type"],
    )
    (row,) = [r for r in rows if r["ZTFID"] == case["id"]]
    hours, mins, secs = map(float, row["RA"].split(":"))
    deg, dm, ds = map(float, row["Dec"].split(":"))
    ra = 15 * (hours + mins / 60 + secs / 3600)
    dec = math.copysign(abs(deg) + dm / 60 + ds / 3600, deg)
    return {
        **row,
        "separation_arcsec": target_distance(case, ra, dec),
        "rows_in_box": len(rows),
    }


def catalogs(case, fetch):
    params = {
        "-source": ",".join(CATALOGS),
        "-c": f"{case['ra']} {case['dec']:+}",
        "-c.rs": "5",
        "-out.max": "50",
        "-out.add": "_r",
        "-out": "OID,Name,ATOID,ASASSN-V,RAJ2000,DEJ2000,Type,Class,Period,Per,fp-LSper,fp-lngfitper,ls-Pday",
    }
    raw = source(case["id"] + ".xml", vsx.SERVICE + "?" + urlencode(params), fetch)
    root = ElementTree.fromstring(raw)
    if any(
        i.get("name", "").lower() == "error"
        or i.get("value", "").upper() in {"ERROR", "OVERFLOW"}
        for i in root.findall(".//{*}INFO")
    ):
        raise ValueError("Catalog lookup failed or overflowed")
    result = {t: {"status": "table unavailable"} for t in CATALOGS}
    for table in root.findall(".//{*}TABLE"):
        cols = [f.attrib["name"] for f in table.findall("{*}FIELD")]
        rows = []
        for tr in table.findall(".//{*}TR"):
            r = dict(zip(cols, [(td.text or "").strip() for td in tr], strict=True))
            offset = target_distance(case, float(r["RAJ2000"]), float(r["DEJ2000"]))
            if offset > 5.02 or abs(offset - float(r["_r"])) > 0.1:
                raise ValueError("Catalog position does not match query")
            rows.append({**r, "separation_arcsec": offset})
        if len(rows) >= 50:
            raise ValueError("Catalog truncated")
        result[table.attrib["name"]] = {"status": "returned", "rows": rows}
    return result


def run(fetch):
    ledger = json.loads((HERE.parent / "results/ledger.json").read_text())
    cases = [r for r in ledger["objects"] if r["selected_for"]]
    # Retrieve all SIMBAD neighbors in one request to conserve archive calls.
    clauses = [
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r['ra']},{r['dec']},0.001388889))=1"
        for r in cases
    ]
    query = (
        "SELECT TOP 100 main_id,otype,ra,dec,pmra,pmdec FROM basic WHERE "
        + " OR ".join(clauses)
    )
    url = (
        pilot.SIMBAD
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    )
    neighbors = csv_rows(
        source("simbad.csv", url, fetch), ["main_id", "otype", "ra", "dec"]
    )
    results = {}
    for case in cases:
        local = []
        for r in neighbors:
            distance = target_distance(case, float(r["ra"]), float(r["dec"]))
            if distance <= 5.001:
                local.append({**r, "separation_arcsec": distance})
        info = {"simbad": local}
        if case["selected_for"] == "onset":
            info["bts"] = bts(case, fetch)
        else:
            info["catalogs"] = catalogs(case, fetch)
        results[case["id"]] = info
    return results


def astrometric_case(fetch):
    r = json.loads((HERE.parent / "astrometry/packet.json").read_text())["case"]
    query = (
        "SELECT TOP 100 main_id,otype,ra,dec,pmra,pmdec FROM basic WHERE "
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r['ra']},{r['dec']},0.002777778))=1"
    )
    url = (
        pilot.SIMBAD
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    )
    objects = csv_rows(
        source(r["id"] + ".simbad.csv", url, fetch), ["main_id", "otype", "ra", "dec"]
    )
    for obj in objects:
        obj["separation_arcsec"] = target_distance(
            r, float(obj["ra"]), float(obj["dec"])
        )
    return {
        r["id"]: {"simbad_10_arcsec": objects, "catalogs_5_arcsec": catalogs(r, fetch)}
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--astrometry", action="store_true")
    args = parser.parse_args()
    result = astrometric_case(args.fetch) if args.astrometry else run(args.fetch)
    print(json.dumps(result, indent=2, allow_nan=False))
