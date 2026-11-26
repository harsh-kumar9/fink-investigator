"""Replay the 2025 catalog cones and the three 2024 TNS registry searches."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import zipfile

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parents[1] / "compass-2025/batch-4/catalog"
catalog = runpy.run_path(str(PRIOR / "check.py"))["catalog"]
catalog.HERE = catalog.catalog.HERE = HERE
archive = catalog.catalog
TNS = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
GUEST = 'tns_marker{"tns_id":0,"type":"guest","name":"Codex astronomy archive query"}'


def tns(cases, fetch):
    receipt = HERE / "tns.csv.json"
    def matches(r):
        a, d = float(r["ra"]), float(r["declination"])
        for c in cases:
            dx, dy = archive.pilot.offsets(c["ra"], c["dec"], a, d)
            dist = math.hypot(dx, dy)
            if dist <= 10 or c["id"].lower() in r["internal_names"].lower():
                yield c["id"], round(dist, 3)
    if fetch and not receipt.exists():
        raw = archive.pilot.curl(TNS, ["-A", GUEST])
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            rows = io.TextIOWrapper(z.open("tns_public_objects.csv"))
            staged = rows.readline().strip()
            reader = csv.DictReader(rows)
            out = io.StringIO()
            writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
            writer.writeheader()
            count = 0
            for row in reader:
                count += 1
                if list(matches(row)):
                    writer.writerow(row)
        excerpt = out.getvalue().encode()
        (HERE / "tns.csv.gz").write_bytes(gzip.compress(excerpt, mtime=0))
        archive.pilot.json_write(receipt, {
            "url": TNS, "staged": staged, "rows_scanned": count,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256_snapshot": hashlib.sha256(raw).hexdigest(),
            "sha256_uncompressed": hashlib.sha256(excerpt).hexdigest(),
        })
    rows = csv.DictReader(io.StringIO(archive.source("tns.csv", TNS, False).decode()))
    out = {c["id"]: [] for c in cases}
    for row in rows:
        for name, distance in matches(row):
            out[name].append({"separation_arcsec": distance} | {k: row[k] for k in
                              ["name_prefix", "name", "internal_names", "type", "discoverydate", "reporting_group"]})
    return out


def run(fetch):
    result = catalog.run(fetch)
    objects = json.loads((HERE.parent / "results/ledger.json").read_text())["objects"]
    cases = [c for c in objects if c["selected_for"]]
    for obj, hits in tns(cases, fetch).items():
        result[obj]["tns_10arcsec_or_alias"] = hits
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
