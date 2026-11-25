"""Replay the complete October anomaly archive. --fetch fills missing source snapshots."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timezone
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import pilot
import vsx

BASE = pilot.ROOT / "experiments/replication/2026-10"
CENSUS = pilot.ROOT / "experiments/october-census"
GAIA_QUERY = (
    "SELECT TOP 2000 source_id,ra,dec,ref_epoch,pmra,pmdec,pmra_error,pmdec_error FROM gaia_dr3_source "
    "WHERE CONTAINS(POINT(ra,dec),CIRCLE(TAP_UPLOAD.targets.my_ra,TAP_UPLOAD.targets.my_dec,0.0166666667))=1 "
    "AND (pmra*pmra+pmdec*pmdec)>=22500"
)


def locked_cohort():
    result = pilot.reconcile_snapshot(BASE / "snapshot")
    if not result["gate_passed"]:
        raise ValueError("The original control gate did not pass")
    _, raw = pilot.snapshot(BASE / "snapshot")
    rows = raw["anomaly.json.gz"]
    all_anchors = pilot.select_alerts(
        rows, len(rows), date(2026, 10, 1), date(2026, 11, 1)
    )
    previous = {r["id"] for r in result["objects"]}
    anchors = [r for r in all_anchors if pilot.oid(r) not in previous]
    manifest = json.loads((CENSUS / "cohort.json").read_text())
    if [
        {"id": pilot.oid(r), "candid": str(r["i:candid"])} for r in anchors
    ] != manifest:
        raise ValueError("Cohort differs from the registered object and packet IDs")
    old_history = pilot.histories(
        [r for n, v in raw.items() if n.startswith("history-") for r in v]
    )
    controls = [
        next(r for r in old_history[c["id"]] if str(r["i:candid"]) == c["candid"])
        for c in result["calibration"]
    ]
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["target", "my_ra", "my_dec"])
    w.writerows([[pilot.oid(r), r["i:ra"], r["i:dec"]] for r in anchors + controls])
    if buf.getvalue() != (CENSUS / "gaia_cones.csv").read_text():
        raise ValueError("Gaia cone centers differ from the registered anchors")
    return result, anchors, controls


def queries(provider, anchors, fetch):
    folder = CENSUS / "sources" / provider
    if fetch:
        folder.mkdir(parents=True, exist_ok=True)
    location = folder / "manifest.json"
    manifest = (
        json.loads(location.read_text()) if location.exists() else {"sources": {}}
    )
    ids = sorted(pilot.oid(r) for r in anchors)
    size = 4 if provider == "fink" else 8
    rows = []
    for i in range(0, len(anchors), size):
        if provider == "fink":
            name = f"history-{i:03}.json.gz"
            url = pilot.FINK + "/objects"
            params = {
                "objectId": ",".join(ids[i : i + size]),
                "columns": pilot.COLUMNS,
                "output-format": "json",
            }
        else:
            name = f"simbad-{i:03}.csv.gz"
            url = pilot.SIMBAD
            params = {
                "QUERY": "SELECT TOP 2000 main_id,otype,ra,dec FROM basic WHERE "
                + pilot.cone_clause(anchors[i : i + size], 3)
            }
        if name not in manifest["sources"] and not fetch:
            raise ValueError(f"missing archived source {provider}/{name}")
        batch = pilot.store_raw(
            folder, manifest, name, url, params, tap=provider == "simbad"
        )
        rows.extend(batch)
    return rows


def gaia(fetch):
    folder = CENSUS / "sources/gaia"
    if fetch:
        folder.mkdir(parents=True, exist_ok=True)
    raw_path = folder / "matches.csv.gz"
    receipt = folder / "manifest.json"
    upload = CENSUS / "gaia_cones.csv"
    request = {
        "url": pilot.IRSA,
        "method": "POST multipart/form-data",
        "QUERY": GAIA_QUERY,
        "FORMAT": "CSV",
        "UPLOAD": "targets,param:coords",
        "coords_sha256": hashlib.sha256(upload.read_bytes()).hexdigest(),
    }
    if not receipt.exists() and fetch:
        raw = pilot.curl(
            pilot.IRSA,
            [
                "-F",
                "UPLOAD=targets,param:coords",
                "-F",
                f"coords=@{upload}",
                "-F",
                "FORMAT=CSV",
                "-F",
                f"QUERY={GAIA_QUERY}",
            ],
        )
        pilot.parse_csv(raw)
        raw_path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "request": request,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            },
        )
    m = json.loads(receipt.read_text())
    raw = gzip.decompress(raw_path.read_bytes())
    if (
        request != m["request"]
        or hashlib.sha256(raw).hexdigest() != m["sha256_uncompressed"]
    ):
        raise ValueError("Gaia multi-cone request or response changed")
    return {r["source_id"]: r for r in pilot.parse_csv(raw)}


def assess(fetch):
    previous, anchors, controls = locked_cohort()
    # The providers are independent. Each sees one request at a time.
    with ThreadPoolExecutor(max_workers=3) as workers:
        phot = workers.submit(queries, "fink", anchors, fetch)
        cats = workers.submit(queries, "simbad", anchors, fetch)
        stars = workers.submit(gaia, fetch)
        history = pilot.histories(phot.result())
        simbad = {r["main_id"]: r for r in cats.result()}
        gaia_rows = stars.result()
    for c, r in zip(previous["calibration"], controls, strict=True):
        if pilot.stellar_matches(r, gaia_rows) != c["matches"]:
            raise ValueError(
                f"Gaia multi-cone screen disagrees with the frozen control {c['id']}"
            )
    cases = pilot.screen_targets(
        anchors, history, gaia_rows, simbad, previous["calibration"]
    )
    vsx.annotate_cases(
        CENSUS / "sources/vsx",
        [
            r
            for r in cases
            if r["status"] == "review"
            or (r["status"] == "broker crossmatch failed" and not r["simbad_neighbors"])
        ],
        (lambda url: pilot.curl(url, [])) if fetch else None,
    )
    return previous, cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    with pilot.locked_directory(CENSUS):
        previous, cases = assess(args.fetch)
        dest = CENSUS / "results"
        dest.mkdir(exist_ok=True)
        pilot.json_write(dest / "remainder.json", cases)
        print(dict(Counter(r["status"] for r in cases)))
        for r in cases:
            if r.get("adjudication") and not vsx.closed_alias(r):
                print(r["id"], r["status"], r["adjudication"], r["ra"], r["dec"])
