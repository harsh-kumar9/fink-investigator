"""Capture or replay a calendar year of Fink's nightly anomaly selections."""

import argparse
from collections import Counter
from datetime import date, datetime
import json
import math
import os
from pathlib import Path

import pilot

VERSION = 1
COLUMNS = pilot.COLUMNS + ",i:jdstarthist"
ARCHIVE_CAP = 400  # greater than the documented 10/night for every month
STRATA = ("young", "other_positive", "other")


def months(year):
    for month in range(1, 13):
        start = date(year, month, 1)
        stop = date(year + (month == 12), month % 12 + 1, 1)
        yield start, stop


def parent_request(start, stop):
    return {
        "start_date": str(start),
        "stop_date": str(stop),
        "n": ARCHIVE_CAP,
        "columns": COLUMNS + ",d:anomaly_score",
        "output-format": "json",
    }


def history_requests(anchors):
    names = sorted(pilot.oid(a) for a in anchors)
    for i in range(0, len(names), 4):
        yield f"history-{i:03d}.json.gz", {
            "objectId": ",".join(names[i : i + 4]),
            "columns": COLUMNS,
            "output-format": "json",
        }


def checked_histories(anchors, response):
    """Retry missing anchors by ID. An unresolved anchor yields an explicit abstention."""
    histories = {}
    retries = {}
    unavailable = []
    for filename, params in history_requests(anchors):
        batch = pilot.histories(response(filename, "/objects", params))
        requested = set(params["objectId"].split(","))
        if set(batch) - requested:
            raise ValueError(f"Unexpected object in {filename}")
        for name in requested:
            histories[name] = batch.get(name, [])
    for anchor in anchors:
        name = pilot.oid(anchor)
        candid = anchor["i:candid"]
        rows = histories[name]
        hits = [row for row in rows if row["i:candid"] == candid]
        if not hits:
            params = {"objectId": name, "columns": COLUMNS, "output-format": "json"}
            filename = f"retry-{name}.json.gz"
            single = pilot.histories(response(filename, "/objects", params))
            if set(single) - {name}:
                raise ValueError(f"Unexpected object in {filename}")
            batch_count = len(rows)
            rows = histories[name] = single.get(name, [])
            hits = [row for row in rows if row["i:candid"] == candid]
            retries[name] = {
                "filename": filename,
                "candid": str(candid),
                "batch_rows": batch_count,
                "single_id_rows": len(rows),
                "single_id_rows_through_anchor": sum(
                    r["i:jd"] <= anchor["i:jd"] for r in rows
                ),
                "recovered": bool(hits),
            }
            if not hits:
                unavailable.append(name)
                histories[name] = []
                continue
        if len(hits) != 1 or any(
            hits[0][k] != anchor[k] for k in COLUMNS.split(",") if k.startswith("i:")
        ):
            raise ValueError(f"Anchor measurement differs between endpoints: {name}")
    return histories, retries, unavailable


def stratum(row):
    first = row["i:jdstarthist"]
    if first is not None and 0 <= row["i:jd"] - first <= 90:
        return "young"
    if str(row["i:isdiffpos"]).lower() in {"t", "true", "1"}:
        return "other_positive"
    return "other"


def spaced(rows, count):
    """Use the pilot's endpoint-inclusive evenly spaced indices."""
    rows = sorted(rows, key=lambda r: (r["i:jd"], pilot.oid(r), r["i:candid"]))
    k = min(count, len(rows))
    if k < 2:
        return rows[:k]
    return [rows[round(i * (len(rows) - 1) / (k - 1))] for i in range(k)]


def select(year, monthly, quota):
    seen = set()
    all_returned_ids = set()
    in_year_candids = set()
    anchors = []
    summaries = []
    for start, stop in months(year):
        raw = monthly[start.month]
        if len(raw) >= ARCHIVE_CAP:
            raise ValueError(f"Response limit reached for {start}")
        best = {}
        by_date = Counter()
        packets = set()
        for row in raw:
            name = pilot.oid(row)
            all_returned_ids.add(name)
            if not math.isfinite(row["i:jd"]):
                raise ValueError(f"Invalid observation time for {name}")
            if not (pilot.utc_jd(start) <= row["i:jd"] < pilot.utc_jd(stop)):
                continue
            candid = row["i:candid"]
            if not isinstance(candid, int) or candid <= 0 or candid in in_year_candids:
                raise ValueError(f"Missing or repeated in-window candid in {start}")
            in_year_candids.add(candid)
            packets.add(candid)
            by_date[str(pilot.day_of(row["i:jd"]))] += 1
            if not math.isfinite(row["d:anomaly_score"]):
                raise ValueError(f"Nonfinite score for {name}")
            key = (row["d:anomaly_score"], row["i:jd"], candid)
            if name not in best or key < best[name][0]:
                best[name] = key, row
        pool = [r for name, (_, r) in best.items() if name not in seen]
        seen.update(best)
        chosen = []
        for group in STRATA:
            chosen.extend(spaced([r for r in pool if stratum(r) == group], quota))
        taken = {pilot.oid(r) for r in chosen}
        remainder = [r for r in pool if pilot.oid(r) not in taken]
        chosen.extend(spaced(remainder, 3 * quota - len(chosen)))
        chosen = spaced(chosen, len(chosen))
        anchors.extend(chosen)
        summaries.append(
            {
                "start": str(start),
                "stop_exclusive": str(stop),
                "raw_packets": len(raw),
                "in_window_packets": len(packets),
                "out_of_window_packets": len(raw) - len(packets),
                "distinct_ids": len(best),
                "new_ids": len(pool),
                "selected_ids": len(chosen),
                "utc_dates": dict(sorted(by_date.items())),
                "new_by_stratum": {
                    g: sum(stratum(r) == g for r in pool) for g in STRATA
                },
                "selected_by_stratum": {
                    g: sum(stratum(r) == g for r in chosen) for g in STRATA
                },
            }
        )
    coverage = {
        k: sum(row[k] for row in summaries)
        for k in (
            "raw_packets",
            "in_window_packets",
            "out_of_window_packets",
            "selected_ids",
        )
    }
    coverage.update(
        year=year,
        per_stratum=quota,
        raw_distinct_ids=len(all_returned_ids),
        distinct_ids=len(seen),
        dates=sum(len(row["utc_dates"]) for row in summaries),
        monthly=summaries,
    )
    return anchors, coverage


def config(year, quota):
    if not 1900 <= year < 9999 or not 1 <= quota <= 4:
        raise ValueError("Choose a calendar year and a per-stratum quota from 1 to 4")
    return {"version": VERSION, "year": year, "per_stratum": quota}


def load(snapshot):
    """Return raw selected anchors, full per-object histories, and calendar coverage, offline."""
    folder = Path(snapshot)
    manifest, raw = pilot.snapshot(folder)
    expected = config(manifest["year"], manifest["per_stratum"])
    if any(manifest[k] != v for k, v in expected.items()):
        raise ValueError("Unsupported compass snapshot configuration")
    sources = manifest["sources"]
    verified = set()

    def source(name, endpoint, params):
        receipt = sources[name]
        if (
            receipt["url"] != pilot.FINK + endpoint
            or receipt["method"] != "POST"
            or receipt["params"] != params
            or not isinstance(raw[name], list)
        ):
            raise ValueError(f"Request mismatch: {name}")
        if datetime.fromisoformat(receipt["retrieved_utc"]).utcoffset() is None:
            raise ValueError(f"Missing retrieval timezone: {name}")
        verified.add(name)
        return raw[name]

    monthly = {
        start.month: source(
            f"anomaly-{start:%Y-%m}.json.gz", "/anomaly", parent_request(start, stop)
        )
        for start, stop in months(expected["year"])
    }
    anchors, coverage = select(expected["year"], monthly, expected["per_stratum"])
    histories, retries, unavailable = checked_histories(anchors, source)
    coverage.update(
        history_retries=retries,
        unavailable_history_ids=unavailable,
        verified_histories=len(anchors) - len(unavailable),
    )
    if verified != set(sources):
        raise ValueError("Unexpected response in compass manifest")
    if manifest.get("anchors") != [[pilot.oid(a), str(a["i:candid"])] for a in anchors]:
        raise ValueError("Selected anchor identities differ from manifest")
    return anchors, histories, coverage


def fetch(folder, year, quota):
    folder = Path(folder)
    expected = config(year, quota)
    if folder.exists():
        raise ValueError(
            f"Completed destination already exists: {folder}. Use replay or a new location."
        )
    partial = folder.with_name(folder.name + ".partial")
    with pilot.locked_directory(partial):
        manifest_path = partial / "manifest.json"
        manifest = (
            json.loads(manifest_path.read_text())
            if manifest_path.exists()
            else {**expected, "sources": {}}
        )
        if any(manifest.get(k) != v for k, v in expected.items()):
            raise ValueError("Partial snapshot has different selection settings")
        monthly = {}
        for start, stop in months(year):
            monthly[start.month] = pilot.store_raw(
                partial,
                manifest,
                f"anomaly-{start:%Y-%m}.json.gz",
                pilot.FINK + "/anomaly",
                parent_request(start, stop),
            )
        anchors, _ = select(year, monthly, quota)
        manifest["anchors"] = [[pilot.oid(a), str(a["i:candid"])] for a in anchors]

        def capture(name, endpoint, params):
            return pilot.store_raw(
                partial, manifest, name, pilot.FINK + endpoint, params
            )

        checked_histories(anchors, capture)
        checkpoint = partial / "manifest.tmp"
        pilot.json_write(checkpoint, manifest)
        os.replace(checkpoint, manifest_path)
        result = load(partial)
        os.rename(partial, folder)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    capture = sub.add_parser(
        "fetch",
        help="Download a new calendar-year snapshot; interrupted downloads resume",
    )
    capture.add_argument("--year", type=int, required=True)
    capture.add_argument("--per-stratum", type=int, default=4, choices=range(1, 5))
    for p in (
        capture,
        sub.add_parser("replay", help="Verify and count a frozen snapshot offline"),
    ):
        p.add_argument("--snapshot", type=Path, required=True)
    args = parser.parse_args()
    result = (
        fetch(args.snapshot, args.year, args.per_stratum)
        if args.command == "fetch"
        else load(args.snapshot)
    )
    print(json.dumps(result[2], indent=2))
