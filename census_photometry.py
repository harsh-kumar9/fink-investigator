"""Test epoch-matched direct ZTF photometry under the frozen October protocol."""

import argparse
from bisect import bisect_left
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from statistics import median
from urllib.parse import urlencode

import census
import census_direct_recovery
import pilot
import vsx

SERVICE = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"


def ranks(values):
    positions = defaultdict(list)
    for i, v in enumerate(sorted(values)):
        positions[v].append(i + 1)
    return [sum(positions[v]) / len(positions[v]) for v in values]


def correlation(x, y):
    x = ranks(x)
    y = ranks(y)
    center = (len(x) + 1) / 2
    xx = sum((v - center) ** 2 for v in x)
    yy = sum((v - center) ** 2 for v in y)
    return (
        sum((a - center) * (b - center) for a, b in zip(x, y)) / (xx * yy) ** 0.5
        if xx and yy
        else None
    )


def acquired(case, fetch):
    url = (
        SERVICE
        + "?"
        + urlencode(
            {
                "POS": f"CIRCLE {case['ra']} {case['dec']} 0.0004166667",
                "BANDNAME": "r",
                "BAD_CATFLAGS_MASK": "65535",
                "NOBS_MIN": "30",
                "TIME": f"-Inf {case['jd']-2400000.5:.8f}",
                "FORMAT": "csv",
            }
        )
    )
    folder = census.CENSUS / "sources/direct"
    if fetch:
        folder.mkdir(exist_ok=True, parents=True)
    receipt = folder / (case["id"] + ".json")
    path = folder / (case["id"] + ".csv.gz")
    recovery = folder / "recovery"
    if not receipt.exists() and (recovery / (case["id"] + ".json")).exists():
        raw, recovery_url = census_direct_recovery.recovered(case, recovery, fetch)
        return list(csv.DictReader(io.StringIO(raw.decode()))), recovery_url
    raw = (
        pilot.curl(url, [])
        if fetch and not receipt.exists()
        else gzip.decompress(path.read_bytes())
    )
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if not {"oid", "expid", "mjd", "mag", "magerr", "catflags", "filtercode"}.issubset(
        reader.fieldnames or []
    ):
        raise ValueError(
            f"unexpected IRSA direct-light-curve response for {case['id']}: {raw[:200]!r}"
        )
    rows = list(reader)
    sha = hashlib.sha256(raw).hexdigest()
    if fetch and not receipt.exists():
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
        raise ValueError("Direct measurement provenance mismatch")
    return rows, url


def compare(case, history, fetch, include_pairs=False):
    direct, url = acquired(case, fetch)
    by_oid = defaultdict(dict)
    for row in direct:
        if (
            row["filtercode"] == "zr"
            and int(row["catflags"]) == 0
            and 0 < float(row["magerr"]) <= 0.2
            and float(row["mjd"]) <= case["jd"] - 2400000.5 + 1e-7
        ):
            by_oid[row["oid"]][row["expid"]] = (float(row["mjd"]), float(row["mag"]))
    usable = [
        r
        for r in history[case["id"]]
        if r["i:jd"] <= case["jd"] + 1e-7
        and r["i:fid"] == 2
        and (r["i:drb"] or 0) >= 0.8
        and 0 < (r["i:sigmapsf"] or 99) <= 0.2
    ]

    def match_series(entries):
        points = sorted(entries.values())
        epochs = [m for m, _ in points]
        edges = []
        for i, r in enumerate(usable):
            mjd = r["i:jd"] - 2400000.5
            j = bisect_left(epochs, mjd - 60 / 86400)
            while j < len(epochs) and epochs[j] <= mjd + 60 / 86400:
                edges.append((abs(mjd - epochs[j]), i, j))
                j += 1
        used_alert = set()
        used_direct = set()
        groups = defaultdict(list)
        for _, i, j in sorted(edges):
            if i in used_alert or j in used_direct:
                continue
            used_alert.add(i)
            used_direct.add(j)
            obs = usable[i]
            sign = 1 if str(obs["i:isdiffpos"]).lower() in {"t", "true", "1"} else -1
            diff = sign * 10 ** (-0.4 * (obs["i:magpsf"] - 20))
            directmag = points[j][1]
            groups[pilot.day_of(obs["i:jd"])].append(
                (diff, 10 ** (-0.4 * (directmag - 20)))
            )
        return {
            day: (median(x for x, _ in vals), median(y for _, y in vals))
            for day, vals in groups.items()
        }

    checked = {
        oid: match_series(entries)
        for oid, entries in by_oid.items()
        if len(entries) >= 30
    }
    source = min(checked, key=lambda k: (-len(checked[k]), k)) if checked else None
    dated = checked.get(source, {})
    pairs = [dated[day] for day in sorted(dated)]
    x = [a for a, b in pairs]
    y = [b for a, b in pairs]
    pos = [20 - 2.5 * math.log10(b) for a, b in pairs if a > 0]
    neg = [20 - 2.5 * math.log10(b) for a, b in pairs if a < 0]
    rho = correlation(x, y) if len(x) >= 20 else None
    gap = median(neg) - median(pos) if min(len(pos), len(neg)) >= 5 else None
    measurement = {
        "id": case["id"],
        "query": url,
        "direct_source_id": source,
        "eligible_series_count": len(checked),
        "first_matched_date": str(min(dated)) if dated else None,
        "last_matched_date": str(max(dated)) if dated else None,
        "matched_dates": len(pairs),
        "positive_dates": len(pos),
        "negative_dates": len(neg),
        "spearman": round(rho, 3) if rho is not None else None,
        "negative_minus_positive_median_direct_mag": (
            round(gap, 3) if gap is not None else None
        ),
        "threshold_supported": rho is not None
        and rho >= 0.7
        and gap is not None
        and gap > 0,
    }

    if include_pairs:
        measurement["paired_medians"] = pairs
    return measurement


def run(fetch):
    prev, cases = census.assess(False)
    _, anchors, _ = census.locked_cohort()
    _, raw = pilot.snapshot(census.BASE / "snapshot")
    h = pilot.histories(
        [r for n, v in raw.items() if n.startswith("history-") for r in v]
    )
    control = next(r for r in prev["objects"] if r["id"] == "ZTF18achfhai")
    calibration = compare(control, h, fetch)
    history = pilot.histories(census.queries("fink", anchors, False))
    targets = [
        r
        for r in cases
        if r.get("adjudication")
        and not vsx.closed_alias(r)
        and r["id"] != "ZTF26acbcwpz"
    ]
    with ThreadPoolExecutor(max_workers=2) as workers:
        scored = list(workers.map(lambda r: compare(r, history, fetch), targets))
    return {
        "control": calibration,
        "calibration_passed": calibration["threshold_supported"],
        "exploration": scored,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    res = run(args.fetch)
    pilot.json_write(census.CENSUS / "results/direct.json", res)
    for r in [res["control"]] + res["exploration"]:
        print(
            r["id"],
            r["matched_dates"],
            r["spearman"],
            r["negative_minus_positive_median_direct_mag"],
            r["threshold_supported"],
        )
