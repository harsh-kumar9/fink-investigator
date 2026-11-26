"""Search historical Pan-STARRS dates for the two-filter state change frozen in the protocol."""

import argparse
from collections import Counter, defaultdict
import csv
from itertools import combinations
import io
import json
import math
from statistics import median
from urllib.parse import urlencode

from priors import STARS, archive


def fetch_table(table, key, params, columns, fetch):
    url = (
        "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/"
        + table
        + ".csv?"
        + urlencode({**params, "columns": "[" + ",".join(columns) + "]", "pagesize": 1000})
    )
    reader = csv.DictReader(io.StringIO(archive.source(key + ".csv", url, fetch).decode()))
    rows = list(reader)
    if reader.fieldnames != columns or len(rows) >= 1000:
        raise ValueError("Truncated or incorrect PS1 response")
    return rows


def visit_summary(star, rows, snr, qmin):
    a, d = STARS[star][:2]
    visits = defaultdict(lambda: defaultdict(list))
    for r in rows:
        t, f, e, q = map(float, [r["obsTime"], r["psfFlux"], r["psfFluxErr"], r["psfQfPerfect"]])
        flags = int(r["infoFlag"])
        separation = math.hypot(*archive.pilot.offsets(a, d, float(r["ra"]), float(r["dec"])))
        if (
            not all(map(math.isfinite, (t, f, e, q, separation)))
            or e <= 0
            or f / e < snr
            or q < qmin
            or separation > 1.2
            or flags & (0x8 | 0x80 | 0x100 | 0x400 | 0x800 | 0x1000 | 0x400000)
        ):
            continue
        mag = -2.5 * math.log10(f / 3631)
        err = 2.5 * e / (math.log(10) * f)
        visits[math.floor(t)]["grizy"[int(r["filterID"]) - 1]].append((mag, err))
    nights = {
        day: {f: (median(v[0] for v in obs), max(v[1] for v in obs)) for f, obs in bands.items()}
        for day, bands in visits.items()
    }
    pairs = []
    all_pairs = 0
    close = 0
    for first, last in combinations(sorted(nights), 2):
        bands = sorted(nights[first].keys() & nights[last].keys())
        if len(bands) < 2:
            continue
        all_pairs += 1
        close += last - first <= 10
        delta = {
            b: (nights[last][b][0] - nights[first][b][0], math.hypot(nights[last][b][1], nights[first][b][1]))
            for b in bands
        }
        for direction in [-1, 1]:
            strong = [
                b for b, (dm, err) in delta.items() if direction * dm >= 2 and direction * dm >= 5 * err
            ]
            if len(strong) < 2:
                continue
            pairs.append(
                {
                    "mjd": [first, last],
                    "dt": last - first,
                    "direction": "rise" if direction < 0 else "fade",
                    "changes_mag": {b: [round(delta[b][0], 3), round(delta[b][1], 3)] for b in strong},
                }
            )
    best = {}
    for band in "grizy":
        avail = sorted((m[band][0], m[band][1], day) for day, m in nights.items() if band in m)
        if avail:
            best[band] = {
                "dates": len(avail),
                "brightest": list(avail[0]),
                "faintest": list(avail[-1]),
                "span_mag": round(avail[-1][0] - avail[0][0], 3),
            }
    return {
        "accepted_exposures": sum(len(v) for fs in visits.values() for v in fs.values()),
        "nights": len(nights),
        "days_with_2filters": sum(len(b) >= 2 for b in nights.values()),
        "distinct_date_pairs_with_2_shared_filters": all_pairs,
        "one_to_ten_day_pairs_with_2_shared_filters": close,
        "qualifying_pairs": sorted(pairs, key=lambda p: p["dt"]),
        "single_filter_ranges": best,
    }


def run(fetch=False):
    result = {}
    for n, (a, d, _, _) in STARS.items():
        mean = fetch_table(
            "mean",
            n + "-mean",
            {"ra": a, "dec": d, "radius": 5 / 3600},
            ["objID", "raMean", "decMean", "nDetections"],
            fetch,
        )
        local = [
            dict(
                r,
                separation_arcsec=round(
                    math.hypot(*archive.pilot.offsets(a, d, float(r["raMean"]), float(r["decMean"]))), 3
                ),
            )
            for r in mean
        ]
        out = {}
        for m in local:
            if m["separation_arcsec"] > 1:
                continue
            obj = m["objID"]
            columns = [
                "objID",
                "detectID",
                "filterID",
                "obsTime",
                "ra",
                "dec",
                "psfFlux",
                "psfFluxErr",
                "psfQfPerfect",
                "infoFlag",
                "infoFlag2",
                "infoFlag3",
                "apFlux",
                "psfMajorFWHM",
                "psfMinorFWHM",
            ]
            rows = fetch_table("detection", obj, {"objID": obj}, columns, fetch)
            if (
                len({r["detectID"] for r in rows}) != len(rows)
                or any(r["objID"] != obj for r in rows)
                or len(rows) != int(m["nDetections"])
            ):
                raise ValueError("Mismatched PS1 source or incomplete series")
            out[obj] = {"strict": visit_summary(n, rows, 10, 0.9), "relaxed": visit_summary(n, rows, 5, 0.85)}
        result[n] = {"sources_5arcsec": local, "measurements": out}
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    print(json.dumps(run(p.parse_args().fetch), indent=2))
