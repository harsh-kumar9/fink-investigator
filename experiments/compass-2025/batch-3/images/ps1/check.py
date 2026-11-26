"""Check Pan-STARRS objects and detections before ASASSN-25ba's 2025 outburst."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import median
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
YEAR = HERE.parents[2]
sys.path.insert(0, str(YEAR.parents[1]))
import compass_archive
import compass_direct
import pilot

NAME = "ZTF25aaohzdc"
BASE = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/"
BANDS = "grizy"


def query(table, key, params, columns, fetch):
    uri = (
        BASE
        + table
        + ".csv?"
        + urlencode(
            {**params, "columns": "[" + ",".join(columns) + "]", "pagesize": 1000}
        )
    )
    records = compass_direct.capture(HERE / "sources", key, uri, set(columns), fetch)
    if len(records) >= 1000:
        raise ValueError("Pan-STARRS result hit the request limit")
    return records


def forced_warps(obj, fetch):
    columns = [
        "objid",
        "uniquepspsfwid",
        "forcedwarpid",
        "detectid",
        "filterid",
        "obstime",
        "fpsfflux",
        "fpsffluxerr",
        "fpsfqfperfect",
        "fpsfchisq",
        "finfoflag",
        "finfoflag2",
        "finfoflag3",
    ]
    query_text = (
        "SELECT TOP 1000 "
        + ",".join(columns)
        + " FROM ps1_dr2.forced_warp_measurement WHERE objid="
        + obj
        + " ORDER BY obstime,uniquepspsfwid"
    )
    url = "https://mast.stsci.edu/vo-tap/api/v0.1/mast_catalogs/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query_text}
    )
    rows = compass_direct.capture(
        HERE / "sources", obj + "-forced", url, set(columns), fetch
    )
    if (
        len(rows) >= 1000
        or any(r["objid"] != obj for r in rows)
        or len({r["uniquepspsfwid"] for r in rows}) != len(rows)
    ):
        raise ValueError("Unexpected forced photometry rows")
    filters = defaultdict(lambda: defaultdict(list))
    for r in rows:
        t, f, e, q = (
            float(r[k]) for k in ("obstime", "fpsfflux", "fpsffluxerr", "fpsfqfperfect")
        )
        if all(map(math.isfinite, (t, f, e, q))) and e > 0 and f != -999 and q >= 0.85:
            band = BANDS[int(r["filterid"]) - 1]
            filters[band][str(pilot.day_of(t + 2400000.5))].append((f, e))
    report = {}
    for band, days in filters.items():
        medians, limits = [], []
        for day, values in sorted(days.items()):
            bright = [(f, e) for f, e in values if f / e >= 5]
            if bright:
                medians.append(
                    (
                        median(-2.5 * math.log10(f / 3631) for f, e in bright),
                        median(2.5 * e / (math.log(10) * f) for f, e in bright),
                        day,
                    )
                )
            else:
                limits.append(
                    median(-2.5 * math.log10(5 * e / 3631) for f, e in values)
                )
        report[band] = {
            "quality_dates": len(days),
            "date_range": [min(days), max(days)],
            "dates_with_SNR5": len(medians),
            "dates_below_SNR5": len(limits),
            "median_5sigma_depth_on_nondetection_dates": (
                round(median(limits), 3) if limits else None
            ),
        }
        if medians:
            lo, hi = min(medians), max(medians)
            report[band].update(
                brightest=[round(lo[0], 3), round(lo[1], 3), lo[2]],
                faintest=[round(hi[0], 3), round(hi[1], 3), hi[2]],
                selected=(
                    len(medians) >= 3
                    and hi[0] - lo[0] >= 1
                    and hi[0] - lo[0] >= 5 * math.hypot(lo[1], hi[1])
                ),
            )
    return {
        "raw_exposures": len(rows),
        "quality_exposures": sum(
            len(values) for nights in filters.values() for values in nights.values()
        ),
        "unique_quality_dates": len(
            {day for nights in filters.values() for day in nights}
        ),
        "bands": report,
    }


def run(fetch=False):
    anchors, _, _ = compass_archive.load(YEAR / "snapshot")
    a = next(v for v in anchors if pilot.oid(v) == NAME)
    params = {"ra": a["i:ra"], "dec": a["i:dec"], "radius": 10 / 3600}

    def sep(ra, dec):
        return math.hypot(*pilot.offsets(float(ra), float(dec), a["i:ra"], a["i:dec"]))

    mean = query(
        "mean",
        "mean",
        params,
        [
            "objID",
            "raMean",
            "decMean",
            "nDetections",
            *[b + "MeanPSFMag" for b in BANDS],
        ],
        fetch,
    )
    stack = query(
        "stack",
        "stack",
        params,
        [
            "objID",
            "raStack",
            "decStack",
            "nDetections",
            "primaryDetection",
            "bestDetection",
            *[b + "PSFMag" for b in BANDS],
        ],
        fetch,
    )
    matches = {}
    for rows, ra, dec, source in [
        (mean, "raMean", "decMean", "mean"),
        (stack, "raStack", "decStack", "stack"),
    ]:
        for r in rows:
            separation = sep(r[ra], r[dec])
            if separation <= 1.5:
                key = r["objID"]
                record = matches.setdefault(key, {"mean": [], "stack": []})
                record[source].append({"offset_arcsec": round(separation, 3), **r})
    for obj, record in matches.items():
        rows = query(
            "detection",
            obj,
            {"objID": obj},
            [
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
            ],
            fetch,
        )
        if any(r["objID"] != obj for r in rows) or len(
            {r["detectID"] for r in rows}
        ) != len(rows):
            raise ValueError("Wrong object or duplicate detection")
        groups = defaultdict(list)
        for r in rows:
            t, f, e, q = map(
                float, (r["obsTime"], r["psfFlux"], r["psfFluxErr"], r["psfQfPerfect"])
            )
            if (
                all(map(math.isfinite, (t, f, e, q)))
                and e > 0
                and f / e >= 5
                and q >= 0.85
                and sep(r["ra"], r["dec"]) <= 1.2
            ):
                groups[
                    (BANDS[int(r["filterID"]) - 1], str(pilot.day_of(t + 2400000.5)))
                ].append((-2.5 * math.log10(f / 3631), 2.5 * e / (math.log(10) * f)))
        record["raw_detections"] = len(rows)
        if record["stack"]:
            record["forced"] = forced_warps(obj, fetch)
        record["quality_nights"] = {}
        for b in BANDS:
            nights = sorted(
                (
                    round(median(v[0] for v in vs), 4),
                    round(median(v[1] for v in vs), 4),
                    day,
                )
                for (band, day), vs in groups.items()
                if band == b
            )
            if not nights:
                continue
            lo, hi = nights[0], nights[-1]
            gap = hi[0] - lo[0]
            record["quality_nights"][b] = {
                "count": len(nights),
                "brightest": list(lo),
                "faintest": list(hi),
                "selected": len(nights) >= 3
                and gap >= 1
                and gap >= 5 * math.hypot(lo[1], hi[1]),
            }
    return {
        "id": NAME,
        "mean_rows_10arcsec": len(mean),
        "stack_rows_10arcsec": len(stack),
        "matches_within_1.5arcsec": matches,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    result = run(args.fetch)
    path = HERE / "results.json"
    if args.fetch:
        pilot.json_write(path, result)
    elif json.loads(path.read_text()) != result:
        raise ValueError("Pan-STARRS results changed")
    print(json.dumps(result, indent=2))
