"""Verify how same-exposure field correction affects the candidate's held-out clock."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

import numpy as np

import check
from check import HERE, BATCH, NAME, compass_direct, pilot


def separation(a, b):
    return math.hypot(
        *pilot.offsets(float(a["ra"]), float(a["dec"]), float(b["ra"]), float(b["dec"]))
    )


def select(case, selected, fetch, per_band):
    red, green = selected[2], selected[1]
    if any(red[k] != green[k] for k in ("field", "ccdid", "qid")):
        raise ValueError("Need different detector queries for the filters")
    columns = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag"
    query = (
        f"SELECT TOP 10000 {columns} FROM ztf_objects_dr24 WHERE "
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{case['ra']},{case['dec']},{90/3600}))=1 "
        f"AND field={red['field']} AND ccdid={red['ccdid']} AND qid={red['qid']} "
        "AND fid IN (1,2) AND nobs>=40"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
        )
    )
    cone = compass_direct.capture(
        HERE / "sources", "field", url, set(columns.split(",")), fetch
    )
    by_id = {r["oid"]: r for r in cone}
    if len(cone) >= 10000 or len(cone) != len(by_id):
        raise ValueError("Truncated or duplicate DR24 source cone")
    cand = {}
    for band, target in [("r", red), ("g", green)]:
        typical = float(by_id[target["oid"]]["meanmag"])
        cand[band] = [
            r
            for r in cone
            if r["fid"] == target["fid"]
            and r["meanmag"]
            and abs(float(r["meanmag"]) - typical) <= 1
            and 5 <= separation(r, case) <= 90.001
        ]
    if per_band:
        ranked = {
            b: sorted(rows, key=lambda r: (separation(r, case), int(r["oid"])))[:3]
            for b, rows in cand.items()
        }
        return [
            {
                "red": r["oid"],
                "green": g["oid"],
                "r_mag": float(r["meanmag"]),
                "g_mag": float(g["meanmag"]),
                "red_distance": round(separation(r, case), 4),
                "green_distance": round(separation(g, case), 4),
            }
            for r, g in zip(ranked["r"], ranked["g"])
        ]
    used = set()
    pairs = []
    for r in sorted(cand["r"], key=lambda x: (separation(x, case), int(x["oid"]))):
        choices = sorted(
            (g for g in cand["g"] if g["oid"] not in used and separation(r, g) <= 1),
            key=lambda g: (separation(r, g), int(g["oid"])),
        )
        if choices:
            g = choices[0]
            used.add(g["oid"])
            pairs.append(
                {
                    "red": r["oid"],
                    "green": g["oid"],
                    "r_mag": float(r["meanmag"]),
                    "g_mag": float(g["meanmag"]),
                    "distance": round(separation(r, case), 4),
                    "rg_separation": round(separation(r, g), 4),
                }
            )
        if len(pairs) == 3:
            break
    return pairs


def analyze(fetch, per_band):
    case = next(
        r
        for r in json.loads((BATCH / "results/ledger.json").read_text())["objects"]
        if r["id"] == NAME
    )
    selected, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    pairs = select(case, selected, fetch, per_band)
    out = {"controls": pairs, "bands": {}}
    if len(pairs) < 2:
        out["status"] = "Fewer than two cross-matched stars passed comparison cuts"
        return out, {}
    ids = [p[k] for p in pairs for k in ("red", "green")]
    params = [("ID", oid) for oid in ids] + [
        ("COLLECTION", "ztf_dr24"),
        ("FORMAT", "csv"),
        ("BAD_CATFLAGS_MASK", "65535"),
    ]
    controls = compass_direct.capture(
        HERE / "sources",
        "per-band-comparison" if per_band else "comparison",
        compass_direct.SERVICE + "?" + urlencode(params),
        {"oid", "expid", "mjd", "mag", "magerr", "catflags", "filtercode"},
        fetch,
    )
    if set(ids) != {r["oid"] for r in controls}:
        raise ValueError("Comparison photometry does not match selected sources")
    timing = json.loads((HERE / "results.json").read_text())["timing"]
    cutoff = check.date.fromisoformat(timing["training_last"])
    f = 1 / timing["models"]["training_best"]["period_days"]
    pre = check.nightly(rows, case["jd"] + 1e-7)
    t, y = np.array([v for d, v in pre["r"].items() if d <= cutoff]).T
    origin = t[0]
    _, coeff = check.fit(t - origin, y, f)
    centers = {}
    for oid in ids:
        vals = [
            float(r["mag"])
            for r in controls
            if r["oid"] == oid
            and pilot.day_of(float(r["mjd"]) + 2400000.5) <= cutoff
            and int(r["catflags"]) == 0
            and 0 < float(r["magerr"]) <= 0.2
            and math.isfinite(float(r["mag"]))
        ]
        if vals:
            centers[oid] = median(vals)
    by_exposure = defaultdict(dict)
    for r in controls:
        if (
            r["oid"] in centers
            and int(r["catflags"]) == 0
            and 0 < float(r["magerr"]) <= 0.2
            and math.isfinite(float(r["mag"]))
        ):
            by_exposure[(r["expid"], r["filtercode"])][r["oid"]] = (
                float(r["mag"]) - centers[r["oid"]]
            )
    corrected = []
    for r in rows:
        ensemble = by_exposure[(r["expid"], r["filtercode"])]
        if len(ensemble) >= 2:
            corrected.append({**r, "mag": float(r["mag"]) - median(ensemble.values())})
    curves = check.nightly(corrected, case["jd"] + 1e-7)
    for b in "rg":
        train = [v for d, v in pre[b].items() if d <= cutoff]
        tt, yy = np.array(train).T
        offset = np.mean(yy - check.design(tt - origin, f) @ coeff) if b == "g" else 0
        dates = [d for d in curves[b] if d > cutoff]
        out["bands"][b] = {
            "target": check.score(dates, pre[b], f, coeff, np.mean(yy), origin, offset),
            "corrected": check.score(
                dates, curves[b], f, coeff, np.mean(yy), origin, offset
            ),
        }
    return out, curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--per-band", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result, _ = analyze(args.fetch, args.per_band)
    name = "field-per-band.json" if args.per_band else "field-results.json"
    pilot.json_write(args.out or HERE / name, result)


if __name__ == "__main__":
    main()
