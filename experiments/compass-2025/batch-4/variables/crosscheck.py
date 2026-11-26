"""Compare calibrated forced flux, individual cameras and registered eruption lags."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json, math
from pathlib import Path
from statistics import median

import numpy as np
from scipy.stats import spearmanr

from measure import HERE, NAMES, baseline, load, pilot
from sky import SRC
from asas import obtain
from recurrence import evaluate


def good(r, start, stop):
    d = pilot.day_of(r["jd"])
    return (
        start <= d <= stop
        and r["quality"] == "G"
        and all(math.isfinite(r[k]) for k in ("flux", "flux_err"))
        and r["flux_err"] > 0
    )


def summarize(series, edge, anchor, green):
    grouped = defaultdict(list)
    for r in series:
        grouped[pilot.day_of(r["jd"])].append(r)
    gr = {
        d: (median(r["jd"] for r in rr), median(r["flux"] for r in rr))
        for d, rr in grouped.items()
    }
    change = baseline.change(baseline.bin_curve(gr, anchor), edge)
    overlap = sorted(gr.keys() & green.keys())
    rho = (
        float(
            spearmanr(
                [gr[d][1] for d in overlap], [-green[d][1] for d in overlap]
            ).statistic
        )
        if len(overlap) >= 3
        else None
    )
    return {
        "change_mjy": {k.replace("mag", "mjy"): v for k, v in change.items()},
        "common_g_dates": len(overlap),
        "common_g_rho": round(rho, 5) if rho is not None else None,
    }


def compact(rows):
    return {
        y: {k: v for k, v in r.items() if k != "per_anchor"}
        for y, r in evaluate(rows).items()
    }


def run(name, fetch):
    case, pre, _ = load(name)
    anchor = date.fromisoformat(case["date"])
    boundary = json.loads((HERE / "results.json").read_text())[name]["state"][
        "boundary_mjd"
    ]
    edge = baseline.EPOCH + timedelta(days=boundary)
    cone, rows = obtain(
        fetch,
        None,
        neighbors=True,
        ra=case["ra"],
        dec=case["dec"],
        sources=SRC,
        stem=name + "-asas-field",
    )
    if not rows:
        return {"coverage": "empty"}
    target = min(cone, key=lambda x: x["separation_arcsec"])["asas_sn_id"]
    groups = defaultdict(list)
    for r in rows:
        if good(r, date(2019, 1, 1), anchor):
            groups[r["asas_sn_id"]].append(r)

    def outside(records):
        return [
            r for r in records if not -180 <= (pilot.day_of(r["jd"]) - edge).days < 180
        ]

    typical = median(r["flux"] for r in outside(groups[target]))
    baselines = {}
    info = []
    for star in sorted(cone, key=lambda v: v["separation_arcsec"]):
        oid = star["asas_sn_id"]
        training = outside(groups[oid])
        if star["separation_arcsec"] < 30 or not training:
            continue
        level = median(r["flux"] for r in training)
        if not 0.5 * typical <= level <= 2 * typical:
            continue
        cams = {
            c: median(r["flux"] for r in training if r["camera"] == c)
            for c in set(r["camera"] for r in training)
        }
        baselines[oid] = cams
        info.append(
            {
                "id": oid,
                "separation_arcsec": star["separation_arcsec"],
                "outside_median_mjy": round(level, 6),
                "catalog_memberships": star["catalog_sources"],
            }
        )
        if len(info) == 3:
            break
    refs = defaultdict(dict)
    for oid, cams in baselines.items():
        for r in groups[oid]:
            if r["camera"] in cams:
                refs[(r["camera"], r["image_id"])][oid] = r["flux"] - cams[r["camera"]]
    common = []
    corrected = []
    for r in groups[target]:
        vals = refs[(r["camera"], r["image_id"])]
        if len(vals) >= 2:
            common.append(r)
            corrected.append({**r, "flux": r["flux"] - median(vals.values())})
    cameras = sorted(set(r["camera"] for r in groups[target]))
    out = {
        "target_asassn_id": target,
        "field_cone_sources": len(cone),
        "comparison_stars": info,
        "matched_target": summarize(common, boundary, anchor, pre["g"]),
        "matched_corrected": summarize(corrected, boundary, anchor, pre["g"]),
        "by_camera": {
            c: summarize(
                [r for r in groups[target] if r["camera"] == c],
                boundary,
                anchor,
                pre["g"],
            )
            for c in cameras
        },
        "lag_windows_days": [[6, 10], [14, 20], [24, 29]],
        "lag_scores": {"asas_forced": compact(groups[target])},
    }
    for b in "rg":
        nights = [
            {"jd": pilot.utc_jd(d), "flux": 3631000 * 10 ** (-0.4 * mag)}
            for d, (t, mag) in pre[b].items()
        ]
        out["lag_scores"][b] = compact(nights)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "crosscheck.json")
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    pilot.json_write(args.out, {n: run(n, args.fetch) for n in NAMES})
