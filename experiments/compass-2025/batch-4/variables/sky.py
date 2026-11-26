"""Fetch and compare same-window ASAS-SN forced flux for the two source positions."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json
import math
from pathlib import Path
from statistics import median

import numpy as np
from scipy.stats import spearmanr

from measure import HERE, NAMES, load, baseline, pilot
from asas import obtain

SRC = HERE / "sources"


def compare(name, fetch):
    case, pre, _ = load(name)
    SRC.mkdir(exist_ok=True)
    cone, raw = obtain(
        fetch, None, ra=case["ra"], dec=case["dec"], stem=name + "-asas", sources=SRC
    )
    if not raw:
        return {"cone": cone, "status": "No catalog photometry"}
    boundary = json.loads((HERE / "results.json").read_text())[name]["state"][
        "boundary_mjd"
    ]
    start = baseline.EPOCH + timedelta(days=boundary)
    anchor = date.fromisoformat(case["date"])
    daily = defaultdict(list)
    for v in raw:
        if (
            v["quality"] != "G"
            or not all(math.isfinite(v[k]) for k in ("jd", "flux", "flux_err"))
            or v["flux_err"] <= 0
        ):
            continue
        d = pilot.day_of(v["jd"])
        if date(2019, 1, 1) <= d <= anchor:
            daily[d].append(v)
    flux = {d: median(v["flux"] for v in vv) for d, vv in daily.items()}
    series = {d: (pilot.utc_jd(d), v) for d, v in flux.items()}
    overlap = sorted(pre["g"].keys() & flux.keys())

    def rank(dd):
        return (
            round(
                float(
                    spearmanr(
                        [flux[d] for d in dd], [-pre["g"][d][1] for d in dd]
                    ).statistic
                ),
                5,
            )
            if len(dd) >= 3
            else None
        )

    bins = baseline.bin_curve(series, anchor)
    contrast = baseline.change(bins, boundary)
    out = {
        "cone": cone,
        "archive_measurements": len(raw),
        "usable_dates": len(flux),
        "common_g_dates": len(overlap),
        "common_g_correlation": rank(overlap),
        "boundary": str(start),
        "bin_median_contrast_mjy": {
            k.replace("mag", "mjy"): v for k, v in contrast.items()
        },
        "windows": {},
    }
    for label, lo, hi in [("before", -180, 0), ("after", 0, 180)]:
        dates = sorted(d for d in flux if lo <= (d - start).days < hi)
        common = sorted(set(dates) & set(overlap))
        out["windows"][label] = {
            "dates": len(dates),
            "common_g_dates": len(common),
            "common_g_correlation": rank(common),
            "q20_q50_q80_mjy": (
                list(
                    np.round(np.quantile([flux[d] for d in dates], [0.2, 0.5, 0.8]), 6)
                )
                if dates
                else None
            ),
            "cameras": sorted({v["camera"] for d in dates for v in daily[d]}),
            "nights": [[str(d), round(flux[d], 6)] for d in dates],
        }
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "sky.json")
    a = p.parse_args()
    pilot.json_write(a.out, {n: compare(n, a.fetch) for n in NAMES})
