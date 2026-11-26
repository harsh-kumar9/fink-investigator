"""Replay timing, state and equal-brightness color checks on the two long histories."""

import argparse
from collections import defaultdict
from datetime import date
import json
import math
from pathlib import Path
from statistics import median
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
sys.path.insert(0, str(BATCH.parent / "batch-3/variables"))
import analyze as baseline
from analyze import previous, compass_direct, pilot

NAMES = ["ZTF17aabutiy", "ZTF18aaadshm"]


def load(name):
    ledger = json.loads((BATCH / "results/ledger.json").read_text())
    case = next(r for r in ledger["objects"] if r["id"] == name)
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / name)
    return (
        case,
        previous.nightly(rows, case["jd"] + 1e-7),
        previous.nightly(rows, math.inf),
    )


def match_color(pre, cutoff, tolerance, delta_t):
    red, green = pre["r"], pre["g"]
    dd = list(red)
    bins = defaultdict(lambda: defaultdict(list))
    counts = defaultdict(lambda: defaultdict(int))
    for i, d in enumerate(dd):
        if i == 0 or i + 1 == len(dd):
            continue
        before, after = dd[i - 1], dd[i + 1]
        if max((after - d).days, (d - before).days) > 15:
            continue
        y0, y1 = red[before][1], red[after][1]
        rt, mag = red[d]
        if not (
            13 <= mag <= 19.5
            and abs(y1 - y0) >= 0.25
            and min(y0, y1) - 0.08 <= mag <= max(y0, y1) + 0.08
        ):
            continue
        branch = "rising" if y1 < y0 else "fading"
        counts[branch]["red_eligible"] += 1
        if d not in green or abs(green[d][0] - rt) > delta_t:
            continue
        counts[branch]["paired_color"] += 1
        if green[d][1] >= 20:
            continue
        counts[branch]["retained"] += 1
        subset = "early" if d <= cutoff else "later"
        bins[(d.year, subset)][branch].append((d, mag, green[d][1] - mag))
    out = []
    for (year, subset), branches in sorted(bins.items()):
        possible = sorted(
            (abs(f[1] - r[1]), r, f)
            for r in branches["rising"]
            for f in branches["fading"]
            if abs(f[1] - r[1]) <= tolerance
        )
        used = set()
        pairs = []
        for delta, r, f in possible:
            if r[0] in used or f[0] in used:
                continue
            used.update([r[0], f[0]])
            pairs.append(
                {
                    "rising": str(r[0]),
                    "fading": str(f[0]),
                    "r_mag_rising": round(r[1], 5),
                    "r_mag_fading": round(f[1], 5),
                    "rising_g_minus_r": round(r[2], 5),
                    "fading_g_minus_r": round(f[2], 5),
                    "difference_mag": round(f[2] - r[2], 5),
                }
            )
        out.append(
            {
                "year": year,
                "subset": subset,
                "pairs": pairs,
                "n": len(pairs),
                "median": median(p["difference_mag"] for p in pairs) if pairs else None,
            }
        )
    return {"coverage": dict(counts), "years": out}


def run():
    out = {}
    for name in NAMES:
        case, pre, full = load(name)
        anchor = date.fromisoformat(case["date"])
        cutoff, scanned, scores, waves = baseline.get_models(pre, name)
        limit = date.fromisoformat(cutoff)
        for key, wave in waves.items():
            for b in "rg":
                later = [d for d in full[b] if d > anchor]
                ref = np.mean([m for d, (_, m) in pre[b].items() if d <= limit])
                origin, f, c = wave[b]
                scores[key][b + "_post_anchor"] = previous.score(
                    later, full[b], f, c, ref, origin
                )
        raw = baseline.state(pre, anchor)
        residuals = {}
        for b in "rg":
            origin, f, c = waves["early_red"][b]
            dd = list(pre[b])
            t, y = baseline.arrays(pre[b], dd)
            residuals[b] = dict(zip(dd, zip(t, y - previous.design(t - origin, f) @ c)))
        item = {
            "anchor": str(anchor),
            "training_through": cutoff,
            "frequencies": scanned,
            "bands": {},
            "clocks": scores,
            "state": raw,
            "color": {
                k: match_color(pre, limit, tol, dt)
                for k, tol, dt in [
                    ("registered", 0.08, 0.12),
                    ("tight_mag", 0.04, 0.12),
                    ("tight_time", 0.08, 0.06),
                ]
            },
        }
        if "boundary_mjd" in raw:
            item["state_after_clock"] = baseline.state(
                residuals, anchor, raw["boundary_mjd"]
            )
        for b in "rg":
            y = [m for _, m in pre[b].values()]
            item["bands"][b] = {
                "first": str(min(full[b])),
                "last": str(max(full[b])),
                "pre_anchor": len(y),
                "p05_p50_p95": list(np.round(np.quantile(y, [0.05, 0.5, 0.95]), 5)),
            }
        out[name] = item
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "results.json")
    args = p.parse_args()
    pilot.json_write(args.out, run())
