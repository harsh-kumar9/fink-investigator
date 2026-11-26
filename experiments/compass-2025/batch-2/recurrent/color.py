"""Match brightening and fading dates at the same direct red magnitude."""

import argparse
import json
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from statistics import median

import numpy as np
import field as field_check

from check import BATCH, HERE, NAME, compass_direct, design, fit, nightly, pilot


def match(curves, model_curves, period, cutoff, gap_days=0.25, mag_tol=0.1):
    red, green = curves["r"], curves["g"]
    t, y = np.array([v for d, v in model_curves["r"].items() if d <= cutoff]).T
    origin = t[0]
    frequency = 1 / period
    _, coeff = fit(t - origin, y, frequency)
    ph = np.linspace(0, 1, 10000, endpoint=False)
    peak = ph[np.argmin(design(ph * period, frequency) @ coeff)]
    groups = defaultdict(lambda: defaultdict(list))
    coverage = defaultdict(Counter)
    for d, (rt, r) in red.items():
        phase = (rt - origin) * frequency
        w = 2 * np.pi * phase
        slope = (2 * np.pi * frequency) * (
            coeff[1] * np.cos(w)
            - coeff[2] * np.sin(w)
            + 2 * coeff[3] * np.cos(2 * w)
            - 2 * coeff[4] * np.sin(2 * w)
        )
        if not 13 <= r <= 15.5 or abs(slope) < 0.01:
            continue
        branch = "fading" if slope > 0 else "brightening"
        coverage[branch]["red_eligible"] += 1
        if d not in green or abs(green[d][0] - rt) > gap_days:
            continue
        coverage[branch]["also_measured_green"] += 1
        if green[d][1] >= 19:
            continue
        coverage[branch]["brighter_than_green_limit"] += 1
        cycle = int(np.floor(phase - peak + 0.5))
        groups[cycle][branch].append((d, r, green[d][1] - r))
    pairs = []
    for cycle, branches in sorted(groups.items()):
        possible = sorted(
            (abs(f[1] - b[1]), b, f)
            for b in branches["brightening"]
            for f in branches["fading"]
            if abs(f[1] - b[1]) <= mag_tol
        )
        used = set()
        for dist, b, f in possible:
            if b[0] in used or f[0] in used:
                continue
            used.update((b[0], f[0]))
            subset = (
                "later"
                if min(b[0], f[0]) > cutoff
                else "training" if max(b[0], f[0]) <= cutoff else "straddles"
            )
            pairs.append(
                {
                    "cycle": cycle,
                    "subset": subset,
                    "brightening_date": str(b[0]),
                    "fading_date": str(f[0]),
                    "r_brightening": round(b[1], 5),
                    "r_fading": round(f[1], 5),
                    "color_brightening": round(b[2], 5),
                    "color_fading": round(f[2], 5),
                    "color_difference": round(f[2] - b[2], 5),
                }
            )
    per_cycle = []
    for cycle in sorted({p["cycle"] for p in pairs}):
        p = [r for r in pairs if r["cycle"] == cycle]
        per_cycle.append(
            {
                "cycle": cycle,
                "n": len(p),
                "median_color_difference": median(r["color_difference"] for r in p),
                "subsets": dict(Counter(r["subset"] for r in p)),
            }
        )
    return {
        "coverage": dict(coverage),
        "n": len(pairs),
        "cycles": per_cycle,
        "median": median(p["color_difference"] for p in pairs) if pairs else None,
        "pairs": pairs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "color.json")
    args = parser.parse_args()
    case = next(
        r
        for r in json.loads((BATCH / "results/ledger.json").read_text())["objects"]
        if r["id"] == NAME
    )
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    pre = nightly(rows, case["jd"] + 1e-7)
    timing = json.loads((HERE / "results.json").read_text())["timing"]
    cutoff = date.fromisoformat(timing["training_last"])
    best = timing["models"]["training_best"]["period_days"]
    out = {
        name: match(pre, pre, period, cutoff, gap, tolerance)
        for name, period, gap, tolerance in [
            ("primary", best, 0.25, 0.1),
            ("tight_r", best, 0.25, 0.05),
            ("catalog_period", 325.0, 0.25, 0.1),
            ("tight_g_time", best, 0.125, 0.1),
        ]
    }
    field = json.loads((HERE / "field-results.json").read_text())
    if field["controls"]:
        raise ValueError(
            "Add the registered field-corrected comparison before interpreting the color test"
        )
    out["field_corrected"] = field["status"]
    _, corrected = field_check.analyze(False, True)
    out["field_corrected_per_band"] = (
        match(corrected, pre, best, cutoff) if corrected else None
    )
    pilot.json_write(args.out, out)


if __name__ == "__main__":
    main()
