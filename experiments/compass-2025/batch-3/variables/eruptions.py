"""Count bracketed red brightenings and well-sampled intervals in MGAB-V243."""

import argparse
from datetime import date, timedelta
import json
from pathlib import Path
from statistics import median
import numpy as np
from analyze import HERE, curves, pilot


def run():
    name = "ZTF17aabshbq"
    case, data = curves(name)
    cutoff = date.fromisoformat(
        json.loads((HERE / "results.json").read_text())[name]["training_last"]
    )
    thresholds = {
        b: np.quantile(
            [y for d, (t, y) in data[b].items() if d <= cutoff], [0.25, 0.75]
        )
        for b in "rg"
    }
    red, green = data["r"], data["g"]
    red_bright, red_faint = thresholds["r"]
    green_bright = thresholds["g"][0]
    groups = []
    current = []
    for d, (_, m) in red.items():
        if m <= red_bright:
            current.append(d)
        if m >= red_faint and current:
            groups.append(current)
            current = []
    if current:
        groups.append(current)
    results = []
    for dd in groups:
        close = sorted(
            {
                d
                for i in range(1, len(dd))
                if (dd[i] - dd[i - 1]).days <= 10
                for d in (dd[i - 1], dd[i])
            }
        )
        g = [
            d
            for d, (_, mag) in green.items()
            if mag <= green_bright and any(abs((d - z).days) <= 2 for z in close)
        ]
        returns = [d for d, (_, m) in red.items() if m >= red_faint]
        before = [d for d in returns if d < dd[0]]
        after = [d for d in returns if d > dd[-1]]
        results.append(
            {
                "first_bright": str(dd[0]),
                "last_bright": str(dd[-1]),
                "bright_dates": len(dd),
                "resolved": bool(close),
                "green_matches": list(map(str, g)),
                "brightest_mag": min(red[d][1] for d in dd),
                "prior_faint": str(before[-1]) if before else None,
                "next_faint": str(after[0]) if after else None,
            }
        )
    intervals = []
    for a, b in zip(results, results[1:]):
        if a["resolved"] and b["resolved"]:
            first, last = map(
                date.fromisoformat, [a["first_bright"], b["first_bright"]]
            )
            dates = [d for d in red if first <= d <= last]
            gap = max((v - u).days for u, v in zip(dates, dates[1:]))
            if gap <= 5 and a["next_faint"] and a["next_faint"] < b["first_bright"]:
                intervals.append(
                    {
                        "start": str(first),
                        "end": str(last),
                        "days": (last - first).days,
                        "largest_observing_gap_days": gap,
                    }
                )
    return {
        "training_last": str(cutoff),
        "training_thresholds_mag": {
            b: {"bright_q25": float(q[0]), "faint_q75": float(q[1])}
            for b, q in thresholds.items()
        },
        "groups": results,
        "well_sampled_intervals": intervals,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "eruptions.json")
    args = p.parse_args()
    pilot.json_write(args.out, run())
