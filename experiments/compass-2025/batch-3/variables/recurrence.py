"""Test two-to-three-week recurrence on unused dates from ASAS-SN."""

from collections import defaultdict
from datetime import date
from itertools import combinations
import argparse
import json
from pathlib import Path
from statistics import median
import math

from asas import obtain
from analyze import HERE, curves, pilot
from interval import NAME, START, END

WINDOWS = [(6, 10), (14, 20), (24, 29)]


def excluded():
    intervals = json.loads((HERE / "eruptions.json").read_text())[
        "well_sampled_intervals"
    ]
    exclude = [(START, END)] + [
        (date.fromisoformat(e["start"]), date.fromisoformat(e["end"]))
        for e in intervals
    ]
    return exclude


def accepted(raw):
    exclude = excluded()
    stop = date.fromisoformat(curves(NAME)[0]["date"])
    keep = []
    for r in raw:
        d = pilot.day_of(r["jd"])
        if (
            r["quality"] == "G"
            and date(2019, 1, 1) <= d <= stop
            and all(math.isfinite(r[k]) for k in ("flux", "flux_err"))
            and r["flux_err"] > 0
            and all(not (first <= d <= last) for first, last in exclude)
        ):
            keep.append(r)
    return keep


def evaluate(raw):
    grouped = defaultdict(list)
    for r in raw:
        grouped[pilot.day_of(r["jd"])].append(r["flux"])
    flux = {d: median(vals) for d, vals in grouped.items()}
    out = {}
    for year in sorted({d.year for d in flux}):
        dd = sorted(d for d in flux if d.year == year)
        lags = defaultdict(lambda: [[], [], []])
        for a, b in combinations(dd, 2):
            lag = (b - a).days
            for j, (lo, hi) in enumerate(WINDOWS):
                if lo <= lag <= hi:
                    lags[a][j].append(abs(flux[a] - flux[b]))
        rows = []
        for d, values in sorted(lags.items()):
            if all(values):
                rows.append(
                    {
                        "date": str(d),
                        "lag_medians": [round(median(v), 6) for v in values],
                    }
                )
        med = (
            [round(median(r["lag_medians"][j] for r in rows), 6) for j in range(3)]
            if rows
            else None
        )
        out[str(year)] = {
            "available_nights": len(dd),
            "common_anchors": len(rows),
            "delta_flux_medians_mjy": med,
            "middle_smaller_than_both": bool(med and med[1] < min(med[0], med[2])),
            "per_anchor": rows,
        }
    return out


def run():
    _, raw = obtain(False, None)
    return {"lag_windows_days": WINDOWS, "years": evaluate(accepted(raw))}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "recurrence.json")
    a = p.parse_args()
    pilot.json_write(a.out, run())
