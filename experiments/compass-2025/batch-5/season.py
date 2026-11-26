"""Check the later seasons of SN 2025xza on original ZTF images."""

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import math
from statistics import stdev, mean
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "batch-4"))
import late
import pilot

POS = (39.4938772, -1.8497861)
FOLDER = HERE / "season-sources"
WINDOWS = {
    "before": ("2023-09-01", "2025-08-31", 12),
    "dec": ("2025-12-01", "2025-12-31", 6),
    "winter": ("2026-01-13", "2026-02-03", 6),
    "post": ("2026-02-04", "2026-03-31", 8),
    "apr_may": ("2026-04-01", "2026-05-31", 6),
    "jun_jul": ("2026-06-01", "2026-07-31", 6),
    "aug_sep": ("2026-08-01", "2026-09-30", 6),
    "oct_nov": ("2026-10-01", "2026-11-26", 6),
}


I_WINDOWS = {
    "before": WINDOWS["before"],
    "dec": WINDOWS["dec"],
    "march": ("2026-03-01", "2026-03-31", 31),
    "summer": ("2026-07-01", "2026-08-31", 62),
}


def select(inventory, band, windows):
    groups = defaultdict(dict)
    for r in inventory:
        if (
            r["field"] != "401"
            or r["filtercode"] != "z" + band
            or r["imgtypecode"] != "o"
            or int(r["infobits"]) >= 33554432
            or not r["maglimit"]
            or not r["seeing"]
            or float(r["maglimit"]) < (19 if band == "i" else 20)
            or float(r["seeing"]) > 3.2
        ):
            continue
        if 4 * (int(r["ccdid"]) - 1) + int(r["qid"]) - 1 != 34:
            continue
        d = str(pilot.day_of(float(r["obsjd"])))
        for key, (start, stop, _) in windows.items():
            if start <= d <= stop:
                order = lambda row: (
                    -float(row["maglimit"]),
                    float(row["obsjd"]),
                    int(row["pid"]),
                )
                if d not in groups[key] or order(r) < order(groups[key][d]):
                    groups[key][d] = r
    sample = []
    for key, (_, _, cap) in windows.items():
        dates = sorted(groups[key])
        n = min(cap, len(dates))
        inds = (
            [round(j * (len(dates) - 1) / (n - 1)) for j in range(n)]
            if n > 1
            else range(n)
        )
        sample.extend((key, dates[i], groups[key][dates[i]]) for i in inds)
    return {k: len(groups[k]) for k in windows}, sample


def nightly_contrasts(pre, post):
    checks = []
    for row in post:
        base = [r for r in pre if abs(r["seeing"] - row["seeing"]) <= 0.4]
        if len(base) < 4:
            continue
        differences = {}
        for key in ["source", *row["blank_flux_ujy"]]:

            def flux(r):
                return r["flux_ujy"] if key == "source" else r["blank_flux_ujy"][key]

            values = [
                flux(r) for r in base if key == "source" or key in r["blank_flux_ujy"]
            ]
            if len(values) < 4:
                continue
            noise = stdev(values)
            error = math.hypot(
                max(row["error_ujy"], noise), noise / math.sqrt(len(values))
            )
            differences[key] = dict(excess=flux(row) - mean(values), error=error)
        checks.append(
            dict(date=row["date"], n_baseline=len(base), differences=differences)
        )
    return checks


def run(fetch, band):
    folder = FOLDER if band == "r" else HERE / "season-i-sources"
    windows = WINDOWS if band == "r" else I_WINDOWS
    inventory, stars = late.tables(fetch, POS, folder, {"r": 2, "i": 3}[band])
    counts, sample = select(inventory, band, windows)
    offsets = [(a, b) for r in [12, 24] for a, b in [(0, r), (0, -r), (r, 0), (-r, 0)]]
    controls = {
        f"{a},{b}": (
            POS[0] + a / 3600 / math.cos(math.radians(POS[1])),
            POS[1] + b / 3600,
        )
        for a, b in offsets
    }
    with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
        rows = list(
            pool.map(
                lambda item: late.measure(item, stars, controls, fetch, POS, folder),
                sample,
            )
        )
    valid = [r for r in rows if r["status"] == "measured"]
    pre = [r for r in valid if r["window"] == "before"]
    comparisons = {}
    for key in windows:
        if key == "before":
            continue
        post = [r for r in valid if r["window"] == key]

        def common(a, b):
            return [
                r for r in a if any(abs(r["seeing"] - s["seeing"]) <= 0.4 for s in b)
            ]

        comparisons[key] = dict(
            nightly=nightly_contrasts(pre, post),
            all=late.contrast(pre, post),
            matched=late.contrast(common(pre, post), common(post, pre)),
        )
    return dict(
        inventory=len(inventory),
        eligible=counts,
        controls=controls,
        measurements=rows,
        comparisons=comparisons,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--band", choices=["r", "i"], default="r")
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    pilot.json_write(
        a.out or HERE / f"results/season-{a.band}.json", run(a.fetch, a.band)
    )
