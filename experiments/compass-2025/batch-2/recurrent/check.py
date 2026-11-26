"""Replay the registered episode and held-out timing checks for the recurrent source."""

import argparse
from collections import defaultdict
from datetime import date
import json
import math
from pathlib import Path
from statistics import median
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
sys.path.insert(0, str(BATCH.parents[2]))
import compass_direct
import pilot

NAME = "ZTF18abmfrjj"


def nightly(rows, stop):
    grouped = defaultdict(dict)
    for row in rows:
        if int(row["catflags"]) or not 0 < float(row["magerr"]) <= 0.2:
            continue
        mjd, hjd, mag = map(float, (row["mjd"], row["hjd"], row["mag"]))
        if not all(map(math.isfinite, (mjd, hjd, mag))) or mjd + 2400000.5 > stop:
            continue
        key = row["filtercode"][-1], pilot.day_of(mjd + 2400000.5)
        expid = row["expid"]
        if expid in grouped[key] and grouped[key][expid] != (hjd, mag):
            raise ValueError("Conflicting records for the same exposure")
        grouped[key][expid] = hjd, mag
    data = {b: {} for b in "gr"}
    for (band, d), exposures in sorted(grouped.items()):
        t, y = zip(*exposures.values())
        data[band][d] = median(t), median(y)
    return data


def episodes(data, typical, anchor):
    red, green = data["r"], data["g"]
    returns = [d for d, (_, mag) in red.items() if mag >= typical["r"] - 0.5]
    bright_green = [d for d, (_, mag) in green.items() if mag <= typical["g"] - 2]
    groups = []
    current = []
    for day, (_, mag) in red.items():
        if mag <= typical["r"] - 2:
            current.append(day)
        elif mag >= typical["r"] - 0.5 and current:
            groups.append(current)
            current = []
    if current:
        groups.append(current)

    output = []
    for dates in groups:
        pre = [d for d in returns if d < dates[0]]
        post = [d for d in returns if d > dates[-1]]
        matches = [
            d for d in bright_green if any(abs((d - r).days) <= 2 for r in dates)
        ]
        output.append(
            {
                "bright_red_dates": list(map(str, dates)),
                "brightest_r": min(red[d][1] for d in dates),
                "resolved": any((b - a).days <= 10 for a, b in zip(dates, dates[1:])),
                "green_matches": list(map(str, matches)),
                "previous_return": str(pre[-1]) if pre else None,
                "next_return": str(post[0]) if post else None,
                "sample": (
                    "pre_anchor"
                    if dates[-1] <= anchor
                    else "crosses_anchor" if dates[0] <= anchor else "post_anchor"
                ),
            }
        )
    return output


def design(t, frequency):
    angle = 2 * np.pi * frequency * t
    return np.column_stack(
        [
            np.ones(len(t)),
            np.sin(angle),
            np.cos(angle),
            np.sin(2 * angle),
            np.cos(2 * angle),
        ]
    )


def fit(t, mag, frequency):
    matrix = design(t, frequency)
    coefficients = np.linalg.lstsq(matrix, mag, rcond=None)[0]
    return np.sum((mag - matrix @ coefficients) ** 2), coefficients


def score(dates, curve, frequency, coefficients, baseline, origin, offset=0):
    if len(dates) < 2:
        return {"n": len(dates)}
    t, y = np.asarray([curve[d] for d in dates]).T
    predicted = design(t - origin, frequency) @ coefficients + offset
    gain = 1 - np.sum((y - predicted) ** 2) / np.sum((y - baseline) ** 2)
    rho = float(spearmanr(y, predicted).statistic)
    return {
        "n": len(dates),
        "start": str(dates[0]),
        "end": str(dates[-1]),
        "error_reduction": round(float(gain), 6),
        "spearman": round(rho, 6),
        "passes": bool(len(dates) >= 30 and gain >= 0.2 and rho >= 0.5),
    }


def timing(pre, all_dates):
    red, green = pre["r"], pre["g"]
    n = len(red) * 2 // 3
    if n < 5 or len(red) - n < 30:
        return {"status": "insufficient red dates"}
    train, test = list(red)[:n], list(red)[n:]
    t, y = np.asarray([red[d] for d in train]).T
    origin = t[0]
    t = t - origin
    span = t[-1]
    frequencies = np.linspace(
        1 / 1000, 1 / 2, math.ceil((1 / 2 - 1 / 1000) * 5 * span) + 1
    )
    errors = [fit(t, y, f)[0] for f in frequencies]
    index = int(np.argmin(errors))
    lo, hi = (
        frequencies[max(index - 1, 0)],
        frequencies[min(index + 1, len(errors) - 1)],
    )
    opt = minimize_scalar(
        lambda f: fit(t, y, f)[0],
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": 1e-13},
    )
    best = min((lo, frequencies[index], opt.x, hi), key=lambda f: fit(t, y, f)[0])
    models = {
        "training_best": best,
        "half_period": best * 2,
        "double_period": best / 2,
        "published_325": 1 / 325,
        "published_328.22": 1 / 328.22,
        "atlas_fourier": 1 / 0.993990,
        "atlas_alternate": 1 / 443.379574,
    }
    for label, f in (("annual", 1 / 365.25), ("sidereal_day", 1 / 0.99726957)):
        models[label + "_plus"] = best + f
        if best != f:
            models[label + "_minus"] = abs(best - f)
    gt = [d for d in green if d <= train[-1]]
    ge = [d for d in green if d > train[-1]]
    out = {}
    for name, f in models.items():
        sse, coeff = fit(t, y, f)
        item = {
            "period_days": round(float(1 / f), 8),
            "training_error": round(float(sse), 6),
            "red": score(test, red, f, coeff, np.mean(y), origin),
        }
        later = [d for d in all_dates["r"] if d > max(red)]
        item["red_after_anchor"] = score(
            later, all_dates["r"], f, coeff, np.mean(y), origin
        )
        if gt:
            time, mag = np.asarray([green[d] for d in gt]).T
            offset = np.mean(mag - design(time - origin, f) @ coeff)
            item["green"] = score(ge, green, f, coeff, np.mean(mag), origin, offset)
            item["green_unshared"] = score(
                [d for d in ge if d not in all_dates["r"]],
                green,
                f,
                coeff,
                np.mean(mag),
                origin,
                offset,
            )
            later = [d for d in all_dates["g"] if d > max(red) and d not in green]
            item["green_after_anchor"] = score(
                later, all_dates["g"], f, coeff, np.mean(mag), origin, offset
            )
        out[name] = item
    return {
        "training_dates": n,
        "training_first": str(train[0]),
        "training_last": str(train[-1]),
        "models": out,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    ledger = json.loads((BATCH / "results/ledger.json").read_text())
    case = next(r for r in ledger["objects"] if r["id"] == NAME)
    assert case["selected_for"] == "recurrent"
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    pre, full = nightly(rows, case["jd"] + 1e-7), nightly(rows, math.inf)
    anchor = date.fromisoformat(case["date"])
    typical = {b: median(x[1] for x in pre[b].values()) for b in "gr"}
    bands = {}
    for b in "gr":
        y = [x[1] for x in pre[b].values()]
        dates = list(full[b])
        bands[b] = {
            "pre_dates": len(y),
            "total_dates": len(dates),
            "first": str(dates[0]),
            "last": str(dates[-1]),
            "typical": typical[b],
            "min_mag_through_anchor": min(y),
            "max_mag_through_anchor": max(y),
            "p05_to_p95_mag": list(np.round(np.quantile(y, [0.05, 0.95]), 5)),
            "maximum_observing_gap_days": max(
                (v - u).days for u, v in zip(dates, dates[1:])
            ),
        }
    result = {
        "source": NAME,
        "anchor": str(anchor),
        "bands": bands,
        "pre_anchor_episodes": episodes(pre, typical, anchor),
        "full_series_episodes": episodes(full, typical, anchor),
        "timing": timing(pre, full),
    }
    pilot.json_write(args.out, result)


if __name__ == "__main__":
    main()
