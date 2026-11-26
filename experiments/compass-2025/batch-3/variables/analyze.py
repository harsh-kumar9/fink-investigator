"""Search for predictable cycles and sustained changes in the two batch-three variables."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json
import math
from pathlib import Path
from statistics import median
import sys

from astropy.timeseries import LombScargle
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
sys.path.insert(0, str(BATCH.parent / "batch-2/recurrent"))
import check as previous
from check import compass_direct, pilot

EPOCH = date(1858, 11, 17)
NAMES = ["ZTF17aabshbq", "ZTF18aazutgw"]
CATALOG_PERIODS = {
    "ZTF18aazutgw": {
        "vsx": 385.1017206,
        "year": 365.25,
        "chen_green": 194.0059046,
        "chen_red": 188.3239431,
        "gaia_dr3": 380.66803631478734,
    }
}


def arrays(curve, dates):
    t, y = np.array([curve[d] for d in dates]).T
    return t, y


def search(t, y):
    fmin, fmax = 1 / 1000, 1 / 0.1
    frequency = np.linspace(
        fmin, fmax, math.ceil((fmax - fmin) * 10 * (t[-1] - t[0])) + 1
    )
    model = LombScargle(t, y, nterms=2)
    power = model.power(frequency, method="fastchi2", assume_regular_frequency=True)
    best = int(np.nanargmax(power))
    lo, hi = frequency[max(0, best - 1)], frequency[min(best + 1, len(frequency) - 1)]
    exact = lambda f: previous.fit(t, y, f)[0]
    opt = minimize_scalar(
        exact, method="bounded", bounds=(lo, hi), options={"xatol": 1e-13}
    )
    return min([lo, frequency[best], opt.x, hi], key=exact), len(frequency)


def get_models(pre, name):
    train = list(pre["r"])[: len(pre["r"]) * 2 // 3]
    cutoff = train[-1]
    t, y = arrays(pre["r"], train)
    origin = t[0]
    best, scanned = search(t - origin, y)
    frequencies = {
        "early_red": best,
        "half_period": 2 * best,
        "double_period": best / 2,
    }
    for label, df in [("annual", 1 / 365.25), ("sidereal_day", 1 / 0.99726957)]:
        frequencies[label + "_plus"] = best + df
        frequencies[label + "_minus"] = abs(best - df)
    frequencies.update({k: 1 / p for k, p in CATALOG_PERIODS.get(name, {}).items()})
    scores = {}
    waves = {}
    for label, f in frequencies.items():
        _, c = previous.fit(t - origin, y, f)
        scores[label] = {"period_days": round(float(1 / f), 9)}
        waves[label] = {}
        for band in "rg":
            tr = [d for d in pre[band] if d <= cutoff]
            te = [d for d in pre[band] if d > cutoff]
            tt, yy = arrays(pre[band], tr)
            predicted = previous.design(tt - origin, f) @ c
            if band == "g":
                scale, offset = np.linalg.lstsq(
                    np.column_stack([predicted, np.ones(len(tr))]), yy, rcond=None
                )[0]
            else:
                offset, scale = 0.0, 1.0
            coefficients = c * scale
            coefficients[0] += offset
            waves[label][band] = (origin, f, coefficients)
            scores[label][band] = previous.score(
                te, pre[band], f, coefficients, np.mean(yy), origin
            )
            if band == "g":
                scores[label]["green_without_red"] = previous.score(
                    [d for d in te if d not in pre["r"]],
                    pre[band],
                    f,
                    coefficients,
                    np.mean(yy),
                    origin,
                )
    return str(cutoff), scanned, scores, waves


def bin_curve(curve, anchor):
    groups = defaultdict(list)
    for d, (_, m) in curve.items():
        day = (d - EPOCH).days
        start = 58000 + 30 * ((day - 58000) // 30)
        if EPOCH + timedelta(days=start + 30) <= anchor + timedelta(days=1):
            groups[start].append(m)
    return {
        start: float(median(values))
        for start, values in groups.items()
        if len(values) >= 3
    }


def change(bins, edge):
    before = [(d, v) for d, v in sorted(bins.items()) if edge - 180 <= d < edge]
    after = [(d, v) for d, v in sorted(bins.items()) if edge <= d < edge + 180]
    if min(len(before), len(after)) < 3:
        return {
            "before_bins": len(before),
            "after_bins": len(after),
            "after_minus_before_mag": None,
        }
    v0, v1 = [[v for _, v in sub] for sub in (before, after)]
    q = median(v1) - median(v0)
    leave = [median(v1) - median(v0[:i] + v0[i + 1 :]) for i in range(len(v0))]
    leave.extend(median(v1[:i] + v1[i + 1 :]) - median(v0) for i in range(len(v1)))
    return {
        "earlier_bin_starts": [str(EPOCH + timedelta(days=d)) for d, _ in before],
        "later_bin_starts": [str(EPOCH + timedelta(days=d)) for d, _ in after],
        "before_bins": len(before),
        "after_bins": len(after),
        "after_minus_before_mag": round(q, 5),
        "leave_one_bin_out_range": [round(min(leave), 5), round(max(leave), 5)],
    }


def state(curves, anchor, edge=None):
    bins = {b: bin_curve(curves[b], anchor) for b in "rg"}
    if edge is None:
        avail = []
        if bins["r"]:
            for e in range(min(bins["r"]), max(bins["r"]) + 30, 30):
                q = change(bins["r"], e)["after_minus_before_mag"]
                if q is not None:
                    avail.append((abs(q), -e))
        if not avail:
            return {"status": "insufficient bins"}
        edge = -max(avail)[1]
    solo = {d: v for d, v in curves["g"].items() if d not in curves["r"]}
    return {
        "boundary": str(EPOCH + timedelta(days=edge)),
        "boundary_mjd": edge,
        "r": change(bins["r"], edge),
        "g": change(bins["g"], edge),
        "green_without_red": change(bin_curve(solo, anchor), edge),
    }


def curves(name):
    cases = json.loads((BATCH / "results/ledger.json").read_text())["objects"]
    case = next(r for r in cases if r["id"] == name)
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / name)
    return case, previous.nightly(rows, case["jd"] + 1e-7)


def run():
    out = {}
    for name in NAMES:
        case, pre = curves(name)
        anchor = date.fromisoformat(case["date"])
        cutoff, scanned, scores, waves = get_models(pre, name)
        raw = state(pre, anchor)
        item = {
            "anchor": str(anchor),
            "training_last": cutoff,
            "frequencies_scanned": scanned,
            "bands": {},
            "clocks": scores,
            "state": raw,
        }
        for b in "rg":
            dd = list(pre[b])
            v = [m for _, m in pre[b].values()]
            item["bands"][b] = {
                "first": str(dd[0]),
                "last": str(dd[-1]),
                "nights": len(dd),
                "p05_p50_p95": list(np.round(np.quantile(v, [0.05, 0.5, 0.95]), 5)),
            }
        if "boundary_mjd" in raw:
            residuals = {}
            for b in "rg":
                orig, f, c = waves["early_red"][b]
                dd = list(pre[b])
                tt, yy = arrays(pre[b], dd)
                expected = previous.design(tt - orig, f) @ c
                residuals[b] = dict(zip(dd, zip(tt, yy - expected)))
            item["blind_clock_residuals_at_same_boundary"] = state(
                residuals, anchor, raw["boundary_mjd"]
            )
        out[name] = item
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "results.json")
    args = p.parse_args()
    pilot.json_write(args.out, run())
