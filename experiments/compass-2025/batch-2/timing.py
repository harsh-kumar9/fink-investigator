"""Predict the later pulsations of LINEAR 15938573 with single- and two-clock fits."""

from pathlib import Path
import json
import sys
from astropy.timeseries import LombScargle
import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from compass_direct import direct_data

NAME = "ZTF18aagstwh"
SPLIT = 59945  # MJD of 2023-01-01
FAMILIES = ("radial", "nearby")


def series():
    objects = json.loads((HERE / "results/ledger.json").read_text())["objects"]
    (case,) = [r for r in objects if r["id"] == NAME and r["selected_for"] == "fading"]
    _, rows, _ = direct_data(case, HERE / "direct" / NAME)
    data = {}
    for band in ("g", "r"):
        nights = {}
        for r in rows:
            if r["filtercode"] != "z" + band or int(r["catflags"]) != 0:
                continue
            mjd, t, y, e = (float(r[k]) for k in ("mjd", "hjd", "mag", "magerr"))
            if not (0 < e <= 0.2 and np.isfinite([t, mjd, y]).all()):
                continue
            day = int(np.floor(mjd))
            candidate = (e, int(r["expid"]), t - 2458000.0, y)
            if day not in nights or candidate < nights[day]:
                nights[day] = candidate
        data[band] = np.array(
            [(d, t, y, e) for d, (e, _, t, y) in sorted(nights.items())]
        )
    return data


def design(t, f, second=None):
    phase = 2 * np.pi * np.outer(t, f * np.arange(1, 6))
    cols = [np.ones(len(t)), np.sin(phase), np.cos(phase)]
    if second is not None:
        theta = 2 * np.pi * t * second
        cols.extend([np.sin(theta), np.cos(theta)])
    return np.column_stack(cols)


def fit(t, y, f, second=None):
    X = design(t, f, second)
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    return coef, y - X @ coef


def peak(t, y, intervals):
    ls = LombScargle(t, y)
    peaks = []
    for lo, hi in intervals:
        n = int(np.ceil((hi - lo) * np.ptp(t) * 12)) + 1
        fs = np.linspace(lo, hi, n)
        power = ls.power(fs, method="fast", assume_regular_frequency=True)
        i = np.argmax(power)
        peaks.append((float(power[i]), fs[i], lo, hi))
    return max(peaks)


def optimize(t, y, intervals, primary=None):
    residual = y if primary is None else fit(t, y, primary)[1]
    _, mid, lo, hi = peak(t, residual, intervals)

    def objective(freq):
        resid = fit(t, y, freq)[1] if primary is None else fit(t, y, primary, freq)[1]
        return np.dot(resid, resid)

    step = 1 / np.ptp(t)
    grid = np.linspace(max(lo, mid - step), min(hi, mid + step), 41)
    i = min(range(len(grid)), key=lambda i: objective(grid[i]))
    lower, upper = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    result = minimize_scalar(
        objective, bounds=(lower, upper), method="bounded", options={"xatol": 1e-13}
    )
    return min([lower, upper, mid, result.x], key=objective)


def scores(data):
    masks = {band: d[:, 0] < SPLIT for band, d in data.items()}
    for band, m in masks.items():
        if min(m.sum() / 80, (~m).sum() / 40) < 1:
            raise ValueError("Insufficient chronological coverage for " + band)
    red = data["r"][masks["r"]]
    t, y = red[:, 1], red[:, 2]
    f = optimize(t, y, [(1, 5)])
    fits = {
        "radial": optimize(t, y, [(0.70 * f, 0.80 * f), (f / 0.80, f / 0.70)], f),
        "nearby": optimize(t, y, [(f - 0.10, f - 1 / 400), (f + 1 / 400, f + 0.10)], f),
    }
    res = {"frequency_per_day": f, "period_days": 1 / f, "families": {}, "bands": {}}
    for key, sec in fits.items():
        df = abs(sec - f)
        res["families"][key] = {
            "frequency_per_day": sec,
            "period_days": 1 / sec,
            "frequency_difference": df,
            "beat_days": 1 / df,
            "distance_to_one_sidereal_day_alias_per_day": abs(df - 1.0027379),
        }
    for band, d in data.items():
        mask = masks[band]
        t, y = d[mask, 1], d[mask, 2]
        tt, yy = d[~mask, 1], d[~mask, 2]
        pred = design(tt, f) @ fit(t, y, f)[0]
        err = np.sum((yy - pred) ** 2)
        bandres = {
            "training_dates": int(mask.sum()),
            "test_dates": int((~mask).sum()),
            "primary_error_reduction": 1 - err / np.sum((yy - y.mean()) ** 2),
            "primary_test_rms_mag": np.sqrt(err / len(yy)),
        }
        for key, sec in fits.items():
            coef = fit(t, y, f, sec)[0]
            err2 = np.sum((yy - design(tt, f, sec) @ coef) ** 2)
            bandres[key] = {
                "additional_error_reduction": 1 - err2 / err,
                "training_secondary_semi_amplitude_mag": np.hypot(*coef[-2:]),
            }
        res["bands"][band] = bandres
    for key in fits:
        res["families"][key]["passes_both_bands"] = all(
            b[key]["additional_error_reduction"] >= 0.10 for b in res["bands"].values()
        )
    return res


def run():
    data = series()
    out = {}
    for key, amplitude in [("stable_injection", 0.0), ("extra_mode_injection", 0.12)]:
        fake = {}
        for seed, (band, original) in enumerate(sorted(data.items()), 31):
            copy = original.copy()
            t = copy[:, 1]
            f = 1.7
            phi = 2 * np.pi * f * t
            copy[:, 2] = (
                17
                + 0.35 * np.sin(phi)
                + 0.13 * np.sin(2 * phi - 0.8)
                + amplitude * np.sin(phi / 0.745 + 0.8)
                + np.random.default_rng(seed).normal(0, copy[:, 3])
            )
            fake[band] = copy
        out[key] = scores(fake)
    out["source"] = scores(data)
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
