"""Test the published MGAB-V488 period against later ZTF observations and injections."""

import argparse
from collections import defaultdict
from datetime import date
from itertools import combinations
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from measure import HERE, BATCH, load, previous, pilot, compass_direct

NAME = "ZTF18aaadshm"
PERIOD = 0.4904450424
AMPLITUDE = {"r": 1.20, "g": 1.39}
PHASES = np.arange(32) * (2 * np.pi / 32)


def amplitude(coef):
    phase = np.arange(2048) / 2048
    return round(float(np.ptp(previous.design(phase, 1) @ coef)), 6)


def metrics(y, predicted, level, minimum):
    if len(y) < 2:
        return {"n": len(y), "passes": False}
    gain = 1 - np.sum((y - predicted) ** 2) / np.sum((y - level) ** 2)
    rho = spearmanr(y, predicted).statistic
    return {
        "n": len(y),
        "mse_gain": round(float(gain), 6),
        "rho": round(float(rho), 6),
        "passes": bool(len(y) >= minimum and gain >= 0.2 and rho >= 0.5),
    }


def blocks(curve, cutoff, anchor, origin):
    dates = list(curve)
    return {
        label: np.array(
            [(t - origin, m) for d in dates if test(d) for t, m in [curve[d]]]
        )
        .reshape(-1, 2)
        .T
        for label, test in [
            ("early", lambda d: d <= cutoff),
            ("later", lambda d: cutoff < d <= anchor),
            ("post_anchor", lambda d: d > anchor),
        ]
    }


def fits(groups, frequency):
    red_t, red_y = groups["r"]["early"]
    _, red = previous.fit(red_t, red_y, frequency)
    out = {}
    for b, sets in groups.items():
        t, y = sets["early"]
        red_shape = previous.design(t, frequency) @ red
        if b == "g":
            scale, off = np.linalg.lstsq(
                np.column_stack([red_shape, np.ones(len(y))]), y, rcond=None
            )[0]
            tied = red * scale
            tied[0] += off
        else:
            tied = red
        _, own = previous.fit(t, y, frequency)
        item = {}
        for window, (tx, yx) in sets.items():
            matrix = previous.design(tx, frequency)
            item[window] = {
                "red_waveform": metrics(yx, matrix @ tied, np.mean(y), 30),
                "band_waveform": metrics(yx, matrix @ own, np.mean(y), 30),
                "band_refit_peak_to_peak_mag": (
                    amplitude(previous.fit(tx, yx, frequency)[1])
                    if len(tx) >= 5
                    else None
                ),
            }
        out[b] = item
    return out


def within_night(rows, cutoff, anchor, origin):
    by_date = defaultdict(dict)
    for row in rows:
        t = float(row["hjd"]) - origin
        mag = float(row["mag"])
        if (
            int(row["catflags"])
            or not 0 < float(row["magerr"]) <= 0.2
            or not all(map(math.isfinite, (mag, t)))
        ):
            continue
        day = pilot.day_of(float(row["mjd"]) + 2400000.5)
        if day <= anchor:
            by_date[row["filtercode"][-1], day][row["expid"]] = (t, mag)
    pairs = defaultdict(list)
    for (band, day), exps in sorted(by_date.items()):
        choices = [
            (t2 - t1, t1, t2, y2 - y1)
            for (t1, y1), (t2, y2) in combinations(sorted(exps.values()), 2)
            if 0.05 <= t2 - t1 <= 0.35
        ]
        if choices:
            _, t1, t2, dy = max(choices)
            pairs[band].append((day, t1, t2, dy))
    out = {}
    for b in "rg":
        records = pairs[b]
        designs = {}
        for label, pred in [
            ("early", lambda d: d <= cutoff),
            ("later", lambda d: d > cutoff),
        ]:
            a = np.array(
                [[t1, t2, dy] for d, t1, t2, dy in records if pred(d)]
            ).reshape(-1, 3)
            w = 2 * np.pi / PERIOD
            dx = np.column_stack(
                [
                    np.sin(w * a[:, 1]) - np.sin(w * a[:, 0]),
                    np.cos(w * a[:, 1]) - np.cos(w * a[:, 0]),
                ]
            )
            designs[label] = (dx, a[:, 2])
        train, y = designs["early"]
        test, ref = designs["later"]
        item = {
            "n_early": len(y),
            "n_later": len(ref),
            "dates": [str(row[0]) for row in records],
        }
        if len(y) >= 2 and len(ref) >= 2:
            coef = np.linalg.lstsq(train, y, rcond=None)[0]
            item.update(
                {
                    "amplitude_mag": round(float(2 * np.linalg.norm(coef)), 6),
                    "gain_on_early": metrics(y, train @ coef, 0, 10),
                    "forecast": metrics(ref, test @ coef, 0, 10),
                    "median_abs_later_pair_difference_mag": round(
                        float(np.median(np.abs(ref))), 6
                    ),
                }
            )
            trials = []
            for phase in PHASES:
                model = 0.5 * AMPLITUDE[b] * np.array([np.cos(phase), np.sin(phase)])
                synthetic = y + train @ model
                fit = np.linalg.lstsq(train, synthetic, rcond=None)[0]
                trial = metrics(ref + test @ model, test @ fit, 0, 10)
                trials.append((trial, 2 * np.linalg.norm(fit)))
            item["injection"] = {
                "phases": len(trials),
                "passes": sum(s["passes"] for s, a in trials),
                "amplitude_range_mag": list(
                    np.round([min(a for s, a in trials), max(a for s, a in trials)], 6)
                ),
            }
        out[b] = item
    return out


def run():
    case, pre, full = load(NAME)
    anchor = date.fromisoformat(case["date"])
    train = list(pre["r"])[: len(pre["r"]) * 2 // 3]
    cutoff = train[-1]
    origin = pre["r"][train[0]][0]
    groups = {b: blocks(full[b], cutoff, anchor, origin) for b in "rg"}
    t, y = groups["r"]["early"]
    freq = 1 / PERIOD
    grid = np.linspace(1 / (PERIOD + 0.002), 1 / (PERIOD - 0.002), 1001)
    best = min(grid, key=lambda f: previous.fit(t, y, f)[0])
    trials = defaultdict(list)
    for phase in PHASES:
        simulated = {
            b: {
                k: (t, y + AMPLITUDE[b] * 0.5 * np.sin(2 * np.pi * freq * t + phase))
                for k, (t, y) in windows.items()
            }
            for b, windows in groups.items()
        }
        result = fits(simulated, freq)
        for b in "rg":
            trials[b].append(result[b])
    injection = {}
    for b, scores in trials.items():
        amps = [s["early"]["band_refit_peak_to_peak_mag"] for s in scores]
        injection[b] = {
            "phases": len(scores),
            "early_amplitude_range_mag": [min(amps), max(amps)],
            "passes": {
                window: {
                    method: sum(s[window][method]["passes"] for s in scores)
                    for method in ("red_waveform", "band_waveform")
                }
                for window in ("later", "post_anchor")
            },
        }
    _, raw, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    return {
        "object": NAME,
        "published_period_days": PERIOD,
        "published_peak_to_peak_mag": AMPLITUDE,
        "t_origin_hjd": origin,
        "training_through": str(cutoff),
        "anchor": str(anchor),
        "exact_period": fits(groups, freq),
        "nearby_red_training_best": {
            "period_days": float(1 / best),
            "scores": fits(groups, best),
        },
        "injection": injection,
        "within_night": within_night(raw, cutoff, anchor, origin),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "pulsation.json")
    pilot.json_write(parser.parse_args().out, run())
