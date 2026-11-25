"""Characterize narrow fades in one field and predict their timing in another."""

import csv
import gzip
import json
from astropy.timeseries import BoxLeastSquares
import numpy as np

from analyze import CENSUS, nights
from epochs import Ephemeris, NAME


def predictions():
    eph = Ephemeris()
    train = eph.data[: eph.n]
    model = BoxLeastSquares(eph.train_t, -train[:, 2])
    families = {
        "red_clock": eph.frequency,
        "lower_daily_alias": 1.0027379 - eph.frequency,
        "upper_daily_alias": 1.0027379 + eph.frequency,
        "half_frequency": eph.frequency / 2,
        "double_frequency": eph.frequency * 2,
    }
    fits = {}
    for name, f in families.items():
        count = int(np.ceil(0.02 * f * np.ptp(eph.train_t) / 0.02)) + 1
        periods = 1 / np.linspace(1.01 * f, 0.99 * f, count)
        durations = np.arange(1, 6) * 0.02 / f
        result = model.power(periods, durations, objective="likelihood", oversample=20)
        i = np.argmax(result.power)
        fits[name] = {
            "period_days": result.period[i],
            "duration_days": result.duration[i],
            "center_hjd": eph.epoch + result.transit_time[i],
            "training_depth_mag": result.depth[i],
            "training_power": result.power[i],
        }
    winner = max(fits, key=lambda k: fits[k]["training_power"])
    return eph, fits, winner


def evaluate(data, fit):
    dt = (data[:, 1] - fit["center_hjd"]) / fit["period_days"]
    inside = np.abs(dt - np.round(dt)) < fit["duration_days"] / (2 * fit["period_days"])
    dim = data[inside]
    bright = data[~inside, 2]
    return {
        "inside": len(dim),
        "outside": len(bright),
        "in_window_mjd_mag": [[t, mag] for t, _, mag in dim],
        "outside_median_mag": float(np.median(bright)) if len(bright) else None,
        "inside_minus_outside_median_mag": (
            float(np.median(dim[:, 2]) - np.median(bright))
            if len(dim) and len(bright)
            else None
        ),
        "outside_p90_minus_p10_mag": (
            float(np.ptp(np.percentile(bright, [10, 90]))) if len(bright) else None
        ),
    }


def run():
    eph, fits, winner = predictions()
    cohort = json.loads((CENSUS / "results/direct.json").read_text())["exploration"]
    case = next(c for c in cohort if c["id"] == NAME)
    with gzip.open(CENSUS / "sources/direct" / (NAME + ".csv.gz"), "rt") as file:
        ids = sorted(
            {r["oid"] for r in csv.DictReader(file) if r["filtercode"] == "zr"}
        )
    used = np.floor(eph.data[:, 0])
    other = {}
    for oid in ids:
        if oid == case["direct_source_id"]:
            continue
        data = nights({**case, "direct_source_id": oid})
        data = data[~np.isin(np.floor(data[:, 0]), used)]
        other[oid] = evaluate(data, fits[winner])
    return {
        "fits": fits,
        "training_winner": winner,
        "original_later": evaluate(eph.data[eph.n :], fits[winner]),
        "other_series_disjoint_dates": other,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
