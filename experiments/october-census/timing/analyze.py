"""Predict held-out ZTF nights and pair historical optical and infrared visits."""

import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

from astropy.timeseries import LombScargle
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
CENSUS = HERE.parent


def nights(case):
    directory = CENSUS / "sources/direct"
    name = case["id"]
    receipt = directory / (name + ".json")
    if receipt.exists():
        meta = json.loads(receipt.read_text())
        path = directory / (name + ".csv.gz")
    else:
        directory = directory / "recovery"
        meta = json.loads((directory / (name + ".json")).read_text())["photometry"]
        path = directory / (name + ".photometry.csv.gz")
    raw = gzip.decompress(path.read_bytes())
    if (
        meta["url"] != case["query"]
        or hashlib.sha256(raw).hexdigest() != meta["sha256_uncompressed"]
    ):
        raise ValueError(f"Unexpected photometry for {name}")
    best = {}
    for row in csv.DictReader(io.StringIO(raw.decode())):
        if (
            row["oid"] != case["direct_source_id"]
            or row["filtercode"] != "zr"
            or int(row["catflags"])
        ):
            continue
        error, mjd, hjd, mag = (float(row[k]) for k in ("magerr", "mjd", "hjd", "mag"))
        if not 0 < error <= 0.2 or not all(map(math.isfinite, (mjd, hjd, mag))):
            continue
        date = math.floor(mjd)
        value = (error, int(row["expid"]), mjd, hjd, mag)
        if date not in best or value < best[date]:
            best[date] = value
    return np.array([best[k][2:] for k in sorted(best)]).reshape(-1, 3)


def search(data):
    n = len(data) * 2 // 3
    out = {"training_nights": n, "test_nights": len(data) - n}
    if n < 60 or len(data) - n < 30:
        return {**out, "status": "insufficient dates"}
    t, y = data[:, 1] - data[0, 1], data[:, 2]
    model = LombScargle(t[:n], y[:n])
    freq, power = model.autopower(
        minimum_frequency=1 / 300,
        maximum_frequency=10,
        samples_per_peak=10,
        method="fast",
    )
    upper = np.searchsorted(freq, 10, side="right")
    index = np.argmax(power[:upper])
    lo = freq[max(0, index - 1)]
    hi = min(10.0, freq[index] + freq[1] - freq[0])
    fit = minimize_scalar(
        lambda f: -model.power(f, method="cython"),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": 1e-13},
    )
    f = max([lo, fit.x, hi], key=lambda v: model.power(v, method="cython"))
    waveform = LombScargle(t[:n], y[:n], nterms=2)
    pred = waveform.model(t[n:], f)
    gain = 1 - np.sum((y[n:] - pred) ** 2) / np.sum((y[n:] - np.mean(y[:n])) ** 2)
    rho = float(spearmanr(y[n:], pred).statistic)
    cycles = float(t[n - 1] * f)
    return {
        **out,
        "training_last_mjd": float(data[n - 1, 0]),
        "test_first_mjd": float(data[n, 0]),
        "test_last_mjd": float(data[-1, 0]),
        "period_days": 1 / f,
        "training_power": float(model.power(f, method="cython")),
        "training_cycles": cycles,
        "test_error_reduction": float(gain),
        "test_spearman": rho,
        "status": (
            "passes"
            if cycles >= 3 and gain >= 0.30 and rho >= 0.60
            else "below threshold"
        ),
    }


def pair_visits(data, source):
    pairs = []
    for visit in source["visits"]:
        t = visit["median_mjd"]
        matched = data[np.abs(data[:, 0] - t) <= 10]
        if len(matched) >= 2:
            pairs.append(
                {
                    "mjd": t,
                    "w1": visit["median_w1"],
                    "r": float(np.median(matched[:, 2])),
                    "optical_nights": len(matched),
                }
            )
    result = {
        "designation": source["designation"],
        "role": source["role"],
        "pairs": pairs,
    }
    if len(pairs) >= 5:

        def rho(points):
            return float(
                spearmanr([p["w1"] for p in points], [p["r"] for p in points]).statistic
            )

        leave_one_out = [rho(pairs[:i] + pairs[i + 1 :]) for i in range(len(pairs))]
        result.update(
            spearman=rho(pairs),
            leave_one_out_range=[min(leave_one_out), max(leave_one_out)],
        )
    return result


def run():
    cohort = json.loads((CENSUS / "results/direct.json").read_text())
    cases = [cohort["control"], *cohort["exploration"]]
    rows = [
        {"id": c["id"], "oid": c["direct_source_id"], **search(nights(c))}
        for c in cases
    ]
    ir_target = next(c for c in cases if c["id"] == "ZTF17aacydyv")
    data = nights(ir_target)
    infrared = json.loads((CENSUS / "infrared/results.json").read_text())
    return {
        "calibration_passed": rows[0]["status"] == "passes",
        "period_search": rows,
        "optical_infrared": [pair_visits(data, s) for s in infrared],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
