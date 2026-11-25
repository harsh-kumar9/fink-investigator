"""Compare fixed candidate frequencies on the same held-out red nights."""

import json
import numpy as np
from scipy.stats import spearmanr

from epochs import Ephemeris, HERE


def run():
    eph = Ephemeris()
    t = eph.data[eph.n :, 1] - eph.epoch
    y = eph.data[eph.n :, 2]
    const = np.mean(eph.data[: eph.n, 2])
    catalog = json.loads((HERE / "catalog/results.json").read_text())["atlas"]
    frequencies = {
        "red_fit": eph.frequency,
        "atlas_lomb_scargle": 1 / float(catalog["fp-LSper"]),
        "atlas_lombscar": 1 / float(catalog["ls-Pday"]),
        "red_minus_sidereal_alias": 1.0027379 - eph.frequency,
        "red_plus_sidereal_alias": 1.0027379 + eph.frequency,
        "double_red_period": eph.frequency / 2,
        "half_red_period": eph.frequency * 2,
    }
    scores = {}
    for name, f in frequencies.items():
        pred = eph.model.model(t, f)
        mse = np.mean((y - pred) ** 2)
        scores[name] = {
            "period_days": 1 / f,
            "rmse_mag": float(np.sqrt(mse)),
            "reduction_vs_constant": float(1 - mse / np.mean((y - const) ** 2)),
            "rho": float(spearmanr(y, pred).statistic),
        }
    phase = np.arange(10001) / 10000
    model = eph.model.model(phase / eph.frequency, eph.frequency)
    return {
        "test_nights": len(t),
        "training_origin_hjd": eph.epoch,
        "training_mean_mag": const,
        "training_model_peak_to_peak_mag": float(np.ptp(model)),
        "model_brightest_phase": float(phase[np.argmin(model)]),
        "model_faintest_phase": float(phase[np.argmax(model)]),
        "constant_rmse_mag": float(np.sqrt(np.mean((y - const) ** 2))),
        "fits": scores,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
