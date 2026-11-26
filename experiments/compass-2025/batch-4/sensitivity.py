"""Propagate per-night flux errors and reference choice into the temporal measurements."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from random import Random
from statistics import median

from comparison import HERE, dist, features, group_forced, inputs, pilot, windows


def quantile(values, q):
    s = sorted(values)
    x = (len(s) - 1) * q
    i = int(x)
    return s[i] + (x - i) * (s[min(i + 1, len(s) - 1)] - s[i])


def scales(red_mag, observed_days):
    # Flat Lambda-CDM. Simpson's rule integrates c/H0/E(z) in comoving Mpc.
    out = {}
    for z in (0.1, 0.25, 0.440, 0.534, 0.623):
        n = 1000
        dz = z / n
        integral = (
            sum(
                (1 if i in (0, n) else 4 if i % 2 else 2)
                / math.sqrt(0.3 * (1 + i * dz) ** 3 + 0.7)
                for i in range(n + 1)
            )
            * dz
            / 3
        )
        dl = (1 + z) * 299792.458 / 70 * integral
        out[str(z)] = {
            "luminosity_distance_mpc": dl,
            "apparent_r_minus_distance_modulus": red_mag - 25 - 5 * math.log10(dl),
            "span_days_divided_by_1_plus_z": observed_days / (1 + z),
        }
    return out


def run():
    raw, _, _ = inputs(False)
    ns = group_forced(raw)
    keep = [r for r in ns if r["mag"] is not None and r["magerr"] <= 0.2]
    controls = json.loads((HERE / "results/comparison.json").read_text())["panels"][
        "2024_bright"
    ]["controls"]
    controls = [r for r in controls if r["all"]]
    rng = Random(19283)
    measurements = []
    near = Counter()
    peaks = Counter()
    for _ in range(2000):
        flux = {(r["day"], r["fid"]): rng.gauss(r["flux"], r["error"]) for r in keep}
        curve = {k: 23.9 - 2.5 * math.log10(v) for k, v in flux.items() if v > 0}
        peak, win = windows(curve, 60790)
        m = features(win)
        if m:
            measurements.append(m)
            near[
                min(controls, key=lambda r: (dist(m, r["all"]), r["id"]))["family"]
            ] += 1
            peaks[peak] += 1
    individual = {}
    for field in sorted({r["field"] for r in raw}):
        v = group_forced([r for r in raw if r["field"] == field])
        c = {
            (r["day"], r["fid"]): r["mag"]
            for r in v
            if r["mag"] is not None and r["magerr"] <= 0.2
        }
        # Compare the same temporal bins across fields.
        day = int(
            json.loads((HERE / "results/comparison.json").read_text())[
                "target_red_peak_mjd_day"
            ]
        )
        color = {
            d - day: g - c[d, 2] for (d, f), g in c.items() if f == 1 and (d, 2) in c
        }
        a = [x for d, x in color.items() if -10 <= d <= 2]
        b = [x for d, x in color.items() if 20 <= d <= 35]
        individual[field] = {
            "early_nights": len(a),
            "late_nights": len(b),
            "early_gr": median(a) if a else None,
            "late_gr": median(b) if b else None,
            "change": median(b) - median(a) if a and b else None,
        }
    measured = json.loads((HERE / "results/onset.json").read_text())["half_peak"]["2"]
    return {
        "draws": 2000,
        "valid": len(measurements),
        "conditional_p16_p50_p84": {
            k: [quantile([m[k] for m in measurements], q) for q in (0.16, 0.5, 0.84)]
            for k in ("decline", "reddening")
        },
        "red_peak_mjd_days": dict(sorted(peaks.items())),
        "closest_class_counts": dict(near),
        "fields": individual,
        "conditional_distance_scales": scales(
            measured["observed_peak"]["mag"],
            measured["crossings_mjd"]["0"]["span_above_days"],
        ),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "results/sensitivity.json")
    args = p.parse_args()
    pilot.json_write(args.out, run())
