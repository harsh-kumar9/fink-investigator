"""Predict April 6 from the first two April nights of the AAVSO outburst light curve."""

import json
import math
import operator
import random
from collections import Counter
from statistics import median, quantiles

from sources import aavso, catalog

FREQUENCIES = [4 + i / 100 for i in range(3201)]
dot = getattr(math, "sumprod", lambda a, b: sum(map(operator.mul, a, b)))


def selected(rows):
    return sorted(
        (r for r in rows if r["by"] == "HMB" and r["band"] == "CV" and r["val"] in {"V", "Z"}
         and r["fainterThan"] == "0" and r["obsType"] == "CCD"
         and 0 < float(r["uncert"]) <= 0.2 and math.isfinite(float(r["mag"]))),
        key=lambda r: float(r["JD"]),
    )


def residualizer(times, nights, slopes):
    pieces = []
    for night in sorted(set(nights)):
        indices = [i for i, value in enumerate(nights) if value == night]
        mean = sum(times[i] for i in indices) / len(indices)
        centered = [times[i] - mean for i in indices]
        pieces.append((indices, centered, dot(centered, centered)))

    def apply(values):
        result = list(values)
        for indices, centered, variance in pieces:
            vals = [values[i] for i in indices]
            level = sum(vals) / len(vals)
            slope = dot(vals, centered) / variance if slopes else 0
            for i, x in zip(indices, centered):
                result[i] -= level + slope * x
        return result
    return apply


def components(times, frequency):
    angle = 2 * math.pi * frequency
    return [math.sin(angle * t) for t in times], [math.cos(angle * t) for t in times]


def design(times, remove, frequency):
    sine, cosine = [remove(v) for v in components(times, frequency)]
    a = math.sqrt(dot(sine, sine))
    u = [x / a for x in sine]
    b = dot(cosine, u)
    centered = [c - b*x for c, x in zip(cosine, u)]
    d = math.sqrt(dot(centered, centered))
    v = [x / d for x in centered]
    return u, v, a, b, d


class Test:
    def __init__(self, times, nights, slopes):
        self.split = nights.index(max(nights))
        self.span = times[-1] - times[self.split]
        self.test_times = times[self.split:]
        self.train_remove = residualizer(times[:self.split], nights[:self.split], slopes)
        self.test_remove = residualizer(self.test_times, nights[self.split:], slopes)
        self.grid = [design(times[:self.split], self.train_remove, f) for f in FREQUENCIES]

    def scan(self, values):
        y = self.train_remove(values[:self.split])
        total = dot(y, y)
        error = [total - dot(u, y)**2 - dot(v, y)**2 for u, v, *_ in self.grid]
        return error

    def fit(self, values, index):
        y = self.train_remove(values[:self.split])
        u, v, a, b, d = self.grid[index]
        cos_coef = dot(v, y) / d
        sin_coef = (dot(u, y) - cos_coef*b) / a
        frequency = FREQUENCIES[index]
        sin1, cos1 = components(self.test_times, frequency)
        test_y = self.test_remove(values[self.split:])
        signal = self.test_remove([sin_coef*s + cos_coef*c for s, c in zip(sin1, cos1)])
        residual = [m - s for m, s in zip(test_y, signal)]
        train_gain = (dot(u, y)**2 + dot(v, y)**2) / dot(y, y)
        heldout_gain = 1 - dot(residual, residual) / dot(test_y, test_y)
        amplitude = math.hypot(sin_coef, cos_coef)
        cycles = self.span * frequency
        return {
            "frequency_per_day": round(frequency, 5), "period_minutes": round(1440/frequency, 4),
            "training_semiamplitude_mag": round(amplitude, 5),
            "training_error_reduction": round(train_gain, 5), "heldout_error_reduction": round(heldout_gain, 5),
            "heldout_cycles_covered": round(cycles, 4),
            "transfer_pass": train_gain >= 0.2 and heldout_gain >= 0.2 and amplitude > 0.02 and cycles >= 1,
        }

    def alternatives(self, values, errors, best):
        frequency = FREQUENCIES[best]
        results = {}
        for label, center, tolerance in [
            ("minus_two_daily", frequency-2, 0.4), ("minus_one_daily", frequency-1, 0.4),
            ("plus_one_daily", frequency+1, 0.4), ("plus_two_daily", frequency+2, 0.4),
            ("half_frequency", frequency/2, 0.00501), ("double_frequency", frequency*2, 0.00501),
        ]:
            if not FREQUENCIES[0] <= center <= FREQUENCIES[-1]:
                continue
            nearby = [i for i, f in enumerate(FREQUENCIES) if abs(f-center) <= tolerance+1e-6]
            ix = min(nearby, key=errors.__getitem__)
            results[label] = {**self.fit(values, ix), "training_error_ratio_to_best": round(errors[ix]/errors[best], 5)}
        return results


def injections(times, errors, test):
    rng = random.Random(1729)
    output = {}
    for amplitude in (0, 0.03, 0.1):
        counts = Counter()
        for _ in range(60):
            phase = rng.uniform(0, math.tau)
            y = [amplitude*math.sin(math.tau*16*t + phase) + rng.gauss(0, err) for t, err in zip(times, errors)]
            scores = test.scan(y)
            best = min(range(len(scores)), key=scores.__getitem__)
            passed = test.fit(y, best)["transfer_pass"]
            close = abs(FREQUENCIES[best]-16) <= 0.3 + 1e-8
            counts["transfer_pass"] += passed
            counts["recovered_frequency"] += close
            counts["transfer_and_recovered"] += passed and close
        output[str(amplitude)] = {"draws": 60, **counts}
    return output


def run():
    all_rows = aavso()
    rows = selected(all_rows)
    jds = [float(r["JD"]) for r in rows]
    nights = [str(catalog.pilot.day_of(t)) for t in jds]
    coverage_cut = {}
    for d in sorted({str(catalog.pilot.day_of(float(r["JD"]))) for r in all_rows}):
        e = [float(r["uncert"]) for r in all_rows if str(catalog.pilot.day_of(float(r["JD"]))) == d]
        coverage_cut[d] = {"received": len(e), "passed": nights.count(d),
                           "min_error": min(e), "median_error": round(median(e), 5), "max_error": max(e)}
    if sorted(set(nights)) != ["2025-04-04", "2025-04-05"] or coverage_cut["2025-04-06"]["passed"] != 0:
        raise ValueError("Coverage differs from the registered two-night amendment")
    times = [t-jds[0] for t in jds]
    mags = [float(r["mag"]) for r in rows]
    errs = [float(r["uncert"]) for r in rows]
    test = Test(times, nights, slopes=True)
    errors = test.scan(mags)
    best = min(range(len(errors)), key=errors.__getitem__)
    result = test.fit(mags, best)
    coverage = {}
    for night in sorted(set(nights)):
        which = [i for i, n in enumerate(nights) if n == night]
        t = [times[i] for i in which]
        gaps = [(b-a)*1440 for a, b in zip(t, t[1:])]
        span = t[-1] - t[0]
        accepted_mags = [mags[i] for i in which]
        quartiles = quantiles(accepted_mags, n=4, method="inclusive")
        coverage[night] = {
            "exposures": len(which), "first_jd": jds[which[0]], "last_jd": jds[which[-1]],
            "span_hours": round(24*span, 3), "cycles_at_peak": round(span*FREQUENCIES[best], 4),
            "cadence_median_minutes": round(median(gaps), 3), "largest_gap_minutes": round(max(gaps), 3),
            "median_mag": round(median(accepted_mags), 4),
            "mag_q25_q75": [round(quartiles[i], 4) for i in (0, 2)],
            "median_check_minus_comp": round(median(float(rows[i]["kmag"]) - float(rows[i]["cmag"]) for i in which), 4),
            "median_err_mag": round(median(errs[i] for i in which), 4),
        }
    other = Test(times, nights, slopes=False)
    scores = other.scan(mags)
    alt_best = min(range(len(scores)), key=scores.__getitem__)
    star_pairs = {(r["compStar1"], r["compStar2"]) for r in rows}
    if len(star_pairs) != 1:
        raise ValueError("Field stars changed among observations")
    check = [float(r["kmag"]) - float(r["cmag"]) for r in rows]
    result.update(
        registered_three_night_test="abstained: no qualifying April 6 exposure",
        two_night_amendment="April 4 selects the frequency; April 5 tests the locked oscillation",
        coverage_before_and_after_quality=coverage_cut,
        nights=coverage, alternatives=test.alternatives(mags, errors, best),
        offset_only_best=other.fit(mags, alt_best), offset_only_at_registered_peak=other.fit(mags, best),
        comp_check_ids=list(next(iter(star_pairs))), field_control_at_registered_peak=test.fit(check, best),
        injections=injections(times, errs, test),
    )
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
