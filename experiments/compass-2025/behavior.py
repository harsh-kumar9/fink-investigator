"""Count resolved two-band eruptions and test a published quiescent timescale."""

from collections import defaultdict
from datetime import date, timedelta
import json
import math
from pathlib import Path
from statistics import mean, median, quantiles
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import compass_direct
from census_photometry import correlation

EPOCH = date(1858, 11, 17)
FT = "ZTF18abodmfy"
CYG = "ZTF18aazenpc"


def curves(name):
    ledger = json.loads((HERE / "results/ledger.json").read_text())
    case = next(r for r in ledger["objects"] if r["id"] == name)
    _, rows, _ = compass_direct.direct_data(case, HERE / "direct" / name)
    anchor = date.fromisoformat(case["date"])
    return anchor, {
        band: compass_direct.nightly(rows, band, case["jd"])[0] for band in "gr"
    }


def eruptions(data):
    red, green = data["r"], data["g"]
    r0, g0 = median(red.values()), median(green.values())
    bright_g = [d for d, mag in green.items() if mag <= g0 - 2]
    groups = []
    points = []
    for d, mag in sorted(red.items()):
        if mag <= r0 - 2:
            points.append(d)
        elif mag >= r0 - 0.5 and points:
            groups.append(points)
            points = []
    if points:
        groups.append(points)
    episodes = []
    for dates in groups:
        supported = any((b - a).days <= 10 for a, b in zip(dates, dates[1:]))
        matches = [d for d in bright_g if any(abs((d - x).days) <= 2 for x in dates)]
        episodes.append(
            {
                "first_bright": str(dates[0]),
                "last_bright": str(dates[-1]),
                "bright_red_dates": len(dates),
                "brightest_r_mag": min(red[d] for d in dates),
                "resolved": supported,
                "green_matches": [str(d) for d in sorted(matches)],
            }
        )
    return {
        "median_r": r0,
        "median_g": g0,
        "red_nights": len(red),
        "green_nights": len(green),
        "groups": episodes,
        "resolved": sum(r["resolved"] for r in episodes),
        "also_green": sum(r["resolved"] and bool(r["green_matches"]) for r in episodes),
    }


def bins(curve, quantile, anchor):
    groups = defaultdict(list)
    for day, mag in curve.items():
        mjd = (day - EPOCH).days
        start = 58000 + (mjd - 58000) // 30 * 30
        groups[start].append((mjd, mag))
    rows = []
    for start, points in sorted(groups.items()):
        if len(points) < 3 or EPOCH + timedelta(days=start + 30) > anchor:
            continue
        t, y = zip(*points)
        first, stop = (EPOCH + timedelta(days=x) for x in (start, start + 30))
        subset = (
            "training"
            if stop <= date(2023, 1, 1)
            else "test" if first >= date(2024, 1, 1) else "unused"
        )
        rows.append(
            {
                "mjd": median(t),
                "mag": quantiles(y, n=100, method="inclusive")[
                    round(100 * quantile) - 1
                ],
                "n": len(t),
                "start_mjd": start,
                "subset": subset,
            }
        )
    return rows


def predict(binned, period):
    def design(row):
        angle = (row["mjd"] - 58000) * 2 * math.pi / period
        return math.sin(angle), math.cos(angle), row["mag"]

    train = [design(r) for r in binned if r["subset"] == "training"]
    test = [design(r) for r in binned if r["subset"] == "test"]
    if len(train) < 3 or len(test) < 2:
        return {"status": "insufficient data"}
    # Remove the intercept, then solve the 2 × 2 sine/cosine normal equations.
    s0, c0, y0 = map(mean, zip(*train))
    centered = [(s - s0, c - c0, y - y0) for s, c, y in train]
    ss = sum(s * s for s, c, y in centered)
    cc = sum(c * c for s, c, y in centered)
    sc = sum(s * c for s, c, y in centered)
    sy = sum(s * y for s, c, y in centered)
    cy = sum(c * y for s, c, y in centered)
    det = ss * cc - sc * sc
    if det <= 1e-12:
        return {"status": "insufficient phase coverage"}
    a, b = (sy * cc - cy * sc) / det, (cy * ss - sy * sc) / det
    offset = y0 - a * s0 - b * c0
    actual = [y for s, c, y in test]
    pred = [offset + a * s + b * c for s, c, y in test]
    baseline = sum((y - y0) ** 2 for y in actual)
    gain = (
        1 - sum((y - p) ** 2 for y, p in zip(actual, pred)) / baseline
        if baseline
        else None
    )
    rho = correlation(actual, pred)
    return {
        "train_bins": len(train),
        "test_bins": len(test),
        "half_amplitude_mag": math.hypot(a, b),
        "coefficients": [offset, a, b],
        "test_error_reduction": gain,
        "test_spearman": rho,
        "passes": len(test) >= 8
        and gain is not None
        and gain >= 0.2
        and rho is not None
        and rho >= 0.5,
    }


def run():
    _, ft = curves(FT)
    anchor, cyg = curves(CYG)
    checks = {}
    for band, curve in cyg.items():
        main = bins(curve, 0.75, anchor)
        checks[band] = {
            "bins": main,
            "main_572_days": predict(main, 572),
            "sensitivity_fixed_period": {
                str(p): predict(main, p) for p in [547, 597, 365.25]
            },
            "sensitivity_quantile": {
                str(q): predict(bins(curve, q, anchor), 572) for q in [0.5, 0.9]
            },
        }
    return {"ft_cam": eruptions(ft), "v1060_cyg": checks}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
