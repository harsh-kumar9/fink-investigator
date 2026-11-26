"""Test MGAB-V960's 2019 sequence on external photometry and other seasons."""

import argparse
from collections import defaultdict
from datetime import date
from itertools import combinations
import json, math
from pathlib import Path
from statistics import median

import numpy as np
from scipy.stats import spearmanr

from analyze_vars import HERE, load, pilot
from duration import SRC
from asas import obtain

NAME = "ZTF18abcckqd"
START, STOP = date(2019, 8, 1), date(2019, 10, 31)
LAGS = [(8, 12), (17, 23), (28, 32)]


def clean(row, anchor):
    return (
        date(2019, 1, 1) <= pilot.day_of(row["jd"]) <= anchor
        and row["quality"] == "G"
        and all(math.isfinite(row[k]) for k in ("flux", "flux_err"))
        and row["flux_err"] > 0
    )


def daily(rows):
    by = defaultdict(list)
    for v in rows:
        by[pilot.day_of(v["jd"])].append(v["flux"])
    return {d: median(v) for d, v in sorted(by.items())}


def lag(values):
    out = {}
    for yr in sorted({d.year for d in values}):
        dd = sorted(d for d in values if d.year == yr)
        bundles = defaultdict(lambda: [[], [], []])
        for a, b in combinations(dd, 2):
            days = (b - a).days
            for i, (lo, hi) in enumerate(LAGS):
                if lo <= days <= hi:
                    bundles[a][i].append(abs(values[a] - values[b]))
        scores = [list(map(median, v)) for v in bundles.values() if all(v)]
        m = np.median(scores, axis=0) if scores else None
        out[str(yr)] = {
            "nights": len(dd),
            "anchors": len(scores),
            "delta_mjy": list(np.round(m, 6)) if scores else None,
            "middle_lower": bool(m is not None and m[1] < min(m[0], m[2])),
        }
    return out


def ephem():
    d = json.loads((HERE / "duration.json").read_text())[NAME]["tests"]["primary"][
        "groups"
    ]
    starts = [
        date.fromisoformat(v["first"])
        for v in d
        if v["kind"] == "short" and START <= date.fromisoformat(v["first"]) <= STOP
    ]
    if len(starts) != 5:
        raise ValueError("Five frozen 2019 brackets expected")
    yy = np.array([(s - START).days for s in starts])
    period, epoch = np.polyfit(np.arange(3), yy[:3], 1)
    return starts, round(float(period), 8), round(float(epoch), 8)


def phases(data, period, epoch, shift):
    aa = defaultdict(lambda: [[], []])
    for d, flux in sorted(data.items()):
        if not START <= d <= STOP:
            continue
        cycle = ((d - START).days - epoch) / period - shift
        bright = int(np.floor(cycle + 0.5))
        dim = int(np.floor(cycle))
        if abs(cycle - bright) * period <= 3:
            aa[bright][0].append((str(d), round(float(flux), 6)))
        elif abs(cycle - (dim + 0.5)) * period <= 3:
            aa[dim][1].append((str(d), round(float(flux), 6)))
    mids = {
        str(k): median(v for d, v in l) - median(v for d, v in r)
        for k, (l, r) in aa.items()
        if l and r
    }
    return {
        "cycles": {
            str(k): {"peak": a, "half_phase": b} for k, (a, b) in sorted(aa.items())
        },
        "complete_cycles": len(mids),
        "cycle_contrasts": {k: round(v, 6) for k, v in mids.items()},
        "median_contrast": round(median(mids.values()), 6) if mids else None,
        "omit_one_range": (
            list(
                np.round(
                    [
                        min(median(v for j, v in mids.items() if j != i) for i in mids),
                        max(median(v for j, v in mids.items() if j != i) for i in mids),
                    ],
                    6,
                )
            )
            if len(mids) >= 2
            else None
        ),
    }


def rank_phases(series, period, epoch):
    base = phases(series, period, epoch, 0)
    grid = [
        phases(series, period, epoch, x / 100)["median_contrast"] for x in range(100)
    ]
    valid = [v for v in grid if v is not None]
    obs = base["median_contrast"]
    return {
        **base,
        "valid_phase_trials": len(valid),
        "shifts_at_least_as_high": (
            sum(v >= obs for v in valid) if obs is not None else None
        ),
    }


def waves(times, f):
    ang = 2 * np.pi * np.array(times) * f
    return np.column_stack([np.ones(len(ang)), np.sin(ang), np.cos(ang)])


def forecast(curves):
    red = curves["r"]
    out = {}
    grid = np.linspace(1 / 26, 1 / 14, 501)
    for yr in sorted({d.year for d in red}):
        dd = sorted(d for d in red if d.year == yr)
        if len(dd) < 40:
            continue
        cutoff = dd[len(dd) // 2 - 1]
        early = dd[: len(dd) // 2]
        later = dd[len(dd) // 2 :]
        if (
            (early[-1] - early[0]).days < 60
            or len(later) < 10
            or (later[-1] - later[0]).days < 40
        ):
            continue
        offset = dd[0]
        x = lambda dates: np.array([(d - offset).days for d in dates], dtype=float)
        t = x(early)
        y = np.array([red[d] for d in early])
        fits = []
        for f in grid:
            mat = waves(t, f)
            coef = np.linalg.lstsq(mat, y, rcond=None)[0]
            fits.append((np.sum((y - mat @ coef) ** 2), f, coef))
        _, f, c = min(fits, key=lambda v: v[0])
        result = {
            "training_through": str(cutoff),
            "n_freq": len(grid),
            "period_days": round(float(1 / f), 6),
            "scores": {},
        }
        for k, s in curves.items():
            tr = sorted(d for d in s if d.year == yr and d <= cutoff)
            te = sorted(d for d in s if d.year == yr and d > cutoff)
            item = {"training_n": len(tr), "test_n": len(te)}
            if len(tr) >= 3 and len(te) >= 3:
                v = np.array([s[d] for d in tr])
                zt = waves(x(tr), f) @ c
                scale, off = np.linalg.lstsq(
                    np.column_stack([zt, np.ones(len(zt))]), v, rcond=None
                )[0]
                pred = scale * (waves(x(te), f) @ c) + off
                yt = np.array([s[d] for d in te])
                gain = 1 - np.sum((yt - pred) ** 2) / np.sum((yt - np.mean(v)) ** 2)
                rho = spearmanr(yt, pred).statistic
                item.update(
                    {
                        "gain": round(float(gain), 5),
                        "rho": round(float(rho), 5),
                        "passes": bool(len(te) >= 10 and gain >= 0.2 and rho >= 0.5),
                    }
                )
            result["scores"][k] = item
        out[str(yr)] = result
    return out


def run(fetch):
    case, pre, _ = load(NAME)
    starts, p, epoch = ephem()
    anchor = date.fromisoformat(case["date"])
    cone, rows = obtain(
        fetch,
        None,
        neighbors=True,
        ra=case["ra"],
        dec=case["dec"],
        stem=NAME + "-asas-field",
        sources=SRC,
    )
    oid = min(cone, key=lambda v: v["separation_arcsec"])["asas_sn_id"]
    by = defaultdict(list)
    for r in rows:
        if clean(r, anchor):
            by[r["asas_sn_id"]].append(r)
    target = by[oid]
    outside = lambda rows: [r for r in rows if pilot.day_of(r["jd"]).year != 2019]
    base = median(r["flux"] for r in outside(target))
    refs = {}
    meta = []
    for star in sorted(cone, key=lambda s: s["separation_arcsec"]):
        v = star["asas_sn_id"]
        other = outside(by[v])
        if star["separation_arcsec"] < 30 or not other:
            continue
        level = median(r["flux"] for r in other)
        if not 0.5 * base <= level <= 2 * base:
            continue
        refs[v] = {
            c: median(r["flux"] for r in other if r["camera"] == c)
            for c in {r["camera"] for r in other}
        }
        meta.append(
            {"id": v, "sep_arcsec": star["separation_arcsec"], "outside_flux": level}
        )
        if len(refs) == 3:
            break
    measures = defaultdict(dict)
    for v, cams in refs.items():
        for r in by[v]:
            if r["camera"] in cams:
                measures[r["camera"], r["image_id"]][v] = r["flux"] - cams[r["camera"]]
    corr = []
    for r in target:
        control = measures[r["camera"], r["image_id"]]
        if len(control) >= 2:
            corr.append({**r, "flux": r["flux"] - median(control.values())})
    field = {str(v): daily(by[v]) for v in refs}
    curves = {
        b: {d: 3631000 * 10 ** (-0.4 * m) for d, (t, m) in pre[b].items()} for b in "rg"
    }
    curves.update({"asas": daily(target), "asas_corrected": daily(corr)})
    field.update(
        {
            "camera_" + cam: daily([r for r in target if r["camera"] == cam])
            for cam in sorted({r["camera"] for r in target})
        }
    )
    return {
        "object": NAME,
        "start": str(START),
        "stop": str(STOP),
        "period_days": p,
        "epoch_days_after_start": epoch,
        "event_starts": list(map(str, starts)),
        "first_three_fit_residuals_days": [
            round((d - START).days - epoch - i * p, 3) for i, d in enumerate(starts)
        ],
        "field": meta,
        "fixed_phase": {
            k: rank_phases(v, p, epoch) for k, v in {**curves, **field}.items()
        },
        "lag_windows_days": LAGS,
        "season_lags": {
            k: lag(v)
            for k, v in {
                **curves,
                **{k: v for k, v in field.items() if not k.startswith("camera_")},
            }.items()
        },
        "season_forecasts": forecast(curves),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "cadence.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
