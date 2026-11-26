"""Measure sampled eruption durations and check the same dates with ASAS-SN."""

import argparse
from collections import defaultdict
from datetime import date
import math
from pathlib import Path
from statistics import median

import numpy as np
from scipy.stats import spearmanr

from analyze_vars import HERE, NAMES, load, pilot
from asas import obtain

SRC = HERE / "sources"


def states(data, levels, delta):
    out = {}
    for d in sorted(data["r"].keys() | data["g"].keys()):
        vals = {b: m[d][1] for b, m in data.items() if d in m}
        bright = sorted(b for b, mag in vals.items() if mag <= levels[b] - delta)
        faint = sorted(b for b, mag in vals.items() if mag >= levels[b] - 0.4)
        kind = (
            "conflict"
            if bright and faint
            else "bright" if bright else "faint" if faint else "middle"
        )
        out[d] = (kind, bright, faint, vals)
    return out


def classify(points, gap):
    faint = [d for d, v in points.items() if v[0] == "faint"]
    groups = []
    active = []
    for d, (kind, br, f, values) in points.items():
        if kind == "faint":
            if active:
                groups.append(active)
                active = []
        if kind == "bright":
            if active and (d - active[-1]).days > gap:
                groups.append(active)
                active = []
            active.append(d)
    if active:
        groups.append(active)
    out = []
    for dd in groups:
        before = [d for d in faint if d < dd[0]]
        after = [d for d in faint if d > dd[-1]]
        l = before[-1] if before else None
        u = after[0] if after else None
        span = (dd[-1] - dd[0]).days
        bound = (u - l).days if l and u else None
        out.append(
            {
                "dates": list(map(str, dd)),
                "bright_bands": ["".join(points[d][1]) for d in dd],
                "first": str(dd[0]),
                "last": str(dd[-1]),
                "previous_faint": str(l) if l else None,
                "next_faint": str(u) if u else None,
                "lower_days": span,
                "upper_days": bound,
                "kind": (
                    "long"
                    if len(dd) >= 4 and span >= 10
                    else "short" if bound is not None and bound <= 7 else "unresolved"
                ),
            }
        )
    return out


def archive(case, fetch):
    SRC.mkdir(exist_ok=True)
    cone, rows = obtain(
        fetch,
        None,
        ra=case["ra"],
        dec=case["dec"],
        stem=case["id"] + "-asas",
        sources=SRC,
    )
    if not cone or min(v["separation_arcsec"] for v in cone) > 3:
        return cone, {}
    oid = min(cone, key=lambda v: v["separation_arcsec"])["asas_sn_id"]
    daily = defaultdict(list)
    for r in rows:
        if (
            r["asas_sn_id"] != oid
            or r["quality"] != "G"
            or not all(math.isfinite(float(r[k])) for k in ("jd", "flux", "flux_err"))
            or r["flux_err"] <= 0
        ):
            continue
        d = pilot.day_of(r["jd"])
        if date(2019, 1, 1) <= d <= date.fromisoformat(case["date"]):
            daily[d].append(r)
    return cone, {
        d: {
            "flux": median(x["flux"] for x in vv),
            "max_error": max(x["flux_err"] for x in vv),
            "cameras": sorted(set(x["camera"] for x in vv)),
        }
        for d, vv in daily.items()
    }


def external(gr, days):
    target = [date.fromisoformat(v) for v in gr["dates"]]
    brackets = [
        date.fromisoformat(gr[k]) for k in ("previous_faint", "next_faint") if gr[k]
    ]
    near = lambda d, vals: any(abs((d - v).days) <= 1 for v in vals)
    br = sorted(d for d in days if near(d, target) and not near(d, brackets))
    fa = sorted(d for d in days if near(d, brackets) and not near(d, target))
    flux = lambda dd: median([days[d]["flux"] for d in dd])
    return {
        "bright_dates": {str(d): days[d] for d in br},
        "bracket_dates": {str(d): days[d] for d in fa},
        "bright_minus_bracket_mjy": (
            round(flux(br) - flux(fa), 6) if br and fa else None
        ),
    }


def run(name, fetch):
    case, pre, _ = load(name)
    levels = {
        b: float(np.quantile([m for t, m in vv.values()], 0.8)) for b, vv in pre.items()
    }
    cone, days = archive(case, fetch)
    common = sorted(pre["g"].keys() & days.keys())
    rho = (
        spearmanr(
            [-pre["g"][d][1] for d in common], [days[d]["flux"] for d in common]
        ).statistic
        if len(common) >= 3
        else None
    )
    out = {
        "levels": levels,
        "asas_cone": cone,
        "asas_nights": len(days),
        "common_g_nights": len(common),
        "common_g_rho": round(float(rho), 5) if rho is not None else None,
        "tests": {},
    }
    for label, cut, gap in [
        ("primary", 1.5, 5),
        ("shallow", 1.25, 5),
        ("deep", 1.75, 5),
        ("tight_gap", 1.5, 4),
    ]:
        points = states(pre, levels, cut)
        groups = classify(points, gap)
        examples = [v for v in groups if v["kind"] != "unresolved"]
        test = {
            "bright_cut": cut,
            "gap": gap,
            "total_groups": len(groups),
            "long": sum(g["kind"] == "long" for g in groups),
            "short": sum(g["kind"] == "short" for g in groups),
            "conflicted_dates": [
                str(d) for d, v in points.items() if v[0] == "conflict"
            ],
        }
        if label == "primary":
            test["groups"] = groups
            for g in examples:
                g["asas"] = external(g, days)
        else:
            test["examples"] = [
                {k: g[k] for k in ("first", "last", "lower_days", "upper_days", "kind")}
                for g in examples
            ]
        out["tests"][label] = test
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "duration.json")
    args = p.parse_args()
    pilot.json_write(args.out, {n: run(n, args.fetch) for n in NAMES})
