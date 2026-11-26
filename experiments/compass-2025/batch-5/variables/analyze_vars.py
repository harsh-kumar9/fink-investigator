"""Measure blind clock transfer, selected state windows and observed bright episodes."""

import argparse
from collections import Counter
from datetime import date
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
sys.path.insert(0, str(BATCH.parent / "batch-3/variables"))
import analyze as prior
from analyze import compass_direct, previous, pilot

NAMES = ["ZTF17aaaoxxi", "ZTF18abcckqd"]


def load(name):
    case = next(
        r
        for r in json.loads((BATCH / "results/ledger.json").read_text())["objects"]
        if r["id"] == name
    )
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / name)
    return (
        case,
        previous.nightly(rows, case["jd"] + 1e-7),
        previous.nightly(rows, math.inf),
    )


def episodes(pre):
    levels = {
        b: float(np.quantile([v for _, v in nights.values()], 0.8))
        for b, nights in pre.items()
        if nights
    }
    dd = list(pre["r"])
    if len(levels) != 2:
        return {"n_days": len(dd), "levels": levels}
    bright_g = {d for d, (_, g) in pre["g"].items() if g <= levels["g"] - 1.5}
    returns = [d for d, (_, r) in pre["r"].items() if r >= levels["r"] - 0.4]
    groups = []
    pending = []
    for d, (_, r) in pre["r"].items():
        if r <= levels["r"] - 1.5:
            pending.append(d)
        elif d in returns and pending:
            groups.append(pending)
            pending = []
    if pending:
        groups.append(pending)
    out = []
    for bb in groups:
        before = [d for d in returns if d < bb[0]]
        after = [d for d in returns if d > bb[-1]]
        gg = sorted(d for d in bright_g if any(abs((r - d).days) <= 2 for r in bb))
        l = before[-1] if before else None
        u = after[0] if after else None
        nearby = [d for d in dd if (l or bb[0]) <= d <= (u or bb[-1])]
        maxgap = max(((b - a).days for a, b in zip(nearby, nearby[1:])), default=None)
        out.append(
            {
                "bright_red_dates": list(map(str, bb)),
                "green_nearby": list(map(str, gg)),
                "brightest_red": min(pre["r"][d][1] for d in bb),
                "sampled_span_days": (bb[-1] - bb[0]).days,
                "prior_faint": str(l) if l else None,
                "next_faint": str(u) if u else None,
                "return_span_days": (u - l).days if l and u else None,
                "max_observed_date_gap": maxgap,
                "censored": l is None or u is None or bool(maxgap and maxgap > 30),
                "two_band_resolved": bool(
                    gg and any((b - a).days <= 10 for a, b in zip(bb, bb[1:]))
                ),
            }
        )
    return {"faint_q80": levels, "groups": out}


def run(name):
    case, pre, full = load(name)
    anchor = date.fromisoformat(case["date"])
    result = {
        "anchor": str(anchor),
        "bands": {},
        "episodes": episodes(pre),
        "state": prior.state(pre, anchor),
    }
    for b in "rg":
        yy = [m for _, m in pre[b].values()]
        result["bands"][b] = {
            "pre_dates": len(yy),
            "all_dates": len(full[b]),
            "nights_by_year": dict(sorted(Counter(d.year for d in pre[b]).items())),
            "q05_q50_q95_mag": (
                list(np.round(np.quantile(yy, [0.05, 0.5, 0.95]), 5)) if yy else None
            ),
        }
    if len(pre["r"]) < 12 or len(pre["g"]) < 8:
        return result
    cutoff, scanned, scores, waves = prior.get_models(pre, name)
    result.update(
        {"training_through": cutoff, "frequencies": scanned, "clocks": scores}
    )
    for key, wave in waves.items():
        for b in "rg":
            later = [d for d in full[b] if d > anchor]
            training = np.mean(
                [m for d, (_, m) in pre[b].items() if d <= date.fromisoformat(cutoff)]
            )
            origin, f, coef = wave[b]
            scores[key][b + "_post_anchor"] = previous.score(
                later, full[b], f, coef, training, origin
            )
    edge = result["state"].get("boundary_mjd")
    if edge is not None:
        resid = {}
        for b in "rg":
            origin, f, coef = waves["early_red"][b]
            t, y = prior.arrays(pre[b], list(pre[b]))
            resid[b] = dict(
                zip(pre[b], zip(t, y - previous.design(t - origin, f) @ coef))
            )
        result["state_after_clock"] = prior.state(resid, anchor, edge)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "results.json")
    pilot.json_write(p.parse_args().out, {name: run(name) for name in NAMES})
