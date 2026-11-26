"""Evaluate alternating strengths on all other years at the discovery-season cut."""

import argparse
from collections import defaultdict
from datetime import date
import json
from pathlib import Path
from statistics import median

import numpy as np

from analyze_vars import HERE, BATCH, compass_direct, load, pilot
from cadence import NAME, START, STOP, ephem


def group_peaks():
    base = json.loads((HERE / "duration.json").read_text())[NAME]["tests"]["primary"][
        "groups"
    ]
    case, pre, _ = load(NAME)
    _, raw, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    errs = defaultdict(list)
    for r in raw:
        if int(r["catflags"]) == 0 and 0 < float(r["magerr"]) <= 0.2:
            errs[
                r["filtercode"][-1], str(pilot.day_of(float(r["mjd"]) + 2400000.5))
            ].append(float(r["magerr"]))
    out = []
    for gr in base:
        dd = [date.fromisoformat(d) for d in gr["dates"]]
        item = {
            k: gr[k]
            for k in (
                "first",
                "last",
                "previous_faint",
                "next_faint",
                "lower_days",
                "upper_days",
                "kind",
            )
        }
        item["n_bright_dates"] = len(dd)
        for b in "rg":
            zz = sorted((pre[b][d][1], d) for d in dd if d in pre[b])
            item[b] = (
                None
                if not zz
                else {
                    "date": str(zz[0][1]),
                    "mag": round(zz[0][0], 5),
                    "max_error": max(errs[b, str(zz[0][1])]),
                    "n": len(zz),
                }
            )
        out.append(item)
    return out


def get_training(groups):
    starts, period, epoch = ephem()
    weak = [next(g for g in groups if g["first"] == str(d)) for d in starts]
    strong = []
    for k in range(4):
        mid = epoch + (k + 0.5) * period
        candidates = [
            g
            for g in groups
            if g["r"]
            and abs((date.fromisoformat(g["r"]["date"]) - START).days - mid) <= 4
        ]
        if len(candidates) != 1:
            raise ValueError("Intervening group is ambiguous")
        strong.extend(candidates)
    threshold = {
        b: (
            median(v[b]["mag"] for v in weak[:3] if v[b])
            + median(v[b]["mag"] for v in strong[:2] if v[b])
        )
        / 2
        for b in "rg"
    }
    return threshold, weak, strong


def eligible(groups, threshold):
    pairs = []
    rr = groups
    for i in range(len(rr) - 1):
        a, b = rr[i : i + 2]
        if not (a["r"] and b["r"] and a["next_faint"] and a["next_faint"] < b["first"]):
            continue
        d, e = map(date.fromisoformat, (a["r"]["date"], b["r"]["date"]))
        if d.year != e.year or not 5 <= (e - d).days <= 15:
            continue
        pairs.append([i, i + 1])
    return np.array(pairs, dtype=int).reshape(-1, 2)


def check(groups, threshold, ignore_discovery):
    g = [v for v in groups if not ignore_discovery or not v["first"].startswith("2019")]
    pairs = eligible(g, threshold)
    strata = defaultdict(list)
    for i, v in enumerate(g):
        key = (v["first"][:4], min(v["r"]["n"], 2) if v["r"] else 0, bool(v["g"]))
        strata[key].append(i)
    strata = [np.array(v, dtype=int) for v in strata.values() if len(v) > 1]
    labels = np.array(
        [[int(v[b]["mag"] <= threshold[b]) if v[b] else -1 for b in "rg"] for v in g],
        int,
    )
    n_red = np.array([v["r"]["n"] if v["r"] else 0 for v in g])
    selections = {
        "r": (0, pairs),
        "g": (1, pairs[np.all(labels[pairs, 1] >= 0, axis=1)]),
        "r_two_nights": (0, pairs[np.all(n_red[pairs] >= 2, axis=1)]),
    }

    def mixed(ax, sub, arr):
        return int(np.sum(arr[sub[:, 0], ax] != arr[sub[:, 1], ax]))

    obs = {k: mixed(b, s, labels) for k, (b, s) in selections.items()}
    rng = np.random.default_rng(20250921)
    trials = {k: [] for k in selections}
    for _ in range(5000):
        perm = np.arange(len(g))
        for ix in strata:
            perm[ix] = rng.permutation(ix)
        v = labels[perm]
        for k, (b, s) in selections.items():
            trials[k].append(mixed(b, s, v))
    scores = {}
    for k, (b, s) in selections.items():
        a = np.array(trials[k])
        scores[k] = {
            "eligible_links": len(s),
            "participating_events": len(set(s.ravel())),
            "mixed_links": obs[k],
            "permutation_mean": round(float(np.mean(a)), 3),
            "permutation_greater_or_equal": int(np.sum(a >= obs[k])),
            "trials": len(a),
            "plus_one_tail": (
                round(float((1 + np.sum(a >= obs[k])) / (len(a) + 1)), 5)
                if len(s)
                else None
            ),
            "by_year": {},
        }
        for yr in sorted({v["first"][:4] for v in g}):
            ss = np.array(
                [pair for pair in s if g[pair[0]]["first"].startswith(yr)], dtype=int
            ).reshape(-1, 2)
            scores[k]["by_year"][yr] = {"links": len(ss), "mixed": mixed(b, ss, labels)}
    return {"n_groups": len(g), "chronological_groups": g, "scores": scores}


def run():
    gr = group_peaks()
    cut, weak, strong = get_training(gr)
    return {
        "cuts": cut,
        "2019_short": weak,
        "2019_intervening": strong,
        "other_years": check(gr, cut, True),
        "all_2019_groups": [
            g for g in gr if "2019-07-15" <= g["first"] <= "2019-10-31"
        ],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "alternating.json")
    pilot.json_write(p.parse_args().out, run())
