"""Check the Mira's frozen ZTF waveforms against independent ASAS-SN forced flux."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json
import math
from pathlib import Path
from statistics import median
import numpy as np
from scipy.stats import spearmanr

from asas import obtain
from analyze import HERE, EPOCH, curves, pilot, arrays, previous

NAME = "ZTF18aazutgw"


def run(fetch):
    case, pre = curves(NAME)
    cone, raw = obtain(fetch, None, ra=case["ra"], dec=case["dec"], stem="mira-asas")
    if not cone:
        return {
            "cone": [],
            "status": "Sky Patrol master catalog returned zero sources within 6 arcseconds",
        }
    bins = defaultdict(list)
    anchor = date.fromisoformat(case["date"])
    for r in raw:
        d = pilot.day_of(r["jd"])
        if (
            r["quality"] == "G"
            and date(2019, 1, 1) <= d <= anchor
            and all(math.isfinite(r[k]) for k in ("jd", "flux", "flux_err"))
            and r["flux_err"] > 0
        ):
            bins[(int(math.floor(r["jd"] - 2400000.5)) - 58000) // 10].append(
                (r["jd"], r["flux"])
            )
    base = json.loads((HERE / "results.json").read_text())[NAME]
    cutoff = date.fromisoformat(base["training_last"])
    splits = {"early": [], "later": []}
    for k, v in sorted(bins.items()):
        start = EPOCH + timedelta(days=58000 + 10 * k)
        end = start + timedelta(days=9)
        if start <= cutoff < end or end > anchor:
            continue
        splits["early" if end <= cutoff else "later"].append(
            [median(t for t, m in v), median(m for t, m in v), len(v)]
        )
    dates = [d for d in pre["r"] if d <= cutoff]
    t, y = arrays(pre["r"], dates)
    scores = {}
    for key, m in base["clocks"].items():
        f = 1 / m["period_days"]
        _, c = previous.fit(t - t[0], y, f)
        scores[key] = {}
        for split, v in splits.items():
            if len(v) < 3:
                continue
            tt, flux, count = np.array(v).T
            rho = float(spearmanr(flux, -previous.design(tt - t[0], f) @ c).statistic)
            scores[key][split] = {
                "bins": len(flux),
                "measurements": int(count.sum()),
                "spearman": round(rho, 5),
                "passes": bool(len(flux) >= 10 and rho >= 0.5),
            }
    return {
        "cone": cone,
        "raw_measurements": len(raw),
        "binned_flux_mjy": splits,
        "scores": scores,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "mira-survey.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
