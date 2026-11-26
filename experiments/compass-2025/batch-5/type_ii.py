"""Compare how much r-band flux remains in the earlier round's fixed Type II panel."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
from statistics import median

from light_curve import HERE, NAME, ZERO, onset, pilot, window_means

PREVIOUS = HERE.parent / "batch-4"


def relative(windows, name, min_baseline=4, max_fraction=None):
    baseline, early, late = (windows[k] for k in ("baseline", "early", name))
    if (
        baseline
        and baseline["n"] >= min_baseline
        and early
        and early["n"] >= 3
        and late
        and late["n"] >= 3
        and early["excess"] > 5 * early["excess_error"]
        and (max_fraction is None or baseline["error"] < max_fraction * early["excess"])
    ):
        return late["excess"] / early["excess"]
    return None


def measure(rows, zero):
    nights = [r for r in onset.group_forced(rows) if r["fid"] == 2]
    w = window_means(nights, zero)
    return dict(windows=w, ratios={k: relative(w, k) for k in ("winter", "late")})


def run(fetch):
    prior = json.loads((PREVIOUS / "results/comparison.json").read_text())
    controls = [
        r for r in prior["panels"]["2024_bright"]["controls"] if r["family"] == "II"
    ]
    folder = PREVIOUS / "comparison-sources"
    registry = gzip.decompress((folder / "tns.csv.gz").read_bytes())
    if (
        hashlib.sha256(registry).hexdigest()
        != json.loads((folder / "tns.json").read_text())["sha256"]
    ):
        raise ValueError("Registry subset checksum mismatch")
    dates = {
        r["name"]: r["discoverydate"]
        for r in csv.DictReader(registry.decode().splitlines())
    }

    def one(r):
        obj = dict(id=r["id"], name=r["tns"], audit=r["audit"])
        if r["audit"]["ok"]:
            discovery = datetime.fromisoformat(dates[r["tns"][2:]]).replace(
                tzinfo=timezone.utc
            )
            zero = discovery.timestamp() / 86400 + 40587
            rows = onset.forced_inputs(
                fetch, r["id"], HERE / "comparison-sources" / r["id"]
            )
            obj.update(clock_mjd=zero, **measure(rows, zero))
        return obj

    with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
        panel = list(pool.map(one, controls))
    target = measure(onset.forced_inputs(False, NAME, HERE / "sn-sources"), ZERO)
    summary = {}
    for k in target["ratios"]:
        ratios = [
            r["ratios"][k] for r in panel if r.get("ratios", {}).get(k) is not None
        ]
        summary[k] = dict(
            eligible=len(ratios),
            median=median(ratios) if ratios else None,
            below_target=sum(v < target["ratios"][k] for v in ratios),
            above_target=sum(v > target["ratios"][k] for v in ratios),
        )
    sensitivity = {}
    for k in target["ratios"]:
        sensitivity[k] = {}
        for n in (1, 2):
            ratios = {
                r["name"]: relative(r["windows"], k, n, 0.05)
                for r in panel
                if "windows" in r
            }
            sensitivity[k][str(n)] = {
                key: v for key, v in ratios.items() if v is not None
            }
    return dict(
        target=target,
        controls=panel,
        comparisons=summary,
        limited_baselines=sensitivity,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "results/type-ii.json")
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch))
