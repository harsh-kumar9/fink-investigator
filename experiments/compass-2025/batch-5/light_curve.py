"""Measure the observed evolution of SN 2025xza from archived signed fluxes."""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from math import hypot, sqrt
from pathlib import Path
from statistics import mean, stdev
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "batch-4"))
import onset
import pilot

NAME = "ZTF25absffwy"
FIRST = datetime(2025, 9, 18, 4, 43, 23, 520000, tzinfo=timezone.utc)
ZERO = FIRST.timestamp() / 86400 + 40587
WINDOWS = {
    "baseline": (-120, -7),
    "early": (5, 25),
    "middle": (30, 60),
    "winter": (75, 110),
    "late": (115, 150),
    "spring": (155, 200),
}


def window_means(nights, zero=ZERO):
    result = {}
    for name, (start, stop) in WINDOWS.items():
        rows = [r for r in nights if start <= r["mjd"] - zero < stop + 1]
        if not rows:
            result[name] = None
            continue
        n = len(rows)
        err = sqrt(sum(r["error"] ** 2 for r in rows)) / n
        if n > 1:
            err = max(err, stdev(r["flux"] for r in rows) / sqrt(n))
        result[name] = dict(
            n=n,
            first=rows[0]["date"],
            last=rows[-1]["date"],
            mean=mean(r["flux"] for r in rows),
            error=err,
        )
    baseline = result["baseline"]
    if baseline:
        for name, w in result.items():
            if w and name != "baseline":
                w.update(
                    excess=w["mean"] - baseline["mean"],
                    excess_error=hypot(w["error"], baseline["error"]),
                )
    return result


def run(fetch):
    forced, alerts, limits = onset.inputs(fetch, NAME, HERE / "sn-sources")
    nights = onset.group_forced(forced)
    by_reference = defaultdict(list)
    for r in forced:
        by_reference[(r["field"], r["rcid"], r["rfid"], r["fid"])].append(r)
    det = onset.difference_detections(alerts)
    _, histories, _ = onset.compass_archive.load(HERE.parent / "snapshot")
    by_candid = {str(r["candid"]): r for r in alerts}
    for r in histories[NAME]:
        d = by_candid[str(r["i:candid"])]
        if (
            d["fid"] != r["i:fid"]
            or abs(d["mjd"] + 2400000.5 - r["i:jd"]) > 1 / 86400
            or abs(d["magpsf"] - r["i:magpsf"]) > 0.001
        ):
            raise ValueError("Mismatch between Fink and ALeRCE alerts")
    pairs = defaultdict(dict)
    for r in nights:
        if r["snr"] >= 5:
            pairs[r["day"]][r["fid"]] = r
    colors = []
    for day, v in sorted(pairs.items()):
        if 1 in v and 2 in v and abs(v[1]["mjd"] - v[2]["mjd"]) <= 0.5:
            colors.append(
                dict(
                    date=v[1]["date"],
                    phase=day - int(ZERO),
                    color=v[1]["mag"] - v[2]["mag"],
                    error=hypot(v[1]["magerr"], v[2]["magerr"]),
                    hours=24 * abs(v[1]["mjd"] - v[2]["mjd"]),
                )
            )
    return dict(
        clock_mjd=ZERO,
        forced_rows=len(forced),
        used=sum(len(r["pids"]) for r in nights),
        fink=len(histories[NAME]),
        alert_exposures=len(det),
        limit_rows=len(limits),
        windows={
            str(fid): window_means([r for r in nights if r["fid"] == fid])
            for fid in [1, 2]
        },
        detector_windows={
            "/".join(key): window_means(onset.group_forced(rows))
            for key, rows in sorted(by_reference.items())
        },
        nights=nights,
        colors=colors,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--out", type=Path, default=HERE / "results/light-curve.json")
    args = parser.parse_args()
    pilot.json_write(args.out, run(args.fetch))
