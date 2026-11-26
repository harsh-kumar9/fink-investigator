"""Measure positive-to-negative centroid shifts within matched calendar years and filters."""

from collections import defaultdict
import json
import math
from statistics import median

from check import (
    HERE,
    RUN,
    compass_archive,
    gaia,
    pilot,
    positive,
    quality,
    relative_position,
)


def run():
    anchors, histories, _ = compass_archive.load(RUN / "snapshot")
    by_id = {pilot.oid(a): a for a in anchors}
    out = []
    for case in json.loads((HERE / "results.json").read_text()):
        if len(case["nearby_5arcsec"]) != 1:
            continue
        name = case["id"]
        a = by_id[name]
        source = case["nearby_5arcsec"][0]["id"]
        star = next(s for s in gaia(a, False) if s["source_id"] == source)
        groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for r in histories[name]:
            if r["i:jd"] > a["i:jd"] or not quality(r) or positive(r) is None:
                continue
            day = pilot.day_of(r["i:jd"])
            groups[(day.year, r["i:fid"])][positive(r)][day].append(
                relative_position(r, star, True)
            )
        bins = []
        for (year, fid), signs in sorted(groups.items()):
            if min(len(signs[b]) for b in (True, False)) < 3:
                continue

            def centroid(sign):
                return [
                    median(median(x[i] for x in day) for day in signs[sign].values())
                    for i in (0, 1)
                ]

            shift = [
                round(p - n, 3)
                for p, n in zip(centroid(True), centroid(False), strict=True)
            ]
            bins.append(
                {
                    "year": year,
                    "filter": fid,
                    "positive_dates": len(signs[True]),
                    "negative_dates": len(signs[False]),
                    "positive_minus_negative_xy_arcsec": shift,
                }
            )
        med = (
            [
                round(
                    median(b["positive_minus_negative_xy_arcsec"][i] for b in bins), 3
                )
                for i in (0, 1)
            ]
            if bins
            else None
        )
        out.append(
            {
                "id": name,
                "gaia": source,
                "bins": bins,
                "equal_bin_median_xy_arcsec": med,
                "length_arcsec": round(math.hypot(*med), 3) if med else None,
            }
        )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
