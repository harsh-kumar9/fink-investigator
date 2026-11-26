"""Inspect the preregistered first and last held-out red pairs in the astrometric case."""

import argparse
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
    sep,
    stamps,
)


def run(fetch):
    packet = json.loads((RUN / "astrometry/packet.json").read_text())
    target = packet["case"]["id"]
    a, h, _ = compass_archive.load(RUN / "snapshot")
    anchor = next(r for r in a if pilot.oid(r) == target)
    folder = HERE / target
    if fetch:
        folder.mkdir(exist_ok=True)
    stars = gaia(anchor, fetch)
    nearby = [s for s in stars if min(sep(anchor, s, False), sep(anchor, s, True)) <= 5]
    near = []
    for s in sorted(nearby, key=lambda s: sep(anchor, s, True)):
        east, north = packet["case"]["bands"]["r"]["test_east_north_arcsec"]
        angle = None
        if s["pmra"] and s["pmdec"]:
            vx, vy = float(s["pmra"]), float(s["pmdec"])
            if math.hypot(vx, vy) and math.hypot(east, north):
                cosine = (east * vx + north * vy) / (
                    math.hypot(vx, vy) * math.hypot(east, north)
                )
                angle = round(math.degrees(math.acos(max(-1, min(1, cosine)))), 2)
        near.append(
            {
                "id": s["source_id"],
                "ra": s["ra"],
                "dec": s["dec"],
                "ref_epoch": s["ref_epoch"],
                "pmra": s["pmra"],
                "pmdec": s["pmdec"],
                "pmra_error": s["pmra_error"],
                "pmdec_error": s["pmdec_error"],
                "pm_angle_to_red_test_vector_deg": angle,
                "g_mag": s["phot_g_mean_mag"],
                "at_anchor_arcsec": round(sep(anchor, s, True), 3),
            }
        )
    look = {str(r["i:candid"]): r for r in h[target]}
    measures = []
    pairs = packet["pairs"]["held_out"]
    for i, pair in enumerate(pairs):
        save_stamp = i in {0, len(pairs) - 1}
        for sign in ("negative", "positive"):
            cids = pair[sign + "_candids"]
            dates = []
            for candid in cids:
                row = look[candid]
                if (
                    row["i:fid"] != 2
                    or not quality(row)
                    or positive(row) != (sign == "positive")
                    or str(pilot.day_of(row["i:jd"])) != pair[sign]
                ):
                    raise ValueError("Pair membership mismatch")
                if save_stamp and candid == min(cids, key=int):
                    stamps([row], row, folder / candid, fetch)
                offsets = {
                    s["source_id"]: [
                        round(k, 3) for k in relative_position(row, s, True)
                    ]
                    for s in nearby
                }
                dates.append(
                    {
                        "candid": candid,
                        "ra": row["i:ra"],
                        "dec": row["i:dec"],
                        "offset_from_gaia_xy_arcsec": offsets,
                    }
                )
            measures.append(
                {"date": pair[sign], "sign": sign, "pair": i + 1, "rows": dates}
            )
    summaries = {}
    for s in nearby:
        labels = {}
        for sign in ("positive", "negative"):
            offsets = [
                [
                    median(
                        r["offset_from_gaia_xy_arcsec"][s["source_id"]][i]
                        for r in d["rows"]
                    )
                    for i in (0, 1)
                ]
                for d in measures
                if d["sign"] == sign
            ]
            labels[sign] = {
                "median_xy_arcsec": [
                    round(median(x[i] for x in offsets), 3) for i in (0, 1)
                ],
                "median_separation_arcsec": round(
                    median(math.hypot(*x) for x in offsets), 3
                ),
            }
        summaries[s["source_id"]] = labels
    moving = pilot.stellar_matches(
        anchor, {s["source_id"]: s for s in stars if s["pmra"] and s["pmdec"]}
    )
    return {
        "id": target,
        "gaia_60arcsec_count": len(stars),
        "gaia_nearby": near,
        "moving_star_tracks": moving,
        "held_out": summaries,
        "epochs": measures,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out")
    args = p.parse_args()
    r = run(args.fetch)
    if args.out:
        pilot.json_write(HERE / args.out, r)
    else:
        print(json.dumps(r, indent=2))
