"""Compare the borderline VSX counterpart to independent epoch-propagated positions."""
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]))
import pilot


def separation(first, second):
    return math.hypot(*pilot.offsets(*first, *second))


def julian_year(jd):
    return 2000 + (jd - 2451545) / 365.25


def propagate(source, epoch):
    dec = float(source["dec"])
    dt = epoch - float(source["ref_epoch"])
    factor = dt / 3600000
    return (float(source["ra"]) + factor * float(source["pmra"]) / math.cos(math.radians(dec)),
            dec + factor * float(source["pmdec"]))


def main():
    for meta in json.loads((HERE / "sources.json").read_text())["files"]:
        if hashlib.sha256((HERE / meta["file"]).read_bytes()).hexdigest() != meta["sha256"]:
            raise ValueError(f"catalog checksum changed for {meta['file']}")
    folder = HERE.parent
    case = next(r for r in json.loads((folder / "results/analysis.json").read_text())["objects"]
                if r["id"] == "ZTF18abdjpbo")
    anchor = (case["ra"], case["dec"])
    vsx = case["vsx_matches"][0]
    variable = (vsx["ra"], vsx["dec"])
    with (HERE / "gaia.csv").open() as f:
        stars = list(csv.DictReader(f))
    epoch = julian_year(case["jd"])
    motions = [s for s in stars if s["pmra"] and s["pmdec"]]
    nearby = sorted(motions, key=lambda s: separation(anchor, propagate(s, epoch)))
    star = nearby[0]
    with (HERE / "ps1.csv").open() as f:
        catalog = list(csv.DictReader(f))
    ps1 = min(catalog, key=lambda s: separation(anchor, (float(s["raMean"]), float(s["decMean"]))))
    _, raw = pilot.snapshot(folder / "snapshot")
    hist = pilot.histories([r for k, rows in raw.items() if k.startswith("history-") for r in rows])
    by_day = defaultdict(list)
    for row in hist[case["id"]]:
        if (row["i:jd"] <= case["jd"] + 1e-7 and (row["i:drb"] or 0) >= 0.8
                and 0 < (row["i:sigmapsf"] or 0) <= 0.2 and str(row["i:isdiffpos"]).lower() in {"t","1","true"}):
            by_day[pilot.day_of(row["i:jd"])].append((row["i:jd"], row["i:ra"], row["i:dec"]))
    nights = [tuple(map(median, zip(*rows))) for rows in by_day.values()]
    near = [separation((ra, dec), propagate(star, julian_year(jd))) for jd, ra, dec in nights]
    fixed = [separation((ra, dec), variable) for jd, ra, dec in nights]
    return {"object_id": case["id"], "candid": case["candid"], "gaia_id": star["source_id"],
            "gaia_sources_in_cone": len(stars), "ps1_objects_in_cone": len(catalog),
            "nearest_ps1_id": ps1["objID"], "quality_positive_nights": len(nights),
            "anchor_separation_gaia_arcsec": round(separation(anchor, propagate(star, epoch)), 3),
            "next_gaia_separation_at_anchor_arcsec": round(separation(anchor, propagate(nearby[1], epoch)), 3),
            "median_nightly_separation_gaia_arcsec": round(median(near), 3),
            "median_nightly_separation_vsx_arcsec": round(median(fixed), 3),
            "vsx_separation_gaia_j2000_arcsec": round(separation(variable, propagate(star, 2000)), 3)}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
