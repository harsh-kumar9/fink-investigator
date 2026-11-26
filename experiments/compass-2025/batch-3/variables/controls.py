"""Challenge two-week recurrence with same-image field stars and lunar/camera terms."""

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import numpy as np
from asas import obtain, RA, DEC
from analyze import BATCH, HERE, pilot
from interval import NAME
from recurrence import accepted, evaluate


def compact(raw):
    return {
        y: {k: v for k, v in r.items() if k != "per_anchor"}
        for y, r in evaluate(raw).items()
    }


def run(fetch):
    cone, raw = obtain(fetch, None, neighbors=True)
    target = min(cone, key=lambda r: r["separation_arcsec"])
    gaia = BATCH / "images" / NAME / "gaia.csv.gz"
    data = gzip.decompress(gaia.read_bytes())
    receipt = json.loads(gaia.with_name("gaia.json").read_text())
    assert hashlib.sha256(data).hexdigest() == receipt["sha256_uncompressed"]
    stars = list(csv.DictReader(data.decode().splitlines()))

    def match(r):
        sep = lambda s: math.hypot(
            *pilot.offsets(float(s["ra"]), float(s["dec"]), r["ra_deg"], r["dec_deg"])
        )
        star = min(stars, key=sep)
        return star, sep(star)

    g = float(match(target)[0]["phot_g_mean_mag"])
    helpers = []
    for r in sorted(cone, key=lambda p: p["separation_arcsec"]):
        star, d = match(r)
        if (
            r["separation_arcsec"] >= 30
            and d <= 2
            and abs(float(star["phot_g_mean_mag"]) - g) <= 1.5
        ):
            helpers.append(
                {
                    "asas_sn_id": r["asas_sn_id"],
                    "gaia_source": star["source_id"],
                    "gaia_g": float(star["phot_g_mean_mag"]),
                    "offset_from_target_arcsec": r["separation_arcsec"],
                    "match_arcsec": round(d, 3),
                }
            )
        if len(helpers) == 3:
            break
    raw = accepted(raw)
    near = {
        r["asas_sn_id"]: [v for v in raw if v["asas_sn_id"] == r["asas_sn_id"]]
        for r in helpers
    }
    target_rows = [v for v in raw if v["asas_sn_id"] == target["asas_sn_id"]]
    fields = defaultdict(dict)
    for source, rows in near.items():
        if not rows:
            continue
        level = median(r["flux"] for r in rows)
        for r in rows:
            key = (r["image_id"], r["camera"])
            fields[key][source] = r["flux"] - level
    corrected = []
    common = []
    for r in target_rows:
        matches = fields[(r["image_id"], r["camera"])]
        if len(matches) >= 2:
            common.append(r)
            corrected.append({**r, "flux": r["flux"] - median(matches.values())})
    seasonal = []
    years = np.array([pilot.day_of(r["jd"]).year for r in target_rows])
    cams = sorted({r["camera"] for r in target_rows})
    phase = (np.array([r["jd"] for r in target_rows]) - 2451545) * 2 * np.pi / 29.53059
    mat = np.column_stack(
        [
            *[
                np.array([r["camera"] == c for r in target_rows], dtype=float)
                for c in cams
            ],
            np.sin(phase),
            np.cos(phase),
            np.sin(2 * phase),
            np.cos(2 * phase),
        ]
    )
    yy = np.array([r["flux"] for r in target_rows])
    for year in sorted(set(years)):
        train = years != year
        coeff = np.linalg.lstsq(mat[train], yy[train], rcond=None)[0]
        for idx in np.where(~train)[0]:
            seasonal.append(
                {**target_rows[idx], "flux": float(yy[idx] - mat[idx] @ coeff)}
            )
    return {
        "comparison_sources": [
            {**h, "usable_measurements": len(near[h["asas_sn_id"]])} for h in helpers
        ],
        "lunar_camera": compact(seasonal),
        "common_target": compact(common),
        "field_corrected_target": compact(corrected),
        "individual_fields": {str(s): compact(rows) for s, rows in near.items()},
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "controls.json")
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch))
