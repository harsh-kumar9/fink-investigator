"""Calibrate the frozen original-image sample with isolated DR2 Pan-STARRS stars."""

import argparse
import math
from statistics import median
from urllib.parse import urlencode
import io
from pathlib import Path

from season import (
    HERE,
    POS,
    select,
    I_WINDOWS,
    WINDOWS,
    FOLDER,
    late,
    pilot,
    nightly_contrasts,
)
from compass_direct import capture
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

FIELDS = [
    "objID",
    "raMean",
    "decMean",
    "nDetections",
    *[b + k for b in "ri" for k in ("MeanPSFMag", "MeanPSFMagErr", "MeanKronMag")],
]


def catalog(fetch):
    url = (
        "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?"
        + urlencode(
            dict(
                ra=POS[0],
                dec=POS[1],
                radius=0.034,
                pagesize=3000,
                columns="[" + ",".join(FIELDS) + "]",
            )
        )
    )
    rows = capture(HERE / "standards-sources", "ps1", url, set(FIELDS), fetch)
    if len(rows) >= 3000 or len({r["objID"] for r in rows}) != len(rows):
        raise ValueError("Invalid or truncated PS1 cone")
    return rows


def references(rows, band):
    stars = []
    for row in rows:
        if not row[band + "MeanPSFMag"]:
            continue
        mag = float(row[band + "MeanPSFMag"])
        if mag > 0:
            pos = (float(row["raMean"]), float(row["decMean"]))
            stars.append((row, pos, mag))
    out = []
    for r, pos, mag in stars:
        err = float(r[band + "MeanPSFMagErr"] or -1)
        kron = float(r[band + "MeanKronMag"] or -999)
        if not (
            15 <= mag <= 19
            and 0 < err <= 0.05
            and abs(mag - kron) <= 0.10
            and int(r["nDetections"]) >= 5
            and 35 <= math.hypot(*pilot.offsets(*POS, *pos)) <= 100
        ):
            continue
        if any(
            q["objID"] != r["objID"]
            and m <= mag + 2
            and math.hypot(*pilot.offsets(*pos, *other)) < 8
            for q, other, m in stars
        ):
            continue
        out.append((r["objID"], pos, mag))
    return sorted(out)


def measure(item, folder, refs, controls, stars, wide_counts=None):
    win, date, row = item
    raw = late.image(row, False, POS, folder)
    result = dict(window=win, date=date, pid=row["pid"], seeing=float(row["seeing"]))
    if raw is None:
        return {**result, "status": "image unavailable"}
    with fits.open(io.BytesIO(raw)) as hdus:
        h = hdus[0].header
        if (
            int(h["EXPID"]) != int(row["expid"])
            or h["FILTER"] != "ZTF_" + row["filtercode"][1:]
        ):
            raise ValueError("Exposure identity changed")
        image = hdus[0].data.astype(float)
        wcs = WCS(h)

        def counts(pos):
            return late.phot.aperture(image, *wcs.all_world2pix([pos], 0)[0])

        bx, by = wcs.all_world2pix([POS], 0)[0]
        xy = wcs.all_world2pix(
            [
                (float(r["raMean"]), float(r["decMean"]))
                for r in stars
                if int(r["nDetections"]) >= 5
            ],
            0,
        )
        blanks = []
        for x in range(14, image.shape[1] - 13, 10):
            for y in range(14, image.shape[0] - 13, 10):
                if (
                    math.hypot(x - bx, y - by) < 38
                    or min(np.hypot(xy[:, 0] - x, xy[:, 1] - y)) < 9
                ):
                    continue
                c = late.phot.aperture(image, x, y)
                if c is not None:
                    blanks.append(c)
        if len(blanks) < 12:
            return {**result, "status": "insufficient sky", "blanks": len(blanks)}
        sky = median(blanks)
        sigma = late.phot.robust_scatter(blanks)
        zero = {}
        for oid, pos, mag in refs:
            c = counts(pos) if wide_counts is None else wide_counts.get(oid)
            if c is not None and c - sky > 5 * sigma:
                zero[oid] = mag + 2.5 * math.log10(c - sky)
        if len(zero) < 3:
            return {**result, "status": "insufficient standards", "stars": len(zero)}
        zp = median(zero.values())
        scale = 10 ** (0.4 * (23.9 - zp))
        signal = counts(POS) - sky
        source = signal * scale
        estimates = {
            oid: v - median(b for key, b in zero.items() if key != oid)
            for oid, v in zero.items()
        }
        err = math.hypot(
            math.sqrt(sigma**2 + max(0.0, signal) / float(h["GAIN"])) * scale,
            source * math.log(10) * 0.4 * late.phot.robust_scatter(list(zero.values())),
        )
        ctrl = {
            key: (counts(pos) - sky) * scale
            for key, pos in controls.items()
            if counts(pos) is not None
        }
        return {
            **result,
            "status": "measured",
            "flux_ujy": source,
            "error_ujy": err,
            "blank_flux_ujy": ctrl,
            "control": None,
            "zero_point": zp,
            "stars": len(zero),
            "excluded_star_residual_mag": estimates,
        }


def run(fetch, wide=False):
    cat = catalog(fetch)
    out = {}
    for band, windows, folder in [
        ("r", WINDOWS, FOLDER),
        ("i", I_WINDOWS, HERE / "season-i-sources"),
    ]:
        inv, _ = late.tables(False, POS, folder, {"r": 2, "i": 3}[band])
        _, sample = select(inv, band, windows)
        controls = {
            f"{x},{y}": (
                POS[0] + x / 3600 / math.cos(math.radians(POS[1])),
                POS[1] + y / 3600,
            )
            for d in [12, 24]
            for x, y in [(0, d), (0, -d), (d, 0), (-d, 0)]
        }
        refs = references(cat, band)
        if wide:
            from field_stars import (
                catalog as field_catalog,
                standards as field_standards,
                tile_counts,
            )
            from concurrent.futures import ThreadPoolExecutor

            stars = field_standards(field_catalog(band, False), band)
            refs = [
                (r["oid"], (float(r["ra"]), float(r["dec"])), float(r["meanmag"]))
                for r in stars
            ]

            def one(item):
                counts = tile_counts(item[2], refs, folder, fetch)
                return measure(item, folder, refs, controls, cat, counts)

            with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
                rows = list(pool.map(one, sample))
        else:
            rows = [measure(item, folder, refs, controls, cat) for item in sample]
        valid = [r for r in rows if r["status"] == "measured"]
        base = [r for r in valid if r["window"] == "before"]
        comparisons = {}
        for w in windows:
            if w == "before":
                continue
            post = [r for r in valid if r["window"] == w]

            def matched(a, b):
                return [
                    r
                    for r in a
                    if any(abs(r["seeing"] - s["seeing"]) <= 0.4 for s in b)
                ]

            pre, after = matched(base, post), matched(post, base)
            common = (
                set.intersection(
                    *(set(r["excluded_star_residual_mag"]) for r in pre + after)
                )
                if pre and after
                else set()
            )
            checks = {
                oid: median(r["excluded_star_residual_mag"][oid] for r in after)
                - median(r["excluded_star_residual_mag"][oid] for r in pre)
                for oid in sorted(common)
            }
            comparisons[w] = dict(
                all=late.contrast(base, post),
                matched=late.contrast(pre, after),
                nightly=nightly_contrasts(base, post),
                star_shifts=checks,
            )
        out[band] = dict(references=refs, measurements=rows, comparisons=comparisons)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--wide", action="store_true")
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    name = "field-calibrated" if args.wide else "ps1-calibrated"
    pilot.json_write(
        args.out or HERE / f"results/{name}.json", run(args.fetch, args.wide)
    )
