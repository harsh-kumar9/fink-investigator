"""Reduce independent TESS pixels and check if the dwarf nova varied in the central aperture."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from zipfile import ZipFile

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.stats import spearmanr

from analyze import HERE, curves, pilot

NAME = "ZTF17aabshbq"
RA, DEC = 99.526908, 45.2491446
BASE = "https://mast.stsci.edu/tesscut/api/v0.1/"
SECTORS = BASE + "sector?" + urlencode({"ra": RA, "dec": DEC, "radius": "0s"})
SOURCE = HERE / "sources"
APERTURES = {
    "center": (3, 3),
    "lower_left": (0, 0),
    "lower_right": (0, 6),
    "upper_left": (6, 0),
    "upper_right": (6, 6),
}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def sector_url(sector):
    return (
        BASE
        + "astrocut?"
        + urlencode(
            {"ra": RA, "dec": DEC, "x": 9, "y": 9, "units": "px", "sector": sector}
        )
    )


def sectors(fetch, local):
    raw_path = SOURCE / "tess-sectors.json"
    receipt = SOURCE / "tess-sectors-receipt.json"
    if not receipt.exists() and (fetch or local):
        raw = (
            (local / "fink-mgab-v243-sectors.json").read_bytes()
            if local
            else pilot.curl(SECTORS, [])
        )
        json.loads(raw)["results"]
        raw_path.write_bytes(raw)
        pilot.json_write(
            receipt,
            {
                "url": SECTORS,
                "sha256": digest(raw),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    r = json.loads(receipt.read_text())
    raw = raw_path.read_bytes()
    assert r["url"] == SECTORS and r["sha256"] == digest(raw)
    values = sorted(
        (int(k["sector"]), int(k["camera"]), int(k["ccd"]))
        for k in json.loads(raw)["results"]
        if int(k["sector"]) <= 93
    )
    assert len(set(k[0] for k in values)) == len(values)
    return values[:4]


def store_array(path, arr):
    stream = io.StringIO()
    np.savetxt(stream, arr, fmt="%.8f", delimiter=",")
    raw = stream.getvalue().encode()
    path.write_bytes(gzip.compress(raw, mtime=0))
    return digest(raw)


def reduction(sector, camera, ccd, fetch, local):
    receipt = SOURCE / f"tess-{sector}.json"
    images = SOURCE / f"tess-{sector}-six-hour.csv.gz"
    apertures = SOURCE / f"tess-{sector}-apertures.csv.gz"
    url = sector_url(sector)
    if not receipt.exists() and (fetch or local):
        zipped = (
            (local / f"fink-mgab-v243-sector{sector}.zip").read_bytes()
            if local
            else pilot.curl(url, [])
        )
        with ZipFile(io.BytesIO(zipped)) as z:
            (name,) = z.namelist()
            raw = z.read(name)
        with fits.open(io.BytesIO(raw)) as h:
            assert all(
                int(h[0].header[k]) == v
                for k, v in [("SECTOR", sector), ("CAMERA", camera), ("CCD", ccd)]
            )
            assert h[1].header["TIMESYS"] == "TDB"
            x, y = WCS(h[2].header).world_to_pixel(SkyCoord(RA, DEC, unit="deg"))
            assert max(abs(x - 4), abs(y - 4)) < 1.5, "Target outside central aperture"
            data = h[1].data
            cube = np.asarray(data["FLUX"], dtype=float)
            assert cube.shape[1:] == (9, 9)
            good = (
                (data["QUALITY"] == 0)
                & np.isfinite(data["TIME"])
                & np.all(np.isfinite(cube), axis=(1, 2))
            )
            flux = cube[good]
            border = np.ones((9, 9), dtype=bool)
            border[1:-1, 1:-1] = False
            flux -= np.median(flux[:, border], axis=1)[:, None, None]
            time = data["TIME"][good] + h[1].header["BJDREFI"] + h[1].header["BJDREFF"]
            six = []
            indices = np.floor(time * 4).astype(int)
            for k in np.unique(indices):
                ok = indices == k
                if ok.sum() >= 5:
                    six.append(
                        [np.median(time[ok]), *np.median(flux[ok], axis=0).ravel()]
                    )
            a = np.column_stack(
                [
                    time,
                    *[
                        flux[:, r : r + 3, c : c + 3].sum(axis=(1, 2))
                        for r, c in APERTURES.values()
                    ],
                ]
            )
            meta = {
                "sector": sector,
                "camera": camera,
                "ccd": ccd,
                "url": url,
                "fits_name": name,
                "fits_sha256": digest(raw),
                "target_xy": [float(x), float(y)],
                "total_cadences": len(data),
                "retained": int(good.sum()),
                "six_hour_sha256": store_array(images, six),
                "apertures_sha256": store_array(apertures, a),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            }
            pilot.json_write(receipt, meta)
    m = json.loads(receipt.read_text())
    assert (m["sector"], m["camera"], m["ccd"], m["url"]) == (sector, camera, ccd, url)
    data = gzip.decompress(images.read_bytes())
    ap = gzip.decompress(apertures.read_bytes())
    assert digest(data) == m["six_hour_sha256"] and digest(ap) == m["apertures_sha256"]
    return m, np.loadtxt(io.BytesIO(data), delimiter=",")


def compare(meta, array, pre):
    t = array[:, 0]
    pixels = array[:, 1:].reshape(-1, 9, 9)
    peak = None
    out = {}
    for b in "rg":
        matches = []
        for d, (t0, y) in pre[b].items():
            idx = int(np.argmin(abs(t - t0)))
            if abs(t[idx] - t0) <= 0.25:
                matches.append((d, idx, y))
        metrics = {"matched_nights": len(matches)}
        if len(matches) >= 3:
            ii = [k for d, k, y in matches]
            mag = np.array([y for d, k, y in matches])
            scores = {}
            for key, (r, c) in APERTURES.items():
                counts = pixels[ii, r : r + 3, c : c + 3].sum(axis=(1, 2))
                scores[key] = round(float(spearmanr(counts, -mag).statistic), 5)
            order = np.argsort(mag)
            n = len(order) // 3
            delta = np.median(pixels[np.array(ii)[order[:n]]], axis=0) - np.median(
                pixels[np.array(ii)[order[-n:]]], axis=0
            )
            yy, xx = np.unravel_index(np.argmax(delta), delta.shape)
            x0, y0 = meta["target_xy"]
            sep = math.hypot(xx - x0, yy - y0)
            metrics.update(
                spearman=scores,
                third_size=n,
                peak_pixel=[int(xx), int(yy)],
                peak_distance_pixels=round(sep, 4),
                passes=bool(
                    len(matches) >= 6 and scores["center"] >= 0.65 and sep <= 1.5
                ),
                matches=[
                    {"date": str(d), "ztf_mag": round(y, 5), "tess_bjd": float(t[k])}
                    for d, k, y in matches
                ],
            )
        out[b] = metrics
    return {
        "sector": meta["sector"],
        "target_xy": meta["target_xy"],
        "six_hour_bins": len(t),
        "bands": out,
    }


def run(fetch, local):
    SOURCE.mkdir(exist_ok=True)
    _, pre = curves(NAME)
    return [
        compare(*reduction(s, c, d, fetch, local), pre)
        for s, c, d in sectors(fetch, local)
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--fetch", action="store_true")
    group.add_argument("--reduce-from", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "tess-results.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch, args.reduce_from))
