"""Archive quality-zero 9x9 TESS pixels and compare them with the ZTF observations."""

import argparse
from datetime import datetime, timezone
import hashlib, io, json
from pathlib import Path
from urllib.parse import urlencode
from zipfile import ZipFile

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import requests

from analyze_vars import HERE, NAMES, load, pilot
import tess as previous_tess

BASE = "https://mast.stsci.edu/tesscut/api/v0.1/"
SRC = HERE / "sources"
SELECTION = {"ZTF17aaaoxxi": [6, 33, 87], "ZTF18abcckqd": [14, 15, 16]}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def sectors(case, fetch):
    prefix = SRC / (case["id"] + "-tess-sectors")
    raw = prefix.with_suffix(".json")
    receipt = prefix.with_suffix(".receipt.json")
    url = (
        BASE
        + "sector?"
        + urlencode({"ra": case["ra"], "dec": case["dec"], "radius": "0s"})
    )
    if fetch and not receipt.exists():
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        response.json()["results"]
        raw.write_bytes(response.content)
        pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256": sha(response.content),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    meta = json.loads(receipt.read_text())
    body = raw.read_bytes()
    if meta["url"] != url or sha(body) != meta["sha256"]:
        raise ValueError("TESS sector inventory changed")
    return {
        int(v["sector"]): (int(v["camera"]), int(v["ccd"]))
        for v in json.loads(body)["results"]
    }


def cube(case, sector, cam, ccd, fetch):
    stem = case["id"] + f"-tess-{sector}"
    receipt = SRC / (stem + ".json")
    compressed = SRC / (stem + ".npz")
    url = (
        BASE
        + "astrocut?"
        + urlencode(
            {
                "ra": case["ra"],
                "dec": case["dec"],
                "x": 9,
                "y": 9,
                "units": "px",
                "sector": sector,
            }
        )
    )
    if fetch and not receipt.exists():
        response = requests.get(url, timeout=180)
        response.raise_for_status()
        with ZipFile(io.BytesIO(response.content)) as z:
            (member,) = z.namelist()
            raw = z.read(member)
        with fits.open(io.BytesIO(raw)) as h:
            header = h[0].header
            table = h[1].data
            if [int(header[k]) for k in ("SECTOR", "CAMERA", "CCD")] != [
                sector,
                cam,
                ccd,
            ] or h[1].header["TIMESYS"] != "TDB":
                raise ValueError("Unexpected sector, detector or time scale")
            x, y = WCS(h[2].header).world_to_pixel(
                SkyCoord(case["ra"], case["dec"], unit="deg")
            )
            if max(abs(x - 4), abs(y - 4)) > 1.5:
                raise ValueError("Unexpected TESS image center")
            flux = np.asarray(table["FLUX"], dtype="<f4")
            t = (
                np.array(table["TIME"], dtype="<f8")
                + h[1].header["BJDREFI"]
                + h[1].header["BJDREFF"]
            )
            mask = (
                (table["QUALITY"] == 0)
                & np.isfinite(t)
                & np.all(np.isfinite(flux), axis=(1, 2))
            )
            if flux.shape[1:] != (9, 9) or not np.any(mask):
                raise ValueError("Missing TESS pixels")
            np.savez_compressed(compressed, bjd=t[mask], flux=flux[mask])
            meta = {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sector": sector,
                "camera": cam,
                "ccd": ccd,
                "fits_name": member,
                "fits_sha256": sha(raw),
                "cube_sha256": sha(compressed.read_bytes()),
                "target_xy": [float(x), float(y)],
                "all_cadences": len(t),
                "retained": int(sum(mask)),
                "first_bjd": float(t[mask][0]),
                "last_bjd": float(t[mask][-1]),
            }
            pilot.json_write(receipt, meta)
    m = json.loads(receipt.read_text())
    if (
        m["url"] != url
        or m["sector"] != sector
        or m["camera"] != cam
        or m["ccd"] != ccd
        or sha(compressed.read_bytes()) != m["cube_sha256"]
    ):
        raise ValueError("TESS archive differs")
    with np.load(compressed, allow_pickle=False) as z:
        t, flux = z["bjd"], z["flux"].astype(float)
    border = np.ones((9, 9), dtype=bool)
    border[1:-1, 1:-1] = False
    flux -= np.median(flux[:, border], axis=1)[:, None, None]
    return m, t, flux


def bin6(t, pix):
    idx = np.floor(t * 4).astype(int)
    return np.array(
        [
            [np.median(t[idx == i]), *np.median(pix[idx == i], axis=0).ravel()]
            for i in np.unique(idx)
            if np.sum(idx == i) >= 5
        ]
    )


def run(name, fetch):
    case, pre, _ = load(name)
    options = sectors(case, fetch)
    out = {"available_sectors": sorted(options), "evaluated": []}
    for s in SELECTION[name]:
        m, t, px = cube(case, s, *options[s], fetch)
        six = bin6(t, px)
        metrics = previous_tess.compare(m, six, pre)
        metrics.update(
            {
                "first_utc_date": str(pilot.day_of(t[0])),
                "last_utc_date": str(pilot.day_of(t[-1])),
                "cadences": len(t),
                "median_cadence_minutes": round(float(np.median(np.diff(t)) * 1440), 2),
            }
        )
        out["evaluated"].append(metrics)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "tess.json")
    args = p.parse_args()
    pilot.json_write(args.out, {n: run(n, args.fetch) for n in NAMES})
