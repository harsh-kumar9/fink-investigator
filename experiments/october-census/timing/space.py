"""Reduce a MAST TESSCut download and replay its fixed-frequency aperture comparison."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.stats import spearmanr

from epochs import Ephemeris, HERE, POS

CUTOUT = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=280.305911&dec=9.7241705&x=9&y=9&units=px&sector=80"
APERTURES = {
    "center": (3, 3),
    "lower_left": (0, 0),
    "lower_right": (0, 6),
    "upper_left": (6, 0),
    "upper_right": (6, 6),
}
DATA = HERE / "tess-apertures.csv.gz"
RECEIPT = HERE / "tess-source.json"


def extract(path):
    with ZipFile(path) as archive:
        (name,) = archive.namelist()
        raw = archive.read(name)
    stream = io.StringIO(newline="")
    with fits.open(io.BytesIO(raw)) as hdus:
        data = hdus[1].data
        hdr = hdus[1].header
        if int(hdus[0].header["SECTOR"]) != 80 or hdr["TIMESYS"] != "TDB":
            raise ValueError("Unexpected TESS sector or time system")
        x, y = WCS(hdus[2].header).world_to_pixel(POS)
        if max(abs(x - 4), abs(y - 4)) > 0.5:
            raise ValueError("Target outside the central pixel")
        cube = np.asarray(data["FLUX"], dtype=float)
        quality = data["QUALITY"] == 0
        good = (
            quality & np.all(np.isfinite(cube), axis=(1, 2)) & np.isfinite(data["TIME"])
        )
        pix = cube[good]
        border = np.ones((9, 9), dtype=bool)
        border[1:-1, 1:-1] = False
        sky = np.median(pix[:, border], axis=1)
        times = data["TIME"][good] + hdr["BJDREFI"] + hdr["BJDREFF"]
        counts = [
            pix[:, r : r + 3, c : c + 3].sum(axis=(1, 2)) - 9 * sky
            for r, c in APERTURES.values()
        ]
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["bjd", *APERTURES])
        writer.writerows(zip(times, *counts))
        meta = {
            "url": CUTOUT,
            "retrieved_utc": datetime.fromtimestamp(
                path.stat().st_mtime, timezone.utc
            ).isoformat(timespec="seconds"),
            "fits_file": name,
            "fits_sha256": hashlib.sha256(raw).hexdigest(),
            "sector": 80,
            "ccd": int(hdus[0].header["CCD"]),
            "camera": int(hdus[0].header["CAMERA"]),
            "target_xy_zero_based": [float(x), float(y)],
            "input_cadences": len(data),
            "quality_zero": int(quality.sum()),
            "retained": int(good.sum()),
        }
    encoded = stream.getvalue().encode()
    DATA.write_bytes(gzip.compress(encoded, mtime=0))
    meta["aperture_csv_sha256"] = hashlib.sha256(encoded).hexdigest()
    RECEIPT.write_text(json.dumps(meta, indent=2) + "\n")


def six_hour_bins():
    raw = gzip.decompress(DATA.read_bytes())
    if (
        hashlib.sha256(raw).hexdigest()
        != json.loads(RECEIPT.read_text())["aperture_csv_sha256"]
    ):
        raise ValueError("Changed TESS aperture file")
    data = np.loadtxt(io.BytesIO(raw), delimiter=",", skiprows=1)
    bins = np.floor(data[:, 0] * 4)
    return np.array(
        [
            np.median(data[bins == key], axis=0)
            for key in np.unique(bins)
            if np.count_nonzero(bins == key) >= 5
        ]
    )


def compare():
    binned = six_hour_bins()
    fit = Ephemeris()
    t = binned[:, 0] - binned[0, 0]
    trend = np.column_stack([np.ones(len(t)), t])

    def residuals(design, y):
        return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]

    metrics = {}
    for i, aperture in enumerate(APERTURES, start=1):
        y = binned[:, i]
        base = residuals(trend, y)
        row = {}
        for name, period in [
            ("ztf", 1 / fit.frequency),
            ("atlas_alternative", 2.408143),
        ]:
            angle = 2 * np.pi * t / period
            model = np.column_stack([trend, np.sin(angle), np.cos(angle)])
            resid = residuals(model, y)
            row[name + "_fractional_error_reduction"] = float(
                1 - (resid @ resid) / (base @ base)
            )
        row["fixed_waveform_rho"] = float(
            spearmanr(base, -fit.mag(binned[:, 0])).statistic
        )
        metrics[aperture] = row
    return {
        "bins": len(t),
        "first_bjd": binned[0, 0],
        "last_bjd": binned[-1, 0],
        "apertures": metrics,
    }


def windows():
    data = six_hour_bins()
    t = data[:, 0]
    result = json.loads((HERE / "dips.json").read_text())
    fit = result["fits"][result["training_winner"]]
    halfwidth = fit["duration_days"] / 2
    spacing = fit["period_days"] / 2
    events = []
    for j in range(
        int(np.ceil((t[0] - fit["center_hjd"]) / spacing)),
        int(np.floor((t[-1] - fit["center_hjd"]) / spacing)) + 1,
    ):
        center = fit["center_hjd"] + j * spacing
        dt = (t - center) / halfwidth
        left, mid, right = (
            (dt >= -3.25) & (dt <= -1.25),
            np.abs(dt) <= 1,
            (dt >= 1.25) & (dt <= 3.25),
        )
        counts = [int(x.sum()) for x in (left, mid, right)]
        event = {
            "predicted_bjd_approx": center,
            "window": "deep" if j % 2 == 0 else "half-cycle",
            "counts_left_in_right": counts,
        }
        if min(counts) >= 2:
            deficits = {}
            for i, name in enumerate(APERTURES, start=1):
                y = data[:, i]
                baseline = np.interp(
                    t[mid],
                    [np.median(t[left]), np.median(t[right])],
                    [np.median(y[left]), np.median(y[right])],
                )
                deficits[name] = float(np.median(baseline - y[mid]))
            event["deficit_electrons_per_sec"] = deficits
        events.append(event)
    return events


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extract", type=Path, help="Reduce a sector-80 MAST TESSCut ZIP"
    )
    parser.add_argument(
        "--windows", action="store_true", help="Test the frozen deep-fade windows"
    )
    args = parser.parse_args()
    if args.extract:
        extract(args.extract)
    print(
        json.dumps(windows() if args.windows else compare(), indent=2, allow_nan=False)
    )
