"""Compare a local original-image crop to the archived IBE crop, offline."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import image_fetch
import onset

NAMES = {"2018-07-16": "historic", "2025-04-24": "april"}
POS = onset.POS


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fixture(path, label, data=None):
    name = path / (label + ".gz")
    if data is not None:
        name.write_bytes(gzip.compress(data, compresslevel=6, mtime=0))
    return gzip.decompress(name.read_bytes())


def run(originals=None, service=None):
    stored = {} if originals else json.loads((HERE / "results.json").read_text())
    sample, _ = onset.choose(onset.tables(False)[0])
    result = {}
    for day, name in NAMES.items():
        (row,) = [r for _, d, r in sample if str(d) == day]
        if originals:
            full = originals[name].read_bytes()
            original_sha = sha(full)
            total = len(full)

            def whole(url, start, end, budget):
                if len(full) > budget:
                    raise ValueError("Full image exceeds byte ceiling")
                return (
                    full,
                    len(full),
                    {
                        "requested_bytes": [start, end],
                        "http_status": 200,
                        "received_bytes": len(full),
                        "sha256": sha(full),
                    },
                )

            whole_fits, whole_receipt = image_fetch.from_original(
                row, POS, read_range=whole
            )
        else:
            original_sha = stored[name]["original_sha256"]
            total = stored[name]["original_bytes"]
        pieces = []

        def sections(url, start, end, budget):
            if (
                not originals
                and [start, end]
                != stored[name]["ranges"][len(pieces)]["requested_bytes"]
            ):
                raise ValueError("Requested source section changed")
            if end - start + 1 > budget:
                raise ValueError("Source section exceeds byte ceiling")
            label = name + ("-header" if start == 0 else "-strip")
            raw = fixture(HERE, label, full[start : end + 1] if originals else None)
            if len(raw) != end - start + 1:
                raise ValueError("Incomplete source section")
            receipt = {
                "requested_bytes": [start, end],
                "http_status": 206,
                "received_bytes": len(raw),
                "sha256": sha(raw),
            }
            pieces.append(receipt)
            return raw, total, receipt

        cut, rec = image_fetch.from_original(row, POS, read_range=sections)
        if originals and cut != whole_fits:
            raise ValueError("Byte-range and full-image crops disagree")
        with fits.open(io.BytesIO(cut)) as f:
            img = f[0].data
            star_xy = WCS(f[0].header).all_world2pix([POS], 0)[0]
            record = {
                "original_url": rec["original_url"],
                "expid": int(row["expid"]),
                "pid": row["pid"],
                "original_bytes": total,
                "original_sha256": original_sha,
                "ranges": pieces,
                "slice_xy_zero_based_half_open": rec["slice_xy_zero_based_half_open"],
                "shape": list(img.shape),
                "target_xy": star_xy.tolist(),
                "crop_sha256": sha(cut),
            }
            if name == "april":
                official = fixture(
                    HERE, "april-service", service.read_bytes() if originals else None
                )
                with fits.open(io.BytesIO(official)) as sf:
                    if sf[0].data.shape != img.shape or not np.array_equal(
                        sf[0].data, img, equal_nan=True
                    ):
                        raise ValueError("Local and service pixels differ")
                    yy, xx = np.mgrid[: img.shape[0], : img.shape[1]]
                    grid = np.column_stack([xx.ravel(), yy.ravel()])
                    original_sky = WCS(f[0].header).all_pix2world(grid, 0)
                    service_sky = WCS(sf[0].header).all_pix2world(grid, 0)
                    max_delta = float(np.max(np.abs(original_sky - service_sky)))
                    if max_delta > 1.0e-9:
                        raise ValueError("Local and service WCS differ")
                    record["service"] = {
                        "url": onset.aperture.cutout_url(row, *POS),
                        "sha256": sha(official),
                        "identical_pixels": img.size,
                        "max_abs_sky_axis_delta_deg": max_delta,
                    }
        result[name] = record
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--capture",
        nargs=3,
        metavar=("JULY_2018_FULL", "APRIL_2025_FULL", "APRIL_SERVICE_CUT"),
        type=Path,
    )
    args = p.parse_args()
    if args.capture:
        historic, april, service = args.capture
        result = run({"historic": historic, "april": april}, service)
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    elif run() != json.loads((HERE / "results.json").read_text()):
        raise ValueError("Fallback proof changed")
    print(
        "Verified both original-image crops and all 26,896 April service pixels and sky coordinates."
    )
