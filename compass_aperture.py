"""Measure a known position on an original ZTF image with same-detector standards."""

import io
import math
from urllib.parse import urlencode

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import pilot


def cutout_url(row, ra, dec, size=165):
    stamp = row["filefracday"]
    file = f"ztf_{stamp}_{int(row['field']):06}_{row['filtercode']}_c{int(row['ccdid']):02}_{row['imgtypecode']}_q{row['qid']}_sciimg.fits"
    path = "/".join([stamp[:4], stamp[4:8], stamp[8:], file])
    params = urlencode(
        {"center": f"{ra},{dec}", "size": f"{size}arcsec", "gzip": "false"}
    )
    return f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{path}?{params}"


def robust_scatter(values):
    values = np.array(values)
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def aperture(img, x, y):
    # Subpixel weights keep the fixed sky position independent of the brightest pixel.
    cx, cy = round(x), round(y)
    if not (12 <= cx < img.shape[1] - 12 and 12 <= cy < img.shape[0] - 12):
        return None
    stamp = img[cy - 11 : cy + 12, cx - 11 : cx + 12]
    yy, xx = np.mgrid[-11:12, -11:12]
    dx, dy = xx + cx - x, yy + cy - y
    r2 = dx**2 + dy**2
    sky = stamp[(r2 >= 49) & (r2 <= 121)]
    if not (np.isfinite(stamp[r2 < 21]).all() and np.isfinite(sky).all()):
        return None
    weights = np.zeros(stamp.shape)
    for i in (-0.375, -0.125, 0.125, 0.375):
        for j in (-0.375, -0.125, 0.125, 0.375):
            weights += ((dx + i) ** 2 + (dy + j) ** 2 <= 16) / 16
    return float(np.sum((stamp[r2 < 21] - np.median(sky)) * weights[r2 < 21]))


def star_list(rows, meta, ra, dec, radius=65):
    field = [
        r
        for r in rows
        if all(r[k] == meta[k] for k in ("field", "ccdid", "qid"))
        and "z" + {"1": "g", "2": "r", "3": "i"}[r["fid"]] == meta["filtercode"]
    ]
    sky = np.array(
        [(float(r["ra"]), float(r["dec"])) for r in field], dtype=float
    ).reshape(-1, 2)
    mags = np.array([float(r["meanmag"]) for r in field])
    standards = []
    for i, r in enumerate(field):
        mag = mags[i]
        if not (
            15 <= mag <= 18.5
            and int(r["nobs"]) >= 30
            and 14 <= math.hypot(*pilot.offsets(ra, dec, *sky[i])) <= radius
        ):
            continue
        if any(
            j != i
            and mags[j] <= mag + 2
            and math.hypot(*pilot.offsets(*sky[i], *sky[j])) < 8
            for j in range(len(field))
        ):
            continue
        standards.append(i)
    return field, sky, standards


def measure(raw, meta, position, star_rows, standard_radius=65):
    with fits.open(io.BytesIO(raw)) as hdus:
        h = hdus[0].header
        image = np.asarray(hdus[0].data, dtype=float)
        if (
            int(h["EXPID"]) != int(meta["expid"])
            or int(h["QID"]) != int(meta["qid"])
            or int(h["CCD_ID"]) != int(meta["ccdid"])
            or int(h["FIELDID"]) != int(meta["field"])
            or h["FILTER"] != "ZTF_" + meta["filtercode"][1]
            or abs(h["OBSMJD"] + 2400000.5 - float(meta["obsjd"])) > 1 / 86400
        ):
            raise ValueError("Image metadata does not match the requested exposure")
        wcs = WCS(h)
        target = wcs.all_world2pix([position], 0)[0]
        # At a quadrant edge IRSA can return an asymmetric cutout.
        field, sky, keep = star_list(star_rows, meta, *position, radius=standard_radius)
        xy = wcs.all_world2pix(sky, 0) if len(sky) else np.empty((0, 2))
        calibrators = []
        for i in keep:
            c = aperture(image, *xy[i])
            if c is not None and c > 0:
                calibrators.append((field[i]["oid"], float(field[i]["meanmag"]), c))
        result = {"standards": len(calibrators), "target_xy": list(target)}
        counts = aperture(image, *target)
        if counts is None or len(calibrators) < 3:
            return {**result, "status": "insufficient usable target or standards"}
        blanks = []
        for y in range(14, image.shape[0] - 13, 10):
            for x in range(14, image.shape[1] - 13, 10):
                if np.hypot(x - target[0], y - target[1]) < 13 or (
                    len(xy) and np.min(np.hypot(xy[:, 0] - x, xy[:, 1] - y)) < 9
                ):
                    continue
                c = aperture(image, x, y)
                if c is not None:
                    blanks.append(c)
        if len(blanks) < 12:
            return {
                **result,
                "status": "insufficient offset apertures",
                "blanks": len(blanks),
            }
        blank = float(np.median(blanks))
        calibrators = [
            (oid, mag, c - blank) for oid, mag, c in calibrators if c > blank
        ]
        result["standards"] = len(calibrators)
        if len(calibrators) < 3:
            return {**result, "status": "insufficient positive standards"}
        anchors = [(oid, mag + 2.5 * np.log10(c)) for oid, mag, c in calibrators]
        zeropoint = np.median([v for _, v in anchors])
        scatter = robust_scatter([v for _, v in anchors])
        signal = counts - blank
        sky_sigma = robust_scatter(blanks)
        sigma = np.sqrt(sky_sigma**2 + max(0.0, signal) / float(h["GAIN"]))
        if not (sigma > 0 and sky_sigma > 0 and float(h["GAIN"]) > 0):
            raise ValueError("Invalid aperture noise")
        sn = signal / sigma
        control = None
        if len(anchors) >= 4:
            oid, mag, c = min(calibrators, key=lambda v: (abs(v[1] - 18), v[0]))
            independent_zp = np.median([v for key, v in anchors if key != oid])
            control = {
                "oid": oid,
                "catalog_mag": mag,
                "measured_mag": independent_zp - 2.5 * np.log10(c),
                "snr": c / np.sqrt(sky_sigma**2 + c / float(h["GAIN"])),
            }
        return {
            **result,
            "status": "measured",
            "blanks": len(blanks),
            "counts": counts,
            "blank_median": float(np.median(blanks)),
            "blank_sigma": sky_sigma,
            "signal": signal,
            "sigma": sigma,
            "snr": sn,
            "zero_point": zeropoint,
            "zero_point_scatter": scatter,
            "standard_residuals": {i: round(v - zeropoint, 4) for i, v in anchors},
            "held_out_star": control,
            "mag": float(zeropoint - 2.5 * np.log10(signal)) if sn >= 5 else None,
            "magerr": (
                float(np.hypot(2.5 / np.log(10) / sn, scatter)) if sn >= 5 else None
            ),
            "five_sigma_limit": float(zeropoint - 2.5 * np.log10(5 * sky_sigma)),
        }
