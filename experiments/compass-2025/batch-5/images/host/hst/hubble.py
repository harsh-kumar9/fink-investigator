"""Extract the two pre-event HAP science planes and compare to propagated Gaia stars."""

import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import field

arc = field.archive
NAME, RA, DEC = field.NAME, field.RA, field.DEC
MAST = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/hst_16691_42_wfc3_uvis_"
FILTERS = ("f555w", "f814w")
arc.SOURCES = HERE / "sources"


def cutout(band, center, fetch):
    url = MAST + band + "_ierg42_drc.fits"
    prefix = arc.SOURCES / band
    saved = prefix.with_suffix(".fits.gz")
    meta = prefix.with_suffix(".json")
    # The science array occupies the first extension in both public visit products.
    hdrs = io.BytesIO(arc.blob(band + "-headers.bin", url, fetch, (0, 86399)))
    primary = fits.Header.fromfile(hdrs)
    head = fits.Header.fromfile(hdrs)
    if (
        head["EXTNAME"] != "SCI"
        or head["BITPIX"] != -32
        or primary["FILTER"].lower() != band
    ):
        raise ValueError("Unexpected HAP mosaic")
    wcs = WCS(head)
    x, y = map(float, wcs.world_to_pixel_values(*center))
    side = 321
    radius = side // 2
    x0, y0 = round(x) - radius, round(y) - radius
    width, height = head["NAXIS1"], head["NAXIS2"]
    if min(x0, y0) < 0 or x0 + side > width or y0 + side > height:
        raise ValueError("Cutout outside mosaic")
    start = hdrs.tell() + 4 * width * y0
    end = start + 4 * width * side - 1
    if fetch and not meta.exists():
        raw = arc.pilot.curl(url, ["--range", f"{start}-{end}"])
        if len(raw) != end - start + 1:
            raise ValueError("Range ignored")
        crop = (
            np.frombuffer(raw, dtype=">f4")
            .reshape(side, width)[:, x0 : x0 + side]
            .copy()
        )
        out = wcs.slice((slice(y0, y0 + side), slice(x0, x0 + side))).to_header()
        for k in [
            "BUNIT",
            "PHOTFLAM",
            "PHOTPLAM",
            "WCSNAME",
            "RMS_RA",
            "RMS_DEC",
            "NMATCHES",
        ]:
            out[k] = head[k]
        out["FILTER"] = primary["FILTER"]
        out["DATE-OBS"] = primary["DATE-OBS"]
        out["EXPTIME"] = primary["EXPTIME"]
        buf = io.BytesIO()
        fits.writeto(buf, crop, out)
        result = buf.getvalue()
        saved.write_bytes(gzip.compress(result, mtime=0))
        arc.pilot.json_write(
            meta,
            {
                "url": url,
                "science_byte_range": [start, end],
                "range_sha256": hashlib.sha256(raw).hexdigest(),
                "cutout_x0_y0": [x0, y0],
                "cutout_sha256": hashlib.sha256(result).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    raw = gzip.decompress(saved.read_bytes())
    m = json.loads(meta.read_text())
    if (
        m["url"] != url
        or m["science_byte_range"] != [start, end]
        or m["cutout_x0_y0"] != [x0, y0]
        or m["cutout_sha256"] != hashlib.sha256(raw).hexdigest()
    ):
        raise ValueError("HST cutout mismatch")
    with fits.open(io.BytesIO(raw)) as h:
        image = h[0].data.copy()
        cutw = WCS(h[0].header)
        if (
            image.shape != (side, side)
            or np.max(
                np.abs(np.array(cutw.world_to_pixel_values(*center)) - [x - x0, y - y0])
            )
            > 1e-5
        ):
            raise ValueError("HST trimmed WCS changed")
    return image, cutw, head, primary


def gaia(fetch):
    cols = "source_id,ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag,ruwe"
    sql = (
        f"SELECT TOP 1000 {cols} FROM gaia_dr3_source WHERE "
        "CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',39.49389185,-1.8497931,0.05))=1"
    )
    url = (
        arc.pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "CSV", "QUERY": sql}
        )
    )
    rows = arc.compass_direct.capture(
        arc.SOURCES, "gaia", url, set(cols.split(",")), fetch
    )
    if len(rows) >= 1000:
        raise ValueError("Truncated Gaia cone")
    return rows


def hap(key, fetch):
    u = MAST + key + "_ierg42_point-cat.ecsv"
    return Table.read(
        io.BytesIO(arc.blob(key + "-point.ecsv", u, fetch)), format="ascii.ecsv"
    )


def arm_aperture(image, wcs, center, cat):
    pix = abs(wcs.proj_plane_pixel_scales()[0].value) * 3600
    cx, cy = map(float, wcs.world_to_pixel_values(*center))

    def aperture(x, y):
        rad = math.ceil(0.46 / pix)
        gx, gy = np.meshgrid(
            np.arange(-rad, rad + 1) + round(x), np.arange(-rad, rad + 1) + round(y)
        )
        d = np.hypot(gx - x, gy - y) * pix
        values = image[gy, gx]
        sky = np.median(values[(d >= 0.25) & (d <= 0.45)])
        return float(np.sum(values[d <= 0.15] - sky)), float(sky)

    flux, sky = aperture(cx, cy)
    points = np.column_stack(wcs.world_to_pixel_values(cat["RA"], cat["DEC"]))
    grid = {}
    for px in np.arange(-4.5, 4.51, 0.6):
        for py in np.arange(-4.5, 4.51, 0.6):
            if not 0.7 <= math.hypot(px, py) <= 4.5:
                continue
            x, y = cx + px / pix, cy + py / pix
            if np.hypot(points[:, 0] - x, points[:, 1] - y).min() * pix < 0.35:
                continue
            f, bg = aperture(x, y)
            if 0.5 * sky <= bg <= 1.5 * sky:
                grid[(float(px), float(py))] = f
    if len(grid) < 15:
        raise ValueError("Too few background controls")
    values = np.array(list(grid.values()))
    med = np.median(values)
    mad = 1.4826 * np.median(abs(values - med))
    return {
        "flux_electrons_per_s": round(flux, 5),
        "controls": len(grid),
        "controls_as_bright": int(np.sum(values >= flux)),
        "background_median_flux": round(float(med), 5),
        "background_robust_sigma": round(float(mad), 5),
        "excess_over_local_median_sigma": round(float((flux - med) / mad), 2),
    }, {pos: f >= flux for pos, f in grid.items()}


def run(fetch=False, plot=False):
    arc.SOURCES.mkdir(exist_ok=True)
    ra, dec, position = arc.centers()
    esa = "https://hst.esac.esa.int/tap-server/tap/sync?"
    fields = "obs_collection,obs_id,instrument_name,t_min,t_max"
    sql = (
        "SELECT TOP 100 " + fields + " FROM ivoa.obscore WHERE "
        "1=CONTAINS(POINT('ICRS',39.49389185,-1.8497931),s_region)"
    )
    rows = arc.compass_direct.capture(
        arc.SOURCES,
        "coverage",
        esa
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": sql}
        ),
        set(fields.split(",")),
        fetch,
    )
    if len(rows) >= 100 or not all(
        any(r["obs_id"] == f"hst_16691_42_wfc3_uvis_{b}_ierg42" for r in rows)
        for b in FILTERS
    ):
        raise ValueError("Missing HST visit")
    total = hap("total", fetch)
    rawstars = gaia(fetch)
    refyear = 2000 + (59821 - 51544.5) / 365.25
    anchors = []
    for s in rawstars:
        if not all(s[k] for k in ("pmra", "pmdec", "ruwe")) or float(s["ruwe"]) > 1.4:
            continue
        dt = refyear - float(s["ref_epoch"])
        a, d = field.offset(
            float(s["ra"]),
            float(s["dec"]),
            float(s["pmra"]) * dt / 1000,
            float(s["pmdec"]) * dt / 1000,
        )
        dx = (np.array(total["RA"]) - a) * 3600 * math.cos(math.radians(d))
        dy = (np.array(total["DEC"]) - d) * 3600
        sep = np.hypot(dx, dy)
        i = np.argmin(sep)
        if sep[i] > 0.8 or any(
            int(total["Flags_" + b][i]) not in (0, 1) for b in FILTERS
        ):
            continue
        anchors.append(
            (s["source_id"], int(total["ID"][i]), float(dx[i]), float(dy[i]))
        )
    if len(anchors) < 3 or len({a[1] for a in anchors}) != len(anchors):
        raise ValueError("Ambiguous Gaia ties")
    residual = np.array([a[2:] for a in anchors])
    shift = np.median(residual, axis=0)
    center = field.offset(
        ra, dec, *shift
    )  # Predicted coordinates on the native HAP frame
    results = {
        "gaia_ties": [list(row) for row in anchors],
        "hst_minus_gaia_east_north_arcsec": shift.tolist(),
        "tie_max_residual_after_translation_arcsec": float(
            np.hypot(*(residual - shift).T).max()
        ),
        "bands": {},
    }
    plots = []
    grids = {}
    for band in FILTERS:
        image, wcs, header, pri = cutout(band, (ra, dec), fetch)
        catalog = hap(band, fetch)
        valid = catalog[
            np.isin(catalog["Flags"], [0, 1])
            & (catalog["MagAp2"] > 0)
            & (catalog["MagErrAp2"] > 0)
        ]
        x = (
            (np.array(valid["RA"]) - center[0])
            * 3600
            * math.cos(math.radians(center[1]))
        )
        y = (np.array(valid["DEC"]) - center[1]) * 3600
        sep = np.hypot(x, y)
        close = []
        for i in np.argsort(sep)[:8]:
            close.append(
                dict(
                    id=int(valid["ID"][i]),
                    sep_arcsec=round(float(sep[i]), 3),
                    AB_0p15arcsec=round(float(valid["MagAp2"][i]), 3),
                    err=round(float(valid["MagErrAp2"][i]), 3),
                    flags=int(valid["Flags"][i]),
                )
            )
        ap, grids[band] = arm_aperture(image, wcs, center, valid)
        zpt = (
            -2.5 * math.log10(header["PHOTFLAM"])
            - 5 * math.log10(header["PHOTPLAM"])
            - 2.408
        )
        ap["empirical_5sigma_AB_0p15arcsec"] = round(
            zpt - 2.5 * math.log10(5 * ap["background_robust_sigma"]), 2
        )
        results["bands"][band] = {
            "date": pri["DATE-OBS"],
            "exposure_seconds": pri["EXPTIME"],
            "wcs": header["WCSNAME"],
            "original_rms_mas": [header["RMS_RA"], header["RMS_DEC"]],
            "nearest_accepted": close,
            "aperture": ap,
        }
        plots.append((band, image, wcs, center, valid[np.argsort(sep)[:8]]))
    common = grids[FILTERS[0]].keys() & grids[FILTERS[1]].keys()
    results["shared_control_positions"] = len(common)
    results["controls_as_bright_in_both"] = sum(
        all(grids[b][p] for b in FILTERS) for p in common
    )
    if plot:
        import matplotlib.pyplot as plt

        f, axes = plt.subplots(1, 2, figsize=(6.7, 3.6), layout="constrained")
        for ax, (band, im, wcs, center, nearby) in zip(axes, plots):
            x, y = map(float, wcs.world_to_pixel_values(*center))
            pix = wcs.proj_plane_pixel_scales()[0].value * 3600
            cut = im[round(y) - 32 : round(y) + 33, round(x) - 32 : round(x) + 33]
            vmin, vmax = np.percentile(cut, [7, 99])
            ax.imshow(im, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
            ax.add_patch(
                plt.Circle(
                    (x, y),
                    position["p90_nightly_radius_arcsec"] / pix,
                    fill=False,
                    color="#ff6e4e",
                    lw=1.2,
                )
            )
            px, py = wcs.world_to_pixel_values(nearby["RA"], nearby["DEC"])
            ax.scatter(
                px,
                py,
                marker="o",
                s=22,
                facecolors="none",
                edgecolors="#65d7fc",
                linewidth=0.7,
            )
            ax.set(
                xlim=(x - 1.35 / pix, x + 1.35 / pix),
                ylim=(y - 1.35 / pix, y + 1.35 / pix),
                xticks=[],
                yticks=[],
                title=band.upper(),
            )
            x0, y0 = x - 1.18 / pix, y - 1.09 / pix
            ax.plot([x0, x0 + 0.5 / pix], [y0, y0], color="white", lw=2)
            ax.text(x0, y0 + 0.05 / pix, "0.5″", color="white", fontsize=8)
        f.suptitle("Hubble’s 2022 view of the transient site", fontsize=13)
        f.savefig(HERE / "field.png", dpi=125)
        plt.close(f)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--plot", action="store_true")
    a = p.parse_args()
    obj = run(a.fetch, a.plot)
    path = HERE / "results.json"
    if a.fetch:
        arc.pilot.json_write(path, obj)
    else:
        arc.compass.check_replay(obj, json.loads(path.read_text()), str(path))
    print(len(obj["gaia_ties"]), "Gaia ties, two Hubble images verified")
