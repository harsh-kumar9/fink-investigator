"""Localize ZTF25absffwy in the pre-event galaxy with two instruments and spectroscopy."""

import argparse
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import texture

HERE = Path(__file__).resolve().parent
YEAR = HERE.parents[2]
sys.path.insert(0, str(YEAR / "batch-4/images/host"))
import check as archive

NAME = "ZTF25absffwy"
RA, DEC = 39.4938772, -1.8497861
VIEW = "https://www.legacysurvey.org/viewer/"
SKY_BOX = {"ralo": 39.479, "rahi": 39.509, "declo": -1.865, "dechi": -1.834}
DATA = HERE / "sources"
archive.SOURCES = DATA
archive.NAME = NAME


def tap(key, table, fields, condition, fetch):
    sql = f"SELECT TOP 500 {fields} FROM {table} WHERE {condition}"
    url = "https://datalab.noirlab.edu/tap/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": sql}
    )
    data = archive.compass_direct.capture(DATA, key, url, set(fields.split(",")), fetch)
    if len(data) >= 500:
        raise ValueError("Query hit row cap")
    return data


def image(name, size, scale, fetch):
    params = {
        "ra": RA,
        "dec": DEC,
        "layer": name,
        "pixscale": scale,
        "size": size,
        "bands": "grz",
    }
    raw = archive.blob(name + ".fits", VIEW + "fits-cutout?" + urlencode(params), fetch)
    with fits.open(io.BytesIO(raw)) as hdus:
        data, hdr = hdus[0].data.copy(), hdus[0].header
        if data.shape != (3, size, size) or hdr["BANDS"] != "grz":
            raise ValueError("Wrong band coverage or image size")
        wcs = WCS(hdr, naxis=2)
        center = wcs.world_to_pixel_values(RA, DEC)
        if max(abs(float(v) - (size - 1) / 2) for v in center) > 2:
            raise ValueError("Incorrect WCS")
    return data, wcs, scale


def sample(images, location):
    data, wcs, scale = images
    x, y = map(float, wcs.world_to_pixel_values(*location))
    ix, iy = round(x), round(y)
    k = math.ceil(2 / scale)
    stamp = data[:, iy - k : iy + k + 1, ix - k : ix + k + 1]
    if (
        stamp.shape != (3, 2 * k + 1, 2 * k + 1)
        or not np.isfinite(stamp).all()
        or np.any(np.std(stamp, axis=(1, 2)) <= 0)
    ):
        return None
    gy, gx = np.mgrid[-k : k + 1, -k : k + 1]
    dist = scale * np.hypot(gx + ix - x, gy + iy - y)
    return (
        np.median(stamp[:, dist <= 1], axis=1),
        np.median(stamp[:, dist <= 0.5], axis=1)
        - np.median(stamp[:, (dist >= 1) & (dist <= 2)], axis=1),
    )


def offset(ra, dec, x, y):
    return ra + x / (3600 * math.cos(math.radians(dec))), dec + y / 3600


def image_test(images, center, host, catalog):
    source_level, contrast = sample(images, center)
    outer = []
    for dx in range(-45, 46, 6):
        for dy in range(-45, 46, 6):
            pos = offset(*center, dx, dy)
            x, _ = archive.pilot.offsets(*pos, float(host["ra"]), float(host["dec"]))
            if (
                abs(x) < 30
                or min(archive.location(*pos, r["ra"], r["dec"]) for r in catalog) <= 5
            ):
                continue
            values = sample(images, pos)
            if values is not None:
                outer.append(values[0])
    if len(outer) < 10:
        raise ValueError("Too few outer field controls")
    local = []
    for angle in np.arange(12) * np.pi / 6:
        pos = offset(*center, 3 * np.cos(angle), 3 * np.sin(angle))
        local.append(sample(images, pos)[1])
    return {
        b: {
            "outer_controls": len(outer),
            "outer_at_or_above_target": int(
                np.sum(np.array(outer)[:, i] >= source_level[i])
            ),
            "median_level_target": round(float(source_level[i]), 7),
            "median_level_outer": round(float(np.median(np.array(outer)[:, i])), 7),
            "core_minus_annulus": round(float(contrast[i]), 7),
            "local_comparisons": len(local),
            "local_contrasts_at_or_above_target": int(
                np.sum(np.array(local)[:, i] >= contrast[i])
            ),
        }
        for i, b in enumerate("grz")
    }


def figure(images, center, host, nearest):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.8), layout="constrained")
    for i, (ax, k, extent, title) in enumerate(
        zip(
            axes,
            ["ls-dr10", "hsc-dr3", "hsc-dr3"],
            [39, 39, 8],
            ["DECam", "Subaru HSC", "Subaru HSC detail"],
        )
    ):
        data, wcs, scale = images[k]
        x, y = map(float, wcs.world_to_pixel_values(*center))
        z = np.median(data, axis=0)
        bg = np.nanpercentile(z, [5, 50])
        z = np.arcsinh((z - bg[1]) / max(1e-12, bg[1] - bg[0]))
        vmin, vmax = np.nanpercentile(z, [4, 99.95])
        ax.imshow(z, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        for (a, b), color, mark in [
            (center, "#ff6347", "+"),
            ((float(host["ra"]), float(host["dec"])), "#36d9e1", "x"),
        ]:
            px, py = wcs.world_to_pixel_values(a, b)
            ax.plot(
                px, py, marker=mark, color=color, markersize=10, markeredgewidth=1.2
            )
        if i == 2:
            px, py = wcs.world_to_pixel_values(
                float(nearest["ra"]), float(nearest["dec"])
            )
            ax.plot(px, py, "o", mfc="none", mec="#edc450", ms=6)
        ax.set(
            xlim=(x - extent / scale, x + extent / scale),
            ylim=(y - extent / scale, y + extent / scale),
            xticks=[],
            yticks=[],
            title=title,
        )
        length = 2 if i == 2 else 10
        x0 = x - 0.84 * extent / scale
        y0 = y - 0.8 * extent / scale
        ax.plot([x0, x0 + length / scale], [y0, y0], color="white", lw=2)
        ax.text(x0, y0 + extent * 0.04 / scale, f"{length}″", color="white", fontsize=8)
    fig.suptitle("NGC 1037 before the 2025 transient", fontsize=14)
    fig.savefig(HERE / "field.png", dpi=140)
    plt.close(fig)


def run(fetch=False, plot=False):
    DATA.mkdir(exist_ok=True)
    ra, dec, position = archive.centers()
    catalog = tap(
        "tractor",
        "ls_dr10.tractor",
        "ls_id,ra,dec,type,shape_r,shape_e1,shape_e2,flux_g,flux_r,flux_z,ref_cat,ref_id",
        "ra BETWEEN 39.479 AND 39.509 AND dec BETWEEN -1.865 AND -1.834",
        fetch,
    )
    sga = json.loads(
        archive.blob(
            "sga.json", VIEW + "sga2020/1/cat.json?" + urlencode(SKY_BOX), fetch
        )
    )
    desi = json.loads(
        archive.blob(
            "desi.json", VIEW + "desi-spec-dr1/1/cat.json?" + urlencode(SKY_BOX), fetch
        )
    )
    if len(sga["rd"]) != 1 or sga["pgc"] != [9973] or not desi["targetid"]:
        raise ValueError("Missing galaxy or spectra")
    targets = [int(i) for i in desi["targetid"]]
    spectra = tap(
        "redshifts",
        "desi_dr1.zpix",
        "targetid,mean_fiber_ra,mean_fiber_dec,z,zerr,zwarn,deltachi2,spectype,survey,program",
        "targetid IN (" + ",".join(str(i) for i in sorted(targets)) + ")",
        fetch,
    )
    if {int(s["targetid"]) for s in spectra} != set(targets):
        raise ValueError("Missing redshift rows")

    def sep(r):
        return archive.location(ra, dec, r["ra"], r["dec"])

    catalog.sort(key=sep)
    host = min(
        catalog, key=lambda r: archive.location(*sga["rd"][0], r["ra"], r["dec"])
    )
    if (
        host["ref_cat"] != "L3"
        or archive.location(*sga["rd"][0], host["ra"], host["dec"]) > 0.3
    ):
        raise ValueError("Ambiguous SGA center")
    spectra.sort(
        key=lambda r: archive.location(
            r["mean_fiber_ra"], r["mean_fiber_dec"], host["ra"], host["dec"]
        )
    )
    if (
        archive.location(
            spectra[0]["mean_fiber_ra"],
            spectra[0]["mean_fiber_dec"],
            host["ra"],
            host["dec"],
        )
        > 1
    ):
        raise ValueError("Spectrum misses nucleus")
    x, y = archive.pilot.offsets(ra, dec, float(host["ra"]), float(host["dec"]))
    images = {
        key: image(key, size, px, fetch)
        for key, size, px in [("ls-dr10", 384, 0.262), ("hsc-dr3", 600, 0.168)]
    }
    if plot:
        figure(images, (ra, dec), host, catalog[0])
    return {
        "id": NAME,
        "ra": ra,
        "dec": dec,
        **position,
        "sga": sga,
        "host_model": host,
        "offset_east_north_arcsec": [round(x, 3), round(y, 3)],
        "nuclear_distance_arcsec": round(sep(host), 3),
        "spectra": spectra,
        "nearest_3_tractor_models": [
            dict(separation_arcsec=round(sep(r), 3), **r) for r in catalog[:3]
        ],
        "local_texture": texture.test(
            images,
            (ra, dec),
            catalog,
            offset,
            archive.location,
            HERE / "texture.png" if plot else None,
        ),
        "pre_event_light": {
            k: image_test(im, (ra, dec), host, catalog) for k, im in images.items()
        },
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--plot", action="store_true")
    a = p.parse_args()
    results = run(a.fetch, a.plot)
    path = HERE / "results.json"
    if a.fetch:
        archive.pilot.json_write(path, results)
    else:
        archive.compass.check_replay(results, json.loads(path.read_text()), str(path))
    print(
        NAME,
        results["nuclear_distance_arcsec"],
        "arcsec from nucleus across",
        results["nights"],
        "nights",
    )
