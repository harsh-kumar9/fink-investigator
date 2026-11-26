"""Replay the Legacy Surveys counterpart, redshifts and matched-flux stellar controls."""

import argparse
import io
import json
import csv
import gzip
import hashlib
from datetime import datetime, timezone
import warnings
from urllib.parse import urlencode

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates

import check as ps1

LS = "https://www.legacysurvey.org/viewer/fits-cutout?"
NERSC = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/"
VIEWER = "https://www.legacysurvey.org/viewer/"
PIX = 0.262
SIZE = 320
GRZ = "grz"


def tractor(release, fetch):
    url = NERSC + f"dr{release}/south/tractor/250/tractor-2500p297.fits"
    prefix = ps1.SOURCES / f"legacy-dr{release}"
    receipt, data_file = prefix.with_suffix(".json"), prefix.with_suffix(".csv.gz")
    if fetch and not receipt.exists():
        raw = ps1.pilot.curl(url, [])
        columns = (
            "release brickid objid ra dec type shape_r shape_r_ivar "
            "flux_g flux_r flux_z flux_ivar_g flux_ivar_r flux_ivar_z maskbits "
            "fracflux_r mjd_min mjd_max nobs_g nobs_r nobs_i nobs_z"
        ).split()
        with fits.open(io.BytesIO(raw)) as hdus:
            data = hdus[1].data
            if release == 9:
                columns.remove("nobs_i")
            use = (
                data["brick_primary"]
                & (data["ra"] > 250.013)
                & (data["ra"] < 250.043)
                & (data["dec"] > 29.750)
                & (data["dec"] < 29.777)
            )
            records = []
            for row in data[use]:
                rec = {c: row[c] for c in columns}
                rec["ls_id"] = (
                    (int(row["release"]) << 40)
                    + (int(row["brickid"]) << (16 if release == 9 else 20))
                    + int(row["objid"])
                )
                rec["dchisq_1"], rec["dchisq_2"] = row["dchisq"][:2]
                records.append(rec)
        csv_text = io.StringIO()
        writer = csv.DictWriter(
            csv_text, fieldnames=[*columns, "ls_id", "dchisq_1", "dchisq_2"]
        )
        writer.writeheader()
        writer.writerows(records)
        subset = csv_text.getvalue().encode()
        data_file.write_bytes(gzip.compress(subset, mtime=0))
        ps1.pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256_full_fits": hashlib.sha256(raw).hexdigest(),
                "sha256_sky_subset_csv": hashlib.sha256(subset).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    subset = gzip.decompress(data_file.read_bytes())
    info = json.loads(receipt.read_text())
    if (
        info["url"] != url
        or info["sha256_sky_subset_csv"] != hashlib.sha256(subset).hexdigest()
    ):
        raise ValueError("Unexpected Tractor source or digest")
    return list(csv.DictReader(io.StringIO(subset.decode())))


def capture(fetch):
    params = {
        "ra": ps1.RA,
        "dec": ps1.DEC,
        "layer": "ls-dr10",
        "pixscale": PIX,
        "size": SIZE,
        "bands": GRZ,
    }
    raw = ps1.blob("legacy-grz.fits", LS + urlencode(params), fetch)
    with fits.open(io.BytesIO(raw)) as hdus:
        data, hdr = hdus[0].data.copy(), hdus[0].header
        if data.shape != (3, SIZE, SIZE) or hdr["BANDS"] != GRZ:
            raise ValueError("Wrong Legacy Surveys image")
        wcs = WCS(hdr).celestial
        x, y = wcs.world_to_pixel_values(ps1.RA, ps1.DEC)
        if max(abs(x - (SIZE - 1) / 2), abs(y - (SIZE - 1) / 2)) > 2:
            raise ValueError("Wrong Legacy Surveys pointing")
    return tractor(10, fetch), tractor(9, fetch), data, wcs


def measure(data, wcs, source, star, neighbors):
    from numpy import round as rnd

    sx, sy = map(
        float, wcs.world_to_pixel_values(float(source["ra"]), float(source["dec"]))
    )
    ux, uy = map(
        float, wcs.world_to_pixel_values(float(star["ra"]), float(star["dec"]))
    )
    radius = 20
    gy, gx = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    phase_x, phase_y = sx - rnd(sx), sy - rnd(sy)
    r = np.hypot(gx - phase_x, gy - phase_y) * PIX
    sky_mask = (r >= 3.3) & (r <= 5)
    apertures = [(r <= d) for d in (0.8, 1.8)]

    def crop(x, y):
        return data[:, int(rnd(y)) + gy, int(rnd(x)) + gx].astype(float)

    def sky_sub(array):
        return array - np.median(array[:, sky_mask], axis=1)[:, None, None]

    def fluxes(array):
        return np.array([array[:, mask].sum(axis=1) for mask in apertures])

    obj = sky_sub(crop(sx, sy))
    observed = fluxes(obj)
    psf = np.array(
        [
            map_coordinates(im, [uy + gy - phase_y, ux + gx - phase_x], order=1)
            for im in data
        ]
    )
    psf = sky_sub(psf)
    matched = psf * (observed[1] / fluxes(psf)[1])[:, None, None]
    positions = np.array(
        [wcs.world_to_pixel_values(float(s["ra"]), float(s["dec"])) for s in neighbors]
    )
    injections = []
    examples = []
    for y in range(radius + 3, SIZE - radius, 40):
        for x in range(radius + 3, SIZE - radius, 40):
            if (
                np.hypot(
                    positions[:, 0] - (x + phase_x), positions[:, 1] - (y + phase_y)
                ).min()
                * PIX
                <= 6
            ):
                continue
            local = crop(x, y)
            if not np.all(np.isfinite(local)) or np.any(
                np.std(local, axis=(1, 2)) <= 0
            ):
                continue
            copy = sky_sub(local + matched)
            inner, outer = fluxes(copy)
            if np.any(outer <= 0):
                raise ValueError("Unexpected nonpositive reference flux")
            injections.append(inner / outer)
            examples.append(copy)
    if len(injections) < 5:
        raise ValueError("Insufficient independent sky positions")
    ratios = np.array(injections)
    measured = observed[0] / observed[1]
    summary = {}
    for i, b in enumerate(GRZ):
        summary[b] = dict(
            source=round(float(measured[i]), 4),
            n=len(ratios),
            controls_percentiles_5_50_95=np.round(
                np.percentile(ratios[:, i], [5, 50, 95]), 4
            ).tolist(),
            controls_at_or_below_target=int(np.sum(ratios[:, i] <= measured[i])),
        )
    # Use the injection nearest the median r concentration for the display.
    ex = np.argmin(abs(ratios[:, 1] - np.median(ratios[:, 1])))
    return summary, obj, examples[ex]


def redshifts(dr9, dr10, fetch):
    # Locate this exact RELEASE, BRICKID, OBJID in the published 10.1 photo-z sweep.
    url = NERSC + "dr10/south/sweep/10.1-photo-z/sweep-250p025-255p030-pz.fits"
    row_index = 1832594
    headers = io.BytesIO(ps1.blob("photo-z-headers.fits", url, fetch, (0, 8639)))
    primary = fits.Header.fromfile(headers)
    tabhead = fits.Header.fromfile(headers)
    row_size = tabhead["NAXIS1"]
    if not (0 <= row_index < tabhead["NAXIS2"]) or headers.tell() != 8640:
        raise ValueError("Unexpected photo-z sweep layout")
    start = headers.tell() + row_size * row_index
    raw_row = ps1.blob("photo-z-row.bin", url, fetch, (start, start + row_size - 1))
    tabhead["NAXIS2"] = 1
    one_row_table = (
        primary.tostring().encode()
        + tabhead.tostring().encode()
        + raw_row
        + bytes((-row_size) % 2880)
    )
    with fits.open(io.BytesIO(one_row_table)) as hdu:
        data = hdu[1].data
        phot = {}
        for name in data.names:
            v = data[name][0]
            phot[name] = v.item() if isinstance(v, np.generic) else v
    if [phot[k] for k in ("RELEASE", "BRICKID", "OBJID")] != [
        int(dr10[k]) for k in ("release", "brickid", "objid")
    ]:
        raise ValueError("Photo-z row belongs to a different source")
    sky_box = {
        "ralo": "250.02629",
        "rahi": "250.02950",
        "declo": "29.76188",
        "dechi": "29.76467",
    }
    old = json.loads(
        ps1.blob(
            "photo-z-dr9.json",
            VIEWER + "photoz-dr9/1/cat.json?" + urlencode(sky_box),
            fetch,
        )
    )
    if len(old["rd"]) != 1 or ps1.location(*old["rd"][0], dr9["ra"], dr9["dec"]) > 0.01:
        raise ValueError("Ambiguous DR9 identity")
    desi = json.loads(
        ps1.blob(
            "desi-dr1.json",
            VIEWER + "desi-spec-dr1/1/cat.json?" + urlencode(sky_box),
            fetch,
        )
    )
    if not isinstance(desi["rd"], list):
        raise ValueError("Missing DESI positions")
    sql = (
        "SELECT s.specObjID,s.ra,s.dec,s.z,s.zWarning FROM "
        "dbo.fGetNearbySpecObjEq(250.02789,29.76327,0.0833334) n JOIN SpecObj s ON n.specObjID=s.specObjID"
    )
    sdss_url = (
        "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch?"
        + urlencode({"cmd": sql, "format": "csv"})
    )
    lines = ps1.blob("sdss-spectra.csv", sdss_url, fetch).decode().splitlines()
    if lines[:2] != ["#Table1", "specObjID,ra,dec,z,zWarning"]:
        raise ValueError("Missing SDSS spectroscopy query")
    return {
        "dr10": phot,
        "dr9": old,
        "desi_dr1_spectrum_positions_in_box": desi["rd"],
        "sdss_dr18_primary_spectra_within_5arcsec": list(csv.DictReader(lines[1:])),
    }


def run(fetch=False, plot=False):
    ra, dec, _ = ps1.centers()
    near, old, data, wcs = capture(fetch)
    sep = lambda row: ps1.location(row["ra"], row["dec"], ra, dec)
    near.sort(key=sep)
    old.sort(key=sep)
    if not near or not old or sep(near[0]) > 1 or sep(old[0]) > 1:
        raise ValueError("Missing optical source")
    ps, stacks, primary = ps1.run(False)
    closest_gaia = ps["gaia_controls"]["r"][0]["objID"]
    control = next(r for r in primary if r["objID"] == closest_gaia)
    nearby_star = min(
        near,
        key=lambda s: ps1.location(
            s["ra"], s["dec"], control["raStack"], control["decStack"]
        ),
    )
    if (
        ps1.location(
            nearby_star["ra"],
            nearby_star["dec"],
            control["raStack"],
            control["decStack"],
        )
        > 0.5
        or nearby_star["type"] != "PSF"
    ):
        raise ValueError("Incorrect stellar control")
    stats, obj, reference = measure(data, wcs, near[0], nearby_star, near)
    results = {
        "dr10": near[0],
        "dr9": old[0],
        "offset_arcsec": round(sep(near[0]), 3),
        "next_dr10_separation_arcsec": round(sep(near[1]), 3),
        "profile_star": nearby_star["ls_id"],
        "concentration_0.8_over_1.8arcsec": stats,
        "redshifts": redshifts(old[0], near[0], fetch),
    }
    if plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 3, figsize=(6, 6), layout="constrained")
        for i, b in enumerate(GRZ):
            for ax, values in zip(axes[1:, i], [obj[i], reference[i]]):
                cut = values[8:-8, 8:-8]
                scale = np.percentile(obj[i][8:-8, 8:-8], 99)
                ax.imshow(
                    cut, origin="lower", cmap="gray", vmin=-scale * 0.12, vmax=scale
                )
            ax = axes[0, i]
            band = ps1.BANDS[i]
            img, pwcs = stacks[band]
            xx, yy = map(round, map(float, pwcs.world_to_pixel_values(ra, dec)))
            cut = img[yy - 13 : yy + 14, xx - 13 : xx + 14]
            bg = np.nanpercentile(cut, [5, 99])
            ax.imshow(cut, origin="lower", cmap="gray", vmin=bg[0], vmax=bg[1])
            for j, label in enumerate(
                ["PS1 " + band, "DR10 " + b, "DR10 " + b + " control"]
            ):
                axes[j, i].set_title(label, fontsize=10)
        for ax in axes.flat:
            ax.set(xticks=[], yticks=[])
        fig.suptitle("Pre-event counterpart and matched-flux field star", fontsize=13)
        fig.savefig(ps1.HERE / "profiles.png", dpi=125)
        plt.close(fig)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FITSFixedWarning, message="PC00.*")
        result = run(args.fetch, args.plot)
    saved = ps1.HERE / "legacy.json"
    if args.fetch:
        ps1.pilot.json_write(saved, result)
    else:
        ps1.compass.check_replay(result, json.loads(saved.read_text()), str(saved))
    print(json.dumps(result, indent=2))
