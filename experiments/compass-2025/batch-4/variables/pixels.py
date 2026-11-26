"""Forced aperture flux at a fixed sky position and field-star calendar."""

import argparse
from collections import defaultdict
import csv
from datetime import date, datetime, timezone
import gzip, hashlib, io, json, math
from pathlib import Path
from statistics import median

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import requests
import time
from scipy.stats import spearmanr

from field import SRC, offset
from measure import HERE, BATCH, baseline, load, pilot

NAME = "ZTF17aabutiy"
MAN = SRC / NAME / "pixel-manifest.json"
TILES = SRC / NAME / "pixel-tiles.npy.gz"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def inputs():
    case, _, _ = load(NAME)
    field = json.loads((HERE / "field.json").read_text())[NAME]
    direct = json.loads((BATCH / "results/direct.json").read_text())[NAME]
    refs = {}
    jobs = []
    catalog = list(csv.DictReader(gzip.open(SRC / NAME / "field.csv.gz", "rt")))
    for b, p in field["primary"].items():
        target = direct["chosen_series"]["1" if b == "g" else "2"]
        same = [
            v
            for v in catalog
            if all(v[k] == target[k] for k in ("fid", "field", "ccdid", "qid"))
        ]
        cands = []
        for v in same:
            if not v["meanmag"] or not (
                10 <= offset(v, case) <= 45 and 14.5 <= float(v["meanmag"]) <= 18
            ):
                continue
            if any(
                a["oid"] != v["oid"]
                and math.hypot(
                    *pilot.offsets(
                        float(a["ra"]), float(a["dec"]), float(v["ra"]), float(v["dec"])
                    )
                )
                <= 5
                for a in same
            ):
                continue
            cands.append(v)
        refs[b] = sorted(cands, key=lambda v: (float(v["meanmag"]), int(v["oid"])))[:5]
        if len(refs[b]) < 2:
            raise ValueError("Insufficient local calibrators")
        exposures = {
            r["expid"]: r
            for r in csv.DictReader(
                gzip.open(SRC / NAME / (b + "-controls.csv.gz"), "rt")
            )
        }
        for span, m in p["coverage"].items():
            for day, v in sorted(m["field_calendar"].items()):
                r = exposures[v["expid"]]
                serial = r["filefracday"]
                field_num = int(r["field"])
                ccd = int(r["ccdid"], 0)
                q = int(r["qid"], 0)
                path = f"{serial[:4]}/{serial[4:8]}/{serial[8:]}/ztf_{serial}_{field_num:06d}_{r['filtercode']}_c{ccd:02d}_o_q{q}_sciimg.fits"
                url = f"https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{path}?center={case['ra']},{case['dec']}&size=120arcsec&gzip=false"
                jobs.append(
                    {
                        "date": day,
                        "span": span,
                        "band": b,
                        "expid": v["expid"],
                        "catalog_mag": v["target_mag"],
                        "field": field_num,
                        "ccd": ccd,
                        "qid": q,
                        "url": url,
                    }
                )
    return case, refs, jobs


def fetch_one(job, case, refs):
    response = requests.get(job["url"], timeout=60)
    response.raise_for_status()
    raw = response.content
    with fits.open(io.BytesIO(raw)) as hdus:
        h = hdus[0].header
        image = hdus[0].data
        if any(
            int(h[k]) != v
            for k, v in [
                ("EXPID", int(job["expid"])),
                ("FIELDID", job["field"]),
                ("CCDID", job["ccd"]),
                ("QID", job["qid"]),
                ("FILTERID", 1 if job["band"] == "g" else 2),
            ]
        ):
            raise ValueError("Image metadata disagree with field calendar")
        wcs = WCS(h).celestial
        stars = [(case["ra"], case["dec"])] + [
            (float(s["ra"]), float(s["dec"])) for s in refs[job["band"]]
        ]
        tiles = []
        centers = []
        for ra, dec in stars:
            x, y = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
            px, py = round(float(x)), round(float(y))
            tile = image[py - 9 : py + 10, px - 9 : px + 10].astype("<f4")
            if tile.shape != (19, 19):
                raise ValueError("Source falls too near edge of cutout")
            tiles.append(tile)
            centers.append([float(x - px + 9), float(y - py + 9)])
        meta = {
            "raw_fits_sha256": sha(raw),
            "seeing_arcsec": float(h["SEEING"]),
            "centers": centers,
        }
    return meta, tiles


def capture(case, refs, jobs):
    saved = {}
    if MAN.exists() and TILES.exists():
        prior = json.loads(MAN.read_text())
        raw = gzip.decompress(TILES.read_bytes())
        if sha(raw) != prior["tile_sha256"]:
            raise ValueError("Previous pixel capture failed checksum")
        array = np.load(io.BytesIO(raw), allow_pickle=False)
        for old in prior["images"]:
            if "tile_index" in old:
                saved[old["url"]] = (old, list(array[old["tile_index"]]))

    def task(job):
        if job["url"] in saved:
            return saved[job["url"]]

        for attempt in range(4):
            try:
                return fetch_one(job, case, refs)
            except requests.RequestException as exc:
                if attempt == 3:
                    return {"error": str(exc)}, []
                time.sleep(2 ** (attempt + 1))

    snapshots = [task(job) for job in jobs]
    metas = []
    tiles = []
    for job, (meta, arr) in zip(jobs, snapshots):
        if arr:
            meta["tile_index"] = len(tiles)
            tiles.append(arr)
        metas.append({**job, **meta})
    if not tiles:
        raise ValueError("No science exposures acquired")
    stream = io.BytesIO()
    np.save(stream, np.array(tiles, dtype="<f4"), allow_pickle=False)
    raw = stream.getvalue()
    TILES.write_bytes(gzip.compress(raw, mtime=0))
    pilot.json_write(
        MAN,
        {
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "ref_ids": {b: [r["oid"] for r in v] for b, v in refs.items()},
            "tile_sha256": sha(raw),
            "images": metas,
        },
    )


def counts(tile, center, radius):
    y, x = np.indices(tile.shape)
    dist = np.hypot(x - center[0], y - center[1])
    ap = tile[dist <= radius]
    sky = tile[(dist >= 5) & (dist <= 8)]
    sky = sky[np.isfinite(sky)]
    if not np.all(np.isfinite(ap)) or len(sky) < 30:
        return np.nan
    m = np.median(sky)
    sigma = 1.4826 * np.median(np.abs(sky - m))
    if sigma > 0:
        sky = sky[np.abs(sky - m) <= 3 * sigma]
    return float(np.sum(ap.astype(float) - np.median(sky)))


def flux_rows(refs, arrays, images, radius):
    out = []
    for job in images:
        if "tile_index" not in job:
            continue
        stars = refs[job["band"]]
        fluxes = np.array([3631000 * 10 ** (-0.4 * float(r["meanmag"])) for r in stars])
        sums = np.array(
            [
                counts(a, c, radius)
                for a, c in zip(arrays[job["tile_index"]], job["centers"], strict=True)
            ]
        )
        scales = [f / s for f, s in zip(fluxes, sums[1:]) if np.isfinite(s) and s > 0]
        if len(scales) < 2 or not np.isfinite(sums[0]):
            continue
        calibrated = []
        for i, s in enumerate(sums[1:]):
            other = [
                f / v
                for j, (f, v) in enumerate(zip(fluxes, sums[1:]))
                if j != i and np.isfinite(v) and v > 0
            ]
            calibrated.append(
                float(s * median(other)) if len(other) >= 2 and np.isfinite(s) else None
            )
        out.append(
            {
                "date": job["date"],
                "span": job["span"],
                "band": job["band"],
                "expid": job["expid"],
                "mjy": round(float(sums[0] * median(scales)), 7),
                "catalog_mag": job["catalog_mag"],
                "seeing": job["seeing_arcsec"],
                "references_mjy": calibrated,
            }
        )
    return out


def measures(rows, edge, anchor, bright, record):
    valid = [r for r in rows if r["catalog_mag"] is not None]
    mag = [r["catalog_mag"] for r in valid]
    flux = [r["mjy"] for r in valid]
    rho = (
        round(float(spearmanr(-np.array(mag), flux).statistic), 5)
        if len(mag) >= 3
        else None
    )
    dd = {date.fromisoformat(r["date"]): (0, r["mjy"]) for r in rows}

    def contrast(curve):
        return {
            k.replace("mag", "mjy"): v
            for k, v in baseline.change(baseline.bin_curve(curve, anchor), edge).items()
        }

    out = {
        "all_dates": len(rows),
        "catalog_dates": len(valid),
        "catalog_rho": rho,
        "change": contrast(dd),
        "windows": {},
    }
    for span in ("before", "after"):
        vs = [r for r in rows if r["span"] == span]
        vals = [r["mjy"] for r in vs]
        if not vals:
            continue
        out["windows"][span] = {
            "n": len(vs),
            "missing_catalog": sum(r["catalog_mag"] is None for r in vs),
            "negative_flux": sum(v < 0 for v in vals),
            "q20_q50_q80_mjy": list(np.round(np.quantile(vals, [0.2, 0.5, 0.8]), 6)),
            "bright_fraction": round(np.mean(np.array(vals) >= bright), 5),
            "median_seeing": median(r["seeing"] for r in vs),
        }
    if record:
        # Compute paired comparisons within exactly the same night and date-bin scheme.
        comparisons = {
            "forced_all": dd,
            "forced_catalogued": {
                date.fromisoformat(r["date"]): (0, r["mjy"]) for r in valid
            },
            "catalogued_flux": {
                date.fromisoformat(r["date"]): (
                    0,
                    3631000 * 10 ** (-0.4 * r["catalog_mag"]),
                )
                for r in valid
            },
        }
        out["matched_bin_medians"] = {}
        for label, data in comparisons.items():
            bins = baseline.bin_curve(data, anchor)
            early = [v for t, v in bins.items() if edge - 180 <= t < edge]
            later = [v for t, v in bins.items() if edge <= t < edge + 180]
            if min(len(early), len(later)) < 3:
                out["matched_bin_medians"][label] = {
                    "before_bins": len(early),
                    "after_bins": len(later),
                }
                continue
            a, z = median(early), median(later)
            out["matched_bin_medians"][label] = {
                "before_bins": len(early),
                "after_bins": len(later),
                "median_before_mjy": round(a, 6),
                "median_after_mjy": round(z, 6),
                "equivalent_dimming_mag": (
                    round(2.5 * math.log10(a / z), 5) if min(a, z) > 0 else None
                ),
            }
        out["dates"] = rows
        out["reference_changes_mjy"] = [
            contrast(
                {
                    date.fromisoformat(r["date"]): (0, r["references_mjy"][i])
                    for r in rows
                    if r["references_mjy"][i] is not None
                }
            ).get("after_minus_before_mjy")
            for i in range(len(rows[0]["references_mjy"]) if rows else 0)
        ]
    return out


def run(fetch):
    c, refs, jobs = inputs()
    if fetch and (
        not MAN.exists()
        or any("error" in r for r in json.loads(MAN.read_text())["images"])
    ):
        capture(c, refs, jobs)
    m = json.loads(MAN.read_text())
    if m["ref_ids"] != {b: [r["oid"] for r in v] for b, v in refs.items()} or len(
        jobs
    ) != len(m["images"]):
        raise ValueError("Science selection changed")
    for job, image in zip(jobs, m["images"], strict=True):
        if not all(image[k] == v for k, v in job.items()):
            raise ValueError("Exposure calendar changed")
    raw = gzip.decompress(TILES.read_bytes())
    if sha(raw) != m["tile_sha256"]:
        raise ValueError("Pixel tiles checksum mismatch")
    array = np.load(io.BytesIO(raw), allow_pickle=False)
    edge = json.loads((HERE / "results.json").read_text())[NAME]["state"][
        "boundary_mjd"
    ]
    field = json.loads((HERE / "field.json").read_text())[NAME]["primary"]
    results = {}
    for radius in (2.2, 1.6, 3.2):
        rows = flux_rows(refs, array, m["images"], radius)
        results[str(radius)] = {}
        for b in "rg":
            vals = [r for r in rows if r["band"] == b]
            bright = 3631000 * 10 ** (
                -0.4 * field[b]["uncorrected"]["outside_q20_q50_q80"][0]
            )
            results[str(radius)][b] = measures(
                vals, edge, date.fromisoformat(c["date"]), bright, radius == 2.2
            )
    return {
        "requested": len(jobs),
        "failed_fetches": [v for v in m["images"] if "error" in v],
        "apertures_pixels": results,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "pixels.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
