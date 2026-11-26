"""Compare the position of ZTF25aanlbls with pre-event optical catalogs and stacks."""

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
YEAR = HERE.parents[2]
sys.path.insert(0, str(YEAR.parents[1]))
import compass
import compass_archive
import compass_direct
import pilot

NAME = "ZTF25aanlbls"
RA, DEC = 250.0278614, 29.763246
BANDS = "gri"
CONE = {"ra": RA, "dec": DEC, "radius": 60 / 3600}
PS1 = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/"
IMAGES = "https://ps1images.stsci.edu/cgi-bin/"
SOURCES = HERE / "sources"


def catalog(table, columns, fetch):
    params = {**CONE, "pagesize": 1000, "columns": "[" + ",".join(columns) + "]"}
    rows = compass_direct.capture(
        SOURCES, table, PS1 + table + ".csv?" + urlencode(params), set(columns), fetch
    )
    if len(rows) >= 1000:
        raise ValueError("Truncated cone")
    return rows


def blob(key, url, fetch, byte_range=None):
    path, receipt = SOURCES / (key + ".gz"), SOURCES / (key + ".json")
    if fetch and not receipt.exists():
        extra = ["--range", f"{byte_range[0]}-{byte_range[1]}"] if byte_range else []
        raw = pilot.curl(url, extra)
        if byte_range and len(raw) != byte_range[1] - byte_range[0] + 1:
            raise ValueError("Server did not honor the requested byte range")
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                **({"byte_range": list(byte_range)} if byte_range else {}),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    raw = gzip.decompress(path.read_bytes())
    info = json.loads(receipt.read_text())
    if (
        info.get("byte_range") != (list(byte_range) if byte_range else None)
        or info["url"] != url
        or info["sha256_uncompressed"] != hashlib.sha256(raw).hexdigest()
    ):
        raise ValueError("Source request or checksum changed: " + key)
    return raw


def load_stacks(fetch):
    common = {"ra": RA, "dec": DEC}
    url = (
        IMAGES
        + "ps1filenames.py?"
        + urlencode({**common, "filters": BANDS, "type": "stack"})
    )
    lines = blob("stacks.txt", url, fetch).decode().splitlines()
    stacks = {}
    for line in lines[1:]:
        row = dict(zip(lines[0].split(), line.split(), strict=True))
        band = row["filter"]
        if band in stacks or row["badflag"] != "0" or row["type"] != "stack":
            raise ValueError("Unexpected skycell row")
        params = {**common, "size": 512, "format": "fits", "red": row["filename"]}
        raw = blob(
            "stack-" + band + ".fits",
            IMAGES + "fitscut.cgi?" + urlencode(params),
            fetch,
        )
        with fits.open(io.BytesIO(raw)) as hdus:
            data, head = hdus[0].data, hdus[0].header
            wcs = WCS(head, naxis=2)
            x, y = wcs.world_to_pixel_values(RA, DEC)
            if data.shape != (512, 512) or max(abs(x - 255.5), abs(y - 255.5)) > 2:
                raise ValueError("Invalid sky coverage")
            stacks[band] = (data.astype(float), wcs)
    if set(stacks) != set(BANDS):
        raise ValueError("Missing stack")
    return stacks


def location(ra, dec, ra0, dec0):
    return math.hypot(*pilot.offsets(float(ra), float(dec), float(ra0), float(dec0)))


def centers():
    _, histories, _ = compass_archive.load(YEAR / "snapshot")
    nights = defaultdict(list)
    for r in histories[NAME]:
        if compass.quality(r) and compass.sign(r) == 1:
            nights[str(pilot.day_of(r["i:jd"]))].append([r["i:ra"], r["i:dec"]])
    xy = np.array([np.median(v, axis=0) for _, v in sorted(nights.items())])
    ra, dec = np.median(xy, axis=0)
    radii = [location(a, b, ra, dec) for a, b in xy]
    return (
        float(ra),
        float(dec),
        {
            "quality_alerts": sum(map(len, nights.values())),
            "nights": len(nights),
            "first": min(nights),
            "last": max(nights),
            "median_nightly_radius_arcsec": round(float(np.median(radii)), 3),
            "p90_nightly_radius_arcsec": round(float(np.percentile(radii, 90)), 3),
        },
    )


def phot(row, band):
    keys = ("PSFMag", "PSFMagErr", "KronMag", "KronMagErr", "psfQfPerfect")
    psf, e1, kron, e2, q = [float(row[band + k]) for k in keys]
    if min(psf, kron, e1, e2) <= 0 or not (0 < q <= 1):
        return None
    return dict(
        psf=round(psf, 3),
        kron=round(kron, 3),
        difference=round(psf - kron, 3),
        psf_err=round(e1, 3),
        kron_err=round(e2, 3),
        unmasked_fraction=round(q, 3),
    )


def gaia_match(row, stars):
    epochs = [float(row[b + "Epoch"]) for b in BANDS if float(row[b + "Epoch"]) > 0]
    if not epochs:
        return None
    year = (np.median(epochs) - 51544.5) / 365.25 + 2000
    candidates = []
    for s in stars:
        x, y = pilot.offsets(
            float(row["raStack"]),
            float(row["decStack"]),
            float(s["ra"]),
            float(s["dec"]),
        )
        if s["pmra"] and s["pmdec"]:
            dt = year - float(s["ref_epoch"])
            x -= float(s["pmra"]) * dt / 1000
            y -= float(s["pmdec"]) * dt / 1000
        candidates.append((math.hypot(x, y), s["source_id"]))
    d, name = min(candidates)
    return name if d <= 0.7 else None


def run(fetch=False):
    ra, dec, astrometry = centers()
    mean = catalog(
        "mean",
        [
            "objID",
            "raMean",
            "decMean",
            "epochMean",
            "nDetections",
            *[b + "MeanPSFMag" for b in BANDS],
        ],
        fetch,
    )
    fields = [
        "objID",
        "raStack",
        "decStack",
        "nDetections",
        "primaryDetection",
        "bestDetection",
    ]
    fields += [
        b + k
        for b in BANDS
        for k in (
            "ra",
            "dec",
            "Epoch",
            "PSFMag",
            "PSFMagErr",
            "KronMag",
            "KronMagErr",
            "psfQfPerfect",
        )
    ]
    stack = catalog("stack", fields, fetch)
    primary = [r for r in stack if r["primaryDetection"] == "1"]
    stacks = load_stacks(fetch)
    with gzip.open(HERE.parent / NAME / "gaia.csv.gz", "rt") as h:
        stars = list(csv.DictReader(h))

    def distance(row):
        return location(row["raStack"], row["decStack"], ra, dec)

    primary.sort(key=distance)
    matches = []
    for row in primary:
        sep = distance(row)
        if sep > 10:
            break
        matches.append(
            dict(
                objID=row["objID"],
                ra=float(row["raStack"]),
                dec=float(row["decStack"]),
                separation_arcsec=round(sep, 3),
                bands={b: phot(row, b) for b in BANDS},
            )
        )
    controls = [(r, gaia_match(r, stars)) for r in primary if 8 <= distance(r) <= 60]
    stellar = {}
    for b in BANDS:
        brightest = next((m["bands"][b]["psf"] for m in matches if m["bands"][b]), 21)
        good = [(r, gid, phot(r, b)) for r, gid in controls if gid and phot(r, b)]
        good = [(r, gid, p) for r, gid, p in good if 15 <= p["psf"] <= 21]
        good.sort(key=lambda v: (abs(v[2]["psf"] - brightest), distance(v[0])))
        stellar[b] = [
            dict(
                objID=r["objID"], gaia=gid, separation_arcsec=round(distance(r), 3), **p
            )
            for r, gid, p in good[:10]
        ]
    near_mean = sorted(
        [(location(r["raMean"], r["decMean"], ra, dec), r) for r in mean],
        key=lambda pair: pair[0],
    )
    answer = {
        "id": NAME,
        "alert_center": {"ra": ra, "dec": dec, **astrometry},
        "mean_rows": len(mean),
        "stack_rows": len(stack),
        "primary_rows": len(primary),
        "primary_within_10arcsec": matches,
        "nearest_mean": [
            [round(d, 3), r["objID"], int(r["nDetections"])] for d, r in near_mean[:3]
        ],
        "gaia_controls": stellar,
    }
    return answer, stacks, primary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    result, stacks, primary = run(args.fetch)
    expected = HERE / "results.json"
    if args.fetch:
        pilot.json_write(expected, result)
    else:
        compass.check_replay(result, json.loads(expected.read_text()), str(expected))
    print(json.dumps(result, indent=2))
