"""Compare PS1 DR2 detections with Sesar's fixed four-band RR Lyrae model."""

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import io
import json
import math
import tarfile
from urllib.parse import urlencode

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

from priors import HERE, POSITIONS, archive, digest, literature

iers.conf.auto_download = False
NAME = "MGAB-V488"
RA, DEC = POSITIONS[NAME]
BANDS = "griz"
MAST = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/"
TEMPLATES = (
    "https://raw.githubusercontent.com/astroML/astroML-data/main/datasets/RRLyr_ugriz_templates.tar.gz"
)


def table(endpoint, label, params, cols, fetch):
    url = (
        MAST
        + endpoint
        + ".csv?"
        + urlencode({**params, "columns": "[" + ",".join(cols) + "]", "pagesize": 1000})
    )
    raw = archive.source("ps1-" + label + ".csv", url, fetch)
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if reader.fieldnames != cols:
        raise ValueError(f"Unexpected {endpoint} response")
    rows = list(reader)
    if len(rows) >= 1000:
        raise ValueError("MAST table truncated")
    return rows


def templates(model, fetch):
    receipt = HERE / "ps1-templates.json"
    if fetch and not receipt.exists():
        raw = archive.pilot.curl(TEMPLATES, [])
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as t:
            subset = {b: t.extractfile(model["T" + b] + b + ".dat").read().decode() for b in BANDS}
        data = json.dumps(subset, sort_keys=True).encode()
        (HERE / "ps1-templates.gz").write_bytes(gzip.compress(data, mtime=0))
        archive.pilot.json_write(
            receipt,
            {
                "url": TEMPLATES,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "sha256_full_response": digest(raw),
                "sha256_uncompressed": digest(data),
                "members": [model["T" + b] + b + ".dat" for b in BANDS],
            },
        )
    data = json.loads(archive.source("ps1-templates", TEMPLATES, False))
    return {b: np.loadtxt(io.StringIO(v)).T for b, v in data.items()}


def sep(ra, dec):
    return math.hypot(*archive.pilot.offsets(ra, dec, RA, DEC))


def collect(rows, snr, clean):
    groups = defaultdict(list)
    for r in rows:
        t, f, e, q = [float(r[k]) for k in ["obsTime", "psfFlux", "psfFluxErr", "psfQfPerfect"]]
        band = "grizy"[int(r["filterID"]) - 1]
        if (
            band in BANDS
            and all(map(math.isfinite, [t, f, e, q]))
            and e > 0
            and f / e >= snr
            and q >= clean
            and sep(float(r["ra"]), float(r["dec"])) <= 1.2
        ):
            groups[(band, math.floor(t))].append(
                [t, -2.5 * math.log10(f / 3631), 2.5 * e / (math.log(10) * f)]
            )
    return groups


def analyze(rows, shapes, model, snr, clean):
    groups = collect(rows, snr, clean)
    keys = sorted(groups, key=lambda key: (key[1], key[0]))
    dates = np.array([d for b, d in keys])
    bands = np.array([b for b, d in keys])
    n = len(keys)
    obs_days = np.unique(dates)
    if n < 1:
        return {"retained_filter_dates": 0}
    y = np.array([np.mean([r[1] for r in groups[k]]) for k in keys])
    # Use the largest individual error in each nightly mean, allowing full intra-night correlation.
    e = np.array([max(r[2] for r in groups[k]) for k in keys])
    bfit = np.unique(bands, return_counts=True)
    eligible = bfit[0][bfit[1] >= 3]
    selected = np.isin(bands, eligible)
    flat = np.array([r[0] for k in keys for r in groups[k]])
    times = Time(
        flat,
        format="mjd",
        scale="utc",
        location=EarthLocation.from_geodetic(-156.2567 * u.deg, 20.7075 * u.deg, 3055 * u.m),
    )
    correction = times.light_travel_time(SkyCoord(RA * u.deg, DEC * u.deg), kind="heliocentric").to_value(
        u.day
    )
    phase = ((flat + 0.5 + correction) / float(model["Per"]) + float(model["phi0"])) % 1
    exposure_group = np.repeat(np.arange(n), [len(groups[k]) for k in keys])
    exposure_band = bands[exposure_group]

    def prediction(shift):
        predicted = np.empty(len(flat))
        for b in np.unique(exposure_band):
            idx = exposure_band == b
            xp, yp = shapes[b]
            predicted[idx] = float(model[b + "Amp"]) * np.interp(phase[idx] + shift, xp, yp, period=1)
        return np.bincount(exposure_group, weights=predicted) / np.bincount(exposure_group)

    def score(pred, training, testing):
        err = []
        for b in eligible:
            train = training & (bands == b)
            test = testing & (bands == b)
            if train.any() and test.any():
                err.extend(y[test] - pred[test] - (y[train] - pred[train]).mean())
        if not err:
            return None
        return {"n": len(err), "rmse_mag": round(float(np.sqrt(np.mean(np.square(err)))), 4)}

    bands_out = {}
    for b in np.unique(bands):
        indices = np.flatnonzero(bands == b)
        a = indices[y[indices].argmin()]
        z = indices[y[indices].argmax()]
        span = y[z] - y[a]
        amp = float(model[b + "Amp"])
        excess = span - amp
        bands_out[b] = {
            "nights": int(indices.size),
            "published_amplitude_mag": amp,
            "brightest": [int(dates[a]), round(y[a], 4), round(e[a], 4)],
            "faintest": [int(dates[z]), round(y[z], 4), round(e[z], 4)],
            "span_mag": round(span, 4),
            "span_excess_over_published_mag": round(excess, 4),
            "span_excess_over_measured_error": round(excess / math.hypot(e[z], e[a]), 2),
        }
    result = {
        "exposures": int(len(flat)),
        "nights": int(obs_days.size),
        "retained_filter_dates": n,
        "phase_bins_of_10": int(len(np.unique(np.floor(phase * 10)))),
        "phase_correction_range_minutes": [
            round(float(correction.min()) * 1440, 2),
            round(float(correction.max()) * 1440, 2),
        ],
        "bands": bands_out,
    }
    if not selected.any():
        return result
    p = prediction(0)
    zero = np.zeros(n)
    mid = obs_days[len(obs_days) // 2]
    train = selected & (dates < mid)
    test = selected & (dates >= mid)
    shift_scores = np.array(
        [score(prediction(shift), selected, selected)["rmse_mag"] for shift in np.arange(100) / 100]
    )
    result.update(
        all_dates={"constant": score(zero, selected, selected), "catalog": score(p, selected, selected)},
        chronological={
            "split_mjd": int(mid),
            "constant": score(zero, train, test),
            "catalog": score(p, train, test),
        },
        phase_shift_grid={
            "n": len(shift_scores),
            "fraction_at_least_as_good_as_catalog": float((shift_scores <= shift_scores[0]).mean()),
            "best_shift_cycles": round(float(shift_scores.argmin()) / 100, 2),
            "best_rmse_mag": float(shift_scores.min()),
        },
    )
    return result


def run(fetch=False):
    model = literature(NAME, ["J/AJ/153/204/table5"], False)["J/AJ/153/204/table5"][0]
    mean = table(
        "mean",
        "mean",
        {"ra": RA, "dec": DEC, "radius": 5 / 3600},
        ["objID", "raMean", "decMean", "nDetections"],
        fetch,
    )
    local = [{**r, "separation_arcsec": round(sep(float(r["raMean"]), float(r["decMean"])), 3)} for r in mean]
    (target,) = [r for r in local if r["separation_arcsec"] <= 1]
    rows = table(
        "detection",
        "detections",
        {"objID": target["objID"]},
        [
            "objID",
            "detectID",
            "filterID",
            "obsTime",
            "ra",
            "dec",
            "psfFlux",
            "psfFluxErr",
            "psfQfPerfect",
            "infoFlag",
            "infoFlag2",
            "infoFlag3",
        ],
        fetch,
    )
    if (
        len({r["detectID"] for r in rows}) != len(rows)
        or any(r["objID"] != target["objID"] for r in rows)
        or len(rows) != int(target["nDetections"])
    ):
        raise ValueError("Incomplete PS1 epoch set")
    shapes = templates(model, fetch)
    return {
        "sources_within_5arcsec": local,
        "mjd_range_all": [min(float(r["obsTime"]) for r in rows), max(float(r["obsTime"]) for r in rows)],
        "primary_snr5_clean85": analyze(rows, shapes, model, 5, 0.85),
        "strict_snr10_clean90": analyze(rows, shapes, model, 10, 0.9),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
