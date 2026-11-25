"""Replay the green-band prediction using the registered October red ephemeris."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
from urllib.parse import urlencode

from astropy.timeseries import LombScargle
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import analyze

RA, DEC = 280.305911, 9.7241705
NAME = "ZTF18abjtfpx"
STOP_MJD = 2461333.6713194 - 2400000.5


def archive(name, url, columns, fetch):
    file = HERE / (name + ".csv.gz")
    receipt = HERE / (name + ".json")
    if fetch and not receipt.exists():
        data = subprocess.check_output(["curl", "-fsSL", "--max-time", "180", "--retry", "1", url])
        reader = csv.DictReader(io.StringIO(data.decode()))
        if not set(columns) <= set(reader.fieldnames or []):
            raise ValueError(f"Unexpected archive response for {name}: {data[:160]!r}")
        file.write_bytes(gzip.compress(data, mtime=0))
        receipt.write_text(json.dumps({"url": url, "sha256": hashlib.sha256(data).hexdigest(),
                                      "retrieved_utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n")
    data = gzip.decompress(file.read_bytes())
    meta = json.loads(receipt.read_text())
    if url != meta["url"] or hashlib.sha256(data).hexdigest() != meta["sha256"]:
        raise ValueError(f"Archive mismatch for {name}")
    return list(csv.DictReader(io.StringIO(data.decode())))


def separation(a, b):
    dlon, dlat = np.radians(np.array(a) - b)
    h = math.sin(dlat/2)**2 + math.cos(math.radians(a[1])) * math.cos(math.radians(b[1])) * math.sin(dlon/2)**2
    return math.degrees(2 * math.asin(math.sqrt(h))) * 3600


def position(row):
    return float(row["ra"]), float(row["dec"])


def select(fetch):
    columns = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag".split(",")
    query = (f"SELECT TOP 10000 {','.join(columns)} FROM ztf_objects_dr24 "
             f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},{120/3600}))=1 "
             "AND fid=1 AND nobs>=30")
    url = "https://irsa.ipac.caltech.edu/TAP/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "CSV", "QUERY": query})
    rows = archive("field", url, columns, fetch)
    if len(rows) >= 10000 or len({r["oid"] for r in rows}) != len(rows):
        raise ValueError("Incomplete or duplicated IRSA field")
    for r in rows:
        r["separation_arcsec"] = separation(position(r), (RA, DEC))
    nearby = sorted((r for r in rows if r["separation_arcsec"] <= 1.5),
                    key=lambda r: (-int(r["nobs"]), r["separation_arcsec"], int(r["oid"])))
    target = nearby[0]
    chosen = [target]
    for r in sorted(rows, key=lambda r: (r["separation_arcsec"], int(r["oid"]))):
        if (15 <= r["separation_arcsec"] <= 120
            and all(r[k] == target[k] for k in ("field", "ccdid", "qid"))
            and r["meanmag"] and abs(float(r["meanmag"]) - float(target["meanmag"])) <= 0.5
            and all(separation(position(old), position(r)) > 3 for old in chosen)):
            chosen.append(r)
            if len(chosen) == 3:
                break
    return nearby, chosen


def green_data(stars, fetch):
    ids = [r["oid"] for r in stars]
    params = [("ID", oid) for oid in ids]
    params += [("COLLECTION", "ztf_dr24"), ("BAD_CATFLAGS_MASK", "65535"),
               ("TIME", f"-Inf {STOP_MJD:.8f}"), ("FORMAT", "csv")]
    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(params)
    rows = archive("green", url, ["oid", "filtercode", "mjd", "hjd", "expid", "mag", "magerr", "catflags"], fetch)
    if set(ids) != {r["oid"] for r in rows}:
        raise ValueError("IRSA response did not cover the exact selected IDs")
    series = {oid: {} for oid in ids}
    for row in rows:
        mjd, hjd, mag, error = (float(row[k]) for k in ("mjd", "hjd", "mag", "magerr"))
        if (row["filtercode"] != "zg" or int(row["catflags"]) or not 0 < error <= 0.2
            or not all(map(math.isfinite, [mjd, hjd, mag])) or mjd > STOP_MJD + 1e-7):
            continue
        day = math.floor(mjd)
        candidate = (error, int(row["expid"]), mjd, hjd, mag)
        best = series[row["oid"]]
        if day not in best or candidate < best[day]:
            best[day] = candidate
    return {oid: np.array([best[k][2:] for k in sorted(best)]).reshape(-1, 3) for oid, best in series.items()}


def validate(data, frequency, model, origin):
    n = len(data)*2//3
    counts = {"training_dates": n, "validation_dates": len(data)-n}
    if len(data) < 3:
        return {**counts, "status": "insufficient dates"}
    m = data[:, 2]
    predicted = model.model(data[:, 1] - origin, frequency)
    offset = np.mean(m[:n] - predicted[:n])
    gain = 1 - np.sum((m[n:] - predicted[n:] - offset)**2) / np.sum((m[n:] - np.mean(m[:n]))**2)
    rho = spearmanr(m[n:], predicted[n:]).statistic
    return {**counts, "test_first_mjd": data[n, 0], "test_last_mjd": data[-1, 0],
            "training_green_minus_red_mag": float(offset), "test_error_reduction": float(gain),
            "test_spearman": float(rho), "all_spearman": float(spearmanr(m, predicted).statistic),
            "status": "insufficient dates" if len(data)-n < 30 else
            ("passes" if gain >= 0.30 and rho >= 0.60 else "below threshold")}


def run(fetch=False):
    nearby, stars = select(fetch)
    series = green_data(stars, fetch)
    direct = json.loads((analyze.CENSUS / "results/direct.json").read_text())
    red = analyze.nights(next(c for c in direct["exploration"] if c["id"] == NAME))
    timing = json.loads((HERE.parent / "results.json").read_text())
    result = next(c for c in timing["period_search"] if c["id"] == NAME)
    n = result["training_nights"]
    origin = red[0, 1]
    freq = 1/result["period_days"]
    model = LombScargle(red[:n, 1]-origin, red[:n, 2], nterms=2)
    red_days = np.floor(red[:, 0])
    rows = []
    for i, star in enumerate(stars):
        data = series[star["oid"]]
        unique = data[~np.isin(np.floor(data[:, 0]), red_days)]
        rows.append({"role": "target" if i == 0 else "comparison", "series": star,
                     "green_only": validate(unique, freq, model, origin),
                     "all_green": validate(data, freq, model, origin)})
    return {"object_id": NAME, "red_period_days": result["period_days"], "red_origin_hjd": origin,
            "red_training_mean_mag": float(model.offset()),
            "red_model_parameters": model.model_parameters(freq).tolist(),
            "eligible_target_oids": [r["oid"] for r in nearby], "results": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2, allow_nan=False))
