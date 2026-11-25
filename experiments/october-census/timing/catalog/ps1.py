"""Check the two catalog frequencies in resolved archival Pan-STARRS photometry."""

from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from statistics import median

import numpy as np

import check

HERE = check.HERE


def read(key, receipt):
    raw = gzip.decompress((HERE / (key + ".csv.gz")).read_bytes())
    if hashlib.sha256(raw).hexdigest() != receipt[key]["sha256_uncompressed"]:
        raise ValueError(f"Changed Pan-STARRS response: {key}")
    return list(csv.DictReader(io.StringIO(raw.decode())))


def fit(rows, period_days):
    days = defaultdict(list)
    for r in rows:
        t, flux, err, quality = (float(r[k]) for k in ("obsTime", "psfFlux", "psfFluxErr", "psfQfPerfect"))
        if all(map(math.isfinite, (t, flux, err, quality))) and err > 0 and flux/err >= 10 and quality >= .9:
            days[(r["filterID"], math.floor(t))].append((t, -2.5*math.log10(flux)))
    bands = sorted({f for f,_ in days if sum(key[0] == f for key in days) >= 5})
    points = [(f, median(t for t,y in rows), median(y for t,y in rows)) for (f,_),rows in days.items() if f in bands]
    result = {"raw_detections":len(rows), "retained_filter_nights":len(days),
              "used_filter_nights":len(points), "used_filters": bands}
    if len(points) < 15:
        return result
    t = np.array([t for f,t,y in points])
    result["mjd_range"] = [float(t.min()), float(t.max())]
    y = np.array([y for f,t,y in points])
    constant = np.array([[float(f == b) for b in bands] for f,t,y in points])
    def error(x):
        return np.sum((y - x @ np.linalg.lstsq(x,y,rcond=None)[0])**2)
    result["fraction_scatter_explained"] = {}
    for key, period in period_days.items():
        angle = 2*np.pi*(t-t.min())/period
        x = np.column_stack((constant, np.cos(angle), np.sin(angle)))
        result["fraction_scatter_explained"][key] = round(float(1 - error(x)/error(constant)), 4)
    return result


def run():
    receipt = json.loads((HERE / "ps1-sources.json").read_text())
    stars = read("ps1-mean", receipt)
    catalog = check.run()
    timing = json.loads((HERE.parent / "results.json").read_text())
    target, = [r for r in timing["period_search"] if r["id"] == "ZTF18abjtfpx"]
    periods = {"ztf":target["period_days"], "atlas_lombscar":float(catalog["atlas"]["ls-Pday"])}
    result = []
    for star in stars:
        sep = check.separation(float(star["raMean"]), float(star["decMean"]))
        if sep > 4:
            continue
        key = star["objID"]
        rows = read(key, receipt)
        if len({r["detectID"] for r in rows}) != int(star["nDetections"]) or len(rows) != int(star["nDetections"]):
            raise ValueError(f"Incomplete or repeated detections for {key}")
        result.append({"objID":key, "offset_arcsec":round(sep,3), **fit(rows, periods)})
    return {"periods_days":periods, "sources":sorted(result, key=lambda r:r["offset_arcsec"])}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
