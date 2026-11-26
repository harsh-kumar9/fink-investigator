"""Measure the unresolved host and transient on original images outside the alert stream."""

import argparse
from collections import defaultdict
from datetime import date
import io
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, median, stdev
import subprocess
from urllib.parse import urlencode
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE.parent / "batch-3"))
import compass_aperture as phot
from compass_direct import capture
import image_fetch
import pilot

POS = (250.0278614, 29.763246)
FOLDER = HERE / "late-sources"
COLS = "field,ccdid,qid,filtercode,obsjd,filefracday,expid,imgtypecode,maglimit,seeing,infobits,pid"
WINDOWS = {
    "baseline": ("2023-01-01", "2024-12-31", 24),
    "peak": ("2025-04-20", "2025-05-02", 8),
    "june": ("2025-06-01", "2025-06-30", 8),
    "july": ("2025-07-01", "2025-07-31", 8),
    "august": ("2025-08-01", "2025-08-31", 8),
    "september": ("2025-09-01", "2025-09-30", 8),
    "control_2026": ("2026-01-01", "2026-11-26", 12),
}


def tables(fetch):
    url = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + urlencode(
        dict(POS=f"{POS[0]},{POS[1]}", ct="csv", COLUMNS=COLS)
    )
    inventory = capture(FOLDER, "inventory", url, set(COLS.split(",")), fetch)
    if len({r["pid"] for r in inventory}) != len(inventory) or not inventory:
        raise ValueError("Duplicate or empty inventory")
    cols = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag"
    query = (
        f"SELECT TOP 2000 {cols} FROM ztf_objects_dr24 WHERE CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{POS[0]},{POS[1]},0.030))=1 AND fid=2"
    )
    stars = capture(
        FOLDER,
        "stars",
        pilot.IRSA
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)),
        set(cols.split(",")),
        fetch,
    )
    if len(stars) >= 2000 or len({r["oid"] for r in stars}) != len(stars):
        raise ValueError("Invalid stellar cone")
    return inventory, stars


def select(inventory, field=None):
    cells = defaultdict(lambda: defaultdict(dict))
    for r in inventory:
        if (
            r["imgtypecode"] != "o"
            or r["filtercode"] != "zr"
            or int(r["infobits"]) >= 33554432
            or not r["seeing"]
            or not r["maglimit"]
            or float(r["seeing"]) > 3.2
            or float(r["maglimit"]) < 20
        ):
            continue
        day = str(pilot.day_of(float(r["obsjd"])))
        key = tuple(int(r[k]) for k in ("field", "ccdid", "qid"))
        for w, (a, b, _) in WINDOWS.items():
            if a <= day <= b:
                days = cells[key][w]
                order = lambda m: (
                    -float(m["maglimit"]),
                    float(m["obsjd"]),
                    int(m["pid"]),
                )
                if day not in days or order(r) < order(days[day]):
                    days[day] = r

    def rank(key):
        v = cells[key]
        return (
            -sum(len(v[k]) for k in ("june", "july", "august", "september")),
            -len(v["baseline"]),
            key,
        )

    key = min(
        (
            k
            for k, v in cells.items()
            if len(v["baseline"]) >= 6 and len(v["control_2026"]) >= 4
        ),
        key=rank,
    )
    if field is not None:
        candidates = [
            k for k in cells if k[0] == field and len(cells[k]["baseline"]) >= 6
        ]
        key = min(candidates, key=rank)
    out = []
    for w, (_, _, n) in WINDOWS.items():
        ds = sorted(cells[key][w])
        n = min(n, len(ds))
        inds = (
            [round(i * (len(ds) - 1) / (n - 1)) for i in range(n)]
            if n > 1
            else list(range(n))
        )
        out.extend((w, ds[i], cells[key][w][ds[i]]) for i in inds)
    return (
        key,
        {str(k): {w: len(ds) for w, ds in v.items()} for k, v in cells.items()},
        out,
    )


def image(row, fetch):
    folder = FOLDER / "cutouts"
    pid = row["pid"]
    path, rec = folder / (pid + ".fits.gz"), folder / (pid + ".json")
    url = phot.cutout_url(row, *POS, size=225)
    if fetch and not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        failed = []
        for mode in ("service", "alternate", "original"):
            try:
                if mode == "service":
                    raw = pilot.curl(url, [])
                    receipt = {"method": mode, "url": url}
                elif mode == "alternate":
                    raw, receipt = image_fetch.from_alternate_cutout(row, POS, size=225)
                else:
                    raw, receipt = image_fetch.from_original(row, POS, size=225)
            except subprocess.CalledProcessError as exc:
                if exc.returncode not in (22, 28):
                    raise
                failed.append({"method": mode, "curl_status": exc.returncode})
                continue
            if not raw.startswith(b"SIMPLE  ="):
                raise ValueError("Invalid image")
            with fits.open(io.BytesIO(raw)) as hdus:
                if hdus[0].header["EXPID"] != int(row["expid"]):
                    raise ValueError("Wrong exposure")
            dest = path.with_suffix(".tmp")
            dest.write_bytes(gzip.compress(raw, compresslevel=6, mtime=0))
            pilot.json_write(
                rec,
                {
                    **receipt,
                    "requested": url,
                    "sha256": hashlib.sha256(raw).hexdigest(),
                    "attempts_before_success": failed,
                    "archived_utc": datetime.now(timezone.utc).isoformat(
                        timespec="seconds"
                    ),
                },
            )
            dest.replace(path)
            break
        else:
            pilot.json_write(
                rec,
                {
                    "method": "unavailable",
                    "requested": url,
                    "failed": failed,
                    "archived_utc": datetime.now(timezone.utc).isoformat(
                        timespec="seconds"
                    ),
                },
            )
    meta = json.loads(rec.read_text())
    if meta["requested"] != url:
        raise ValueError("Changed image request")
    if meta["method"] == "unavailable":
        if path.exists():
            raise ValueError("Contradictory receipt")
        return None
    raw = gzip.decompress(path.read_bytes())
    if hashlib.sha256(raw).hexdigest() != meta["sha256"]:
        raise ValueError("Changed cutout")
    return raw


def sky_controls(stars, ref):
    _, sky, _ = phot.star_list(stars, ref, *POS, radius=95)
    trial = [
        (x, y)
        for y in (-45, -25, 0, 25, 45)
        for x in (-45, -25, 0, 25, 45)
        if (x, y) != (0, 0)
    ]
    result = {}
    for x, y in trial:
        ra = POS[0] + x / 3600 / math.cos(math.radians(POS[1]))
        dec = POS[1] + y / 3600
        if all(math.hypot(*pilot.offsets(ra, dec, *point)) >= 12 for point in sky):
            result[f"{x},{y}"] = (ra, dec)
    return result


def measure(item, stars, controls, fetch):
    win, d, r = item
    out = dict(
        window=win,
        date=d,
        pid=r["pid"],
        mjd=float(r["obsjd"]) - 2400000.5,
        seeing=float(r["seeing"]),
        field=int(r["field"]),
        limit=float(r["maglimit"]),
    )
    raw = image(r, fetch)
    if raw is None:
        return {**out, "status": "unavailable"}
    v = phot.measure(raw, r, POS, stars, standard_radius=95)
    out.update(status=v["status"], standards=v["standards"])
    if v["status"] != "measured":
        return out
    factor = 10 ** (0.4 * (23.9 - v["zero_point"]))
    out.update(
        flux_ujy=v["signal"] * factor,
        error_ujy=math.hypot(
            v["sigma"] * factor,
            v["signal"] * factor * math.log(10) * 0.4 * v["zero_point_scatter"],
        ),
        control=v["held_out_star"],
        blank_flux_ujy={},
        five_sigma_limit=v["five_sigma_limit"],
    )
    with fits.open(io.BytesIO(raw)) as hdus:
        wcs = WCS(hdus[0].header)
        for key, pos in controls.items():
            x, y = wcs.all_world2pix([pos], 0)[0]
            counts = phot.aperture(hdus[0].data.astype(float, copy=False), x, y)
            if counts is not None:
                out["blank_flux_ujy"][key] = (counts - v["blank_median"]) * factor
    return out


def summarize(v):
    if len(v) < 4:
        return None
    flux = [r["flux_ujy"] for r in v]
    return dict(
        n=len(v),
        mean_flux=mean(flux),
        error=max(
            stdev(flux) / math.sqrt(len(v)),
            math.sqrt(sum(r["error_ujy"] ** 2 for r in v)) / len(v),
        ),
        mean_seeing=mean(r["seeing"] for r in v),
    )


def contrast(base, v):
    a, b = summarize(base), summarize(v)
    if not a or not b:
        return {"baseline": a, "window": b}
    delta = b["mean_flux"] - a["mean_flux"]
    err = math.hypot(a["error"], b["error"])
    blanks = {}
    keys = set.intersection(*(set(r["blank_flux_ujy"]) for r in base + v))
    for key in sorted(keys):
        x = [r["blank_flux_ujy"][key] for r in base]
        y = [r["blank_flux_ujy"][key] for r in v]
        diff = mean(y) - mean(x)
        error = math.sqrt(stdev(x) ** 2 / len(x) + stdev(y) ** 2 / len(y))
        blanks[key] = {"delta_ujy": diff, "z": diff / error if error > 0 else None}

    def accepted_controls(rows):
        by_star = defaultdict(list)
        for row in rows:
            ctrl = row["control"]
            if ctrl and ctrl["snr"] >= 5:
                by_star[ctrl["oid"]].append(ctrl["measured_mag"] - ctrl["catalog_mag"])
        return by_star

    pre, post = accepted_controls(base), accepted_controls(v)
    checks = {
        oid: {
            "baseline_dates": len(pre[oid]),
            "window_dates": len(post[oid]),
            "median_change_mag": median(post[oid]) - median(pre[oid]),
        }
        for oid in sorted(pre.keys() & post.keys())
    }
    return {
        "held_out_standards": checks,
        "baseline": a,
        "window": b,
        "delta_ujy": delta,
        "error_ujy": err,
        "z": delta / err,
        "blank_controls": blanks,
    }


def run(fetch=False, field=None):
    inventory, stars = tables(fetch)
    detector, coverage, sample = select(inventory, field)
    controls = sky_controls(stars, sample[0][2])
    with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
        measurements = list(
            pool.map(lambda x: measure(x, stars, controls, fetch), sample)
        )
    valid = [r for r in measurements if r["status"] == "measured"]
    base = [r for r in valid if r["window"] == "baseline"]
    result = {}
    for w in WINDOWS:
        if w == "baseline":
            continue
        v = [r for r in valid if r["window"] == w]
        common = lambda a, b: [
            r for r in a if any(abs(r["seeing"] - s["seeing"]) <= 0.4 for s in b)
        ]
        result[w] = {
            "all": contrast(base, v),
            "matched_seeing": contrast(common(base, v), common(v, base)),
        }
    return {
        "detector": detector,
        "inventory_exposures": len(inventory),
        "qualifying_dates": coverage,
        "sky_controls": controls,
        "measurements": measurements,
        "contrasts": result,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--field", type=int, choices=[633])
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    out = args.out or HERE / (
        "results/late-other.json" if args.field else "results/late.json"
    )
    pilot.json_write(out, run(args.fetch, args.field))
