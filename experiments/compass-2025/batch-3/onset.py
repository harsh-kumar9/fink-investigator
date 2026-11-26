"""Reduce the ASASSN-25ba observing windows from signed flux at a fixed position."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
import compass_aperture as aperture
import compass_direct as direct
from image_fetch import from_alternate_cutout, from_original
import pilot

POS = (309.8640412, -9.5602805)
NAME = "ZTF25aaohzdc"
FOLDER = HERE / "onset-sources"
COLUMNS = "field,ccdid,qid,filtercode,obsjd,filefracday,expid,imgtypecode,maglimit,seeing,infobits,pid"


def tables(fetch):
    params = {"POS": f"{POS[0]},{POS[1]}", "ct": "csv", "COLUMNS": COLUMNS}
    url = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?" + urlencode(
        params
    )
    rows = direct.capture(FOLDER, "inventory", url, set(COLUMNS.split(",")), fetch)
    if (
        not rows
        or len({r["pid"] for r in rows}) != len(rows)
        or any(
            abs(float(r["in_ra"]) - POS[0]) > 1e-7
            or abs(float(r["in_dec"]) - POS[1]) > 1e-7
            for r in rows
        )
    ):
        raise ValueError("Invalid image inventory")
    cols = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag"
    query = f"SELECT TOP 1000 {cols} FROM ztf_objects_dr24 WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{POS[0]},{POS[1]},0.024))=1 AND fid IN (1,2)"
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    stars = direct.capture(
        FOLDER,
        "stars",
        pilot.IRSA + "?" + urlencode(params),
        set(cols.split(",")),
        fetch,
    )
    if not 0 < len(stars) < 1000 or len({s["oid"] for s in stars}) != len(stars):
        raise ValueError("Invalid field-star cone")
    return rows, stars


def choose(rows):
    by_date = defaultdict(list)
    for r in rows:
        if (
            r["imgtypecode"] == "o"
            and r["filtercode"] in ("zg", "zr")
            and int(r["infobits"]) < 33554432
            and r["seeing"]
            and float(r["seeing"]) <= 3.5
            and r["maglimit"]
            and float(r["maglimit"]) >= 19
        ):
            by_date[pilot.day_of(float(r["obsjd"]))].append(r)

    def best(rows):
        return min(
            rows,
            key=lambda r: (-float(r["maglimit"]), float(r["obsjd"]), int(r["pid"])),
        )

    sample = []
    months = defaultdict(list)
    for d, v in sorted(by_date.items()):
        if date(2025, 4, 1) <= d <= date(2025, 6, 15):
            for f in sorted({r["filtercode"] for r in v}):
                sample.append(
                    ("event", d, best([r for r in v if r["filtercode"] == f]))
                )
        elif d.year == 2026:
            sample.append(("following_season", d, best(v)))
        elif 2018 <= d.year <= 2024:
            months[d.year, d.month].append(d)
    for ds in months.values():
        d = min(ds, key=lambda d: (abs(d.day - 15), d))
        sample.append(("earlier_months", d, best(by_date[d])))
    return sorted(sample, key=lambda v: float(v[2]["obsjd"])), by_date


def image_record(item, stars, fetch):
    window, d, r = item
    url = aperture.cutout_url(r, *POS)
    folder = FOLDER / "cutouts"
    path = folder / (r["pid"] + ".fits.gz")
    failure = folder / (r["pid"] + ".unavailable.json")
    info = {
        "window": window,
        "utc_date": str(d),
        "obsjd": float(r["obsjd"]),
        "band": r["filtercode"][1],
        "pid": r["pid"],
        "expid": r["expid"],
        "seeing": float(r["seeing"]),
        "global_limit": float(r["maglimit"]),
        "url": url,
    }
    receipt = path.with_name(r["pid"] + ".original.json")
    if fetch and not path.exists():
        folder.mkdir(parents=True, exist_ok=True)
        data = None
        if not failure.exists():
            try:
                data = pilot.curl(url, [])
            except subprocess.CalledProcessError as exc:
                if exc.returncode not in (22, 28):
                    raise
                pilot.json_write(
                    failure,
                    {
                        "url": url,
                        "curl_status": exc.returncode,
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(
                            timespec="seconds"
                        ),
                    },
                )
        if data is None:
            failed_cutout = json.loads(failure.read_text())
            if failed_cutout["url"] != url:
                raise ValueError("Failed request does not match inventory")
            try:
                data, original = from_alternate_cutout(r, POS)
            except subprocess.CalledProcessError as exc:
                if exc.returncode not in (22, 28):
                    raise
                try:
                    data, original = from_original(r, POS)
                except subprocess.CalledProcessError as original_error:
                    if original_error.returncode not in (22, 28):
                        raise
                    failed_cutout.update(
                        alternate_curl_status=exc.returncode,
                        original_curl_status=original_error.returncode,
                        original_http_status=int(original_error.stdout or 0),
                        original_url=url.split("?")[0],
                        last_attempt_utc=datetime.now(timezone.utc).isoformat(
                            timespec="seconds"
                        ),
                    )
                    pilot.json_write(failure, failed_cutout)
                    return {**info, "status": "cutout unavailable"}
            pilot.json_write(receipt, {**original, "ibe_cutout_failure": failed_cutout})
        if not data.startswith(b"SIMPLE  ="):
            raise ValueError("Invalid FITS cutout: " + url)
        temp = path.with_suffix(".tmp")
        temp.write_bytes(gzip.compress(data, compresslevel=5, mtime=0))
        temp.replace(path)
        failure.unlink(missing_ok=True)
    if not path.exists() and json.loads(failure.read_text())["url"] == url:
        return {**info, "status": "cutout unavailable"}
    raw = gzip.decompress(path.read_bytes())
    digest = hashlib.sha256(raw).hexdigest()
    if receipt.exists():
        original = json.loads(receipt.read_text())
        if (
            original["sha256_cutout"] != digest
            or original["original_url"] != url.split("?")[0]
            or original["method"] not in ("local_crop", "alternate_ibe_cutout")
        ):
            raise ValueError("Original image receipt differs from cutout")
        info["source"] = (
            "local crop of original quadrant"
            if original["method"] == "local_crop"
            else "alternate IRSA cutout"
        )
    else:
        info["source"] = "IRSA cutout"
    summary = aperture.measure(raw, r, POS, stars)
    return {
        **info,
        "sha256": digest,
        **summary,
    }


def run(fetch):
    rows, stars = tables(fetch)
    sample, dates = choose(rows)
    with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
        measures = list(pool.map(lambda item: image_record(item, stars, fetch), sample))
    summary = {}
    for window in ("event", "earlier_months", "following_season"):
        v = [r for r in measures if r["window"] == window]
        tested = [
            r for r in v if r["status"] == "measured" and r["five_sigma_limit"] >= 18
        ]
        summary[window] = {
            "selected_exposures": len(v),
            "selected_dates": len({r["utc_date"] for r in v}),
            "measured": sum(r["status"] == "measured" for r in v),
            "deep_enough": len(tested),
            "bright": [
                r["pid"] for r in tested if r["mag"] is not None and r["mag"] <= 18
            ],
            "first": min((r["utc_date"] for r in tested), default=None),
            "last": max((r["utc_date"] for r in tested), default=None),
        }
        control = [r["held_out_star"] for r in v if r.get("held_out_star")]
        summary[window]["control"] = {
            "tested": len(control),
            "within_0.3_and_snr5": sum(
                bool(c["snr"] >= 5 and abs(c["measured_mag"] - c["catalog_mag"]) <= 0.3)
                for c in control
            ),
        }
        if window == "following_season" and tested:
            days = sorted(r["obsjd"] for r in tested)
            ranges = {
                "first_to_last_observation": (days[0], days[-1]),
                "calendar_window": (
                    pilot.utc_jd(date(2026, 1, 1)),
                    pilot.utc_jd(date(2026, 11, 26)),
                ),
            }
            summary[window]["start_time_coverage"] = {
                key: {
                    str(length): coverage(days, a, b, length)
                    for length in (1, 3, 7, 14)
                }
                for key, (a, b) in ranges.items()
            }
    return {
        "inventory_exposures": len(rows),
        "qualified_dates_per_year": dict(
            sorted(Counter(str(d.year) for d in dates).items())
        ),
        "summary": summary,
        "photometry": measures,
    }


def coverage(times, start, end, duration):
    """Length of the union of plateau starts for which an exposure sees the source."""
    covered = 0.0
    right = start
    for t in sorted(times):
        a, b = max(start, t - duration), min(end, t)
        if b > max(right, a):
            covered += b - max(right, a)
            right = b
    return covered / (end - start)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
