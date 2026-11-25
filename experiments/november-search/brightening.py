"""Test the r-selected state change on separate g-band exposures and field stars."""

import argparse
from collections import defaultdict
import csv
from datetime import date, datetime, timedelta, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import median, quantiles
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import census_photometry
import pilot

ALERT_RA, ALERT_DEC = 92.6614871, 20.6345684
EPOCH = date(1858, 11, 17)
START = date(2021, 1, 1)
SPLIT = date(2023, 1, 1)
END = date(2026, 1, 1)


def capture(name, url, required, fetch):
    folder = HERE / "sources/brightening"
    meta = folder / (name + ".json")
    path = folder / (name + ".csv.gz")
    fresh = fetch and not meta.exists()
    if fresh:
        folder.mkdir(parents=True, exist_ok=True)
        raw = pilot.curl(url, [])
    else:
        raw = gzip.decompress(path.read_bytes())
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"Missing columns for {name}: {raw[:200]!r}")
    rows = list(reader)
    digest = hashlib.sha256(raw).hexdigest()
    if fresh:
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            meta,
            {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": digest,
            },
        )
    receipt = json.loads(meta.read_text())
    if receipt["sha256_uncompressed"] != digest or receipt["url"] != url:
        raise ValueError(f"{name} request or checksum mismatch")
    return rows


def distance(row, ra, dec):
    return math.hypot(*pilot.offsets(float(row["ra"]), float(row["dec"]), ra, dec))


def sources(fetch):
    q = (
        "SELECT TOP 2000 oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag FROM ztf_objects_dr24 "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ALERT_RA},{ALERT_DEC},{120/3600}))=1 "
        "AND field=563 AND ccdid=12 AND qid=1 AND fid=1 AND nobs>=100"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
    )
    rows = capture(
        "cone",
        url,
        {"oid", "ra", "dec", "meanmag", "nobs", "fid", "field", "ccdid", "qid"},
        fetch,
    )
    if not 0 < len(rows) < 2000 or any(
        (r["field"], r["ccdid"], r["qid"], r["fid"]) != ("563", "12", "1", "1")
        or int(r["nobs"]) < 100
        for r in rows
    ):
        raise ValueError("Unexpected cone coverage or detector")
    for r in rows:
        r["distance_arcsec"] = distance(r, ALERT_RA, ALERT_DEC)
    target = min(rows, key=lambda r: (r["distance_arcsec"], r["oid"]))
    if target["distance_arcsec"] > 1.5:
        raise ValueError("No green source within 1.5 arcsec")
    ordered = sorted(
        rows,
        key=lambda r: (
            abs(float(r["meanmag"]) - float(target["meanmag"])),
            r["distance_arcsec"],
            r["oid"],
        ),
    )
    controls = []
    for r in ordered:
        if (
            10 <= r["distance_arcsec"] <= 120
            and abs(float(r["meanmag"]) - float(target["meanmag"])) <= 1
            and all(distance(r, float(c["ra"]), float(c["dec"])) >= 5 for c in controls)
        ):
            controls.append(r)
        if len(controls) == 5:
            break
    if len(controls) < 3:
        raise ValueError("Insufficient comparison stars")
    return target, controls


def shift(series):
    early = [v for d, v in series.items() if d < SPLIT]
    late = [v for d, v in series.items() if d >= SPLIT]
    if not early or not late:
        return None
    return round(median(early) - median(late), 4)


def nightly(by_day):
    return {d: median(m) for d, m in by_day.items()}


def red_reference():
    base = HERE / "sources/direct/ZTF18abqclig"
    raw = gzip.decompress(base.with_suffix(".csv.gz").read_bytes())
    if (
        hashlib.sha256(raw).hexdigest()
        != json.loads(base.with_suffix(".json").read_text())["sha256_uncompressed"]
    ):
        raise ValueError("Red-band source checksum mismatch")
    by_day = defaultdict(list)
    for row in csv.DictReader(io.StringIO(raw.decode())):
        day = EPOCH + timedelta(days=int(float(row["mjd"])))
        if (
            row["oid"] == "563212100016993"
            and row["filtercode"] == "zr"
            and int(row["catflags"]) == 0
            and 0 < float(row["magerr"]) <= 0.2
            and START <= day < END
        ):
            by_day[day].append(float(row["mag"]))
    return nightly(by_day)


def annual_medians(series):
    summary = {}
    for year in range(2021, 2026):
        mags = [v for d, v in series.items() if d.year == year]
        summary[year] = {
            "dates": len(mags),
            "mag": round(median(mags), 4) if mags else None,
            "interquartile_mag": (
                [round(v, 4) for v in quantiles(mags, n=4, method="inclusive")[::2]]
                if len(mags) > 1
                else None
            ),
        }
    return summary


def run(fetch):
    target, controls = sources(fetch)
    ids = [r["oid"] for r in [target] + controls]
    options = [("ID", i) for i in ids]
    options += [
        ("COLLECTION", "ztf_dr24"),
        ("TIME", f"{(START-EPOCH).days} {(END-EPOCH).days}"),
        ("BAD_CATFLAGS_MASK", "65535"),
        ("FORMAT", "csv"),
    ]
    url = census_photometry.SERVICE + "?" + urlencode(options)
    rows = capture(
        "green",
        url,
        {"oid", "expid", "mag", "magerr", "mjd", "filtercode", "catflags"},
        fetch,
    )
    if {r["oid"] for r in rows} != set(ids) or any(
        r["filtercode"] != "zg" for r in rows
    ):
        raise ValueError("Unexpected green series")
    curves = {i: {} for i in ids}
    for r in rows:
        day = EPOCH + timedelta(days=int(float(r["mjd"])))
        if (
            START <= day < END
            and int(r["catflags"]) == 0
            and 0 < float(r["magerr"]) <= 0.2
        ):
            point = (day, float(r["mag"]))
            if curves[r["oid"]].setdefault(r["expid"], point) != point:
                raise ValueError("Conflicting rows for one exposure")
    comparison = {i: median(m for _, m in curves[i].values()) for i in ids[1:]}
    # Use the comparison stars to count available visits even when the target
    # itself was too faint for the catalog's photometry selection.
    reference_visits = defaultdict(list)
    for i in comparison:
        for exp, (day, _) in curves[i].items():
            reference_visits[exp].append(day)
    visit_dates = {days[0] for days in reference_visits.values() if len(days) >= 3}
    target_days, corrected_days = defaultdict(list), defaultdict(list)
    paired = {i: defaultdict(list) for i in comparison}
    control_days = {i: defaultdict(list) for i in comparison}
    for exp, (day, mag) in curves[ids[0]].items():
        shared = {
            i: curves[i][exp][1]
            for i in comparison
            if exp in curves[i] and curves[i][exp][0] == day
        }
        if len(shared) < 3:
            continue
        target_days[day].append(mag)
        corrected_days[day].append(
            mag - median(m - comparison[i] for i, m in shared.items())
        )
        for i, m in shared.items():
            paired[i][day].append(mag - m)
            control_days[i][day].append(m)
    corrected = nightly(corrected_days)
    per_month = {
        m: shift({d: v for d, v in corrected.items() if d.month == m})
        for m in range(1, 13)
    }
    available = [s for s in per_month.values() if s is not None]
    cs = {
        i: {
            "paired_target_shift_mag": shift(nightly(paired[i])),
            "control_shift_mag": shift(nightly(control_days[i])),
        }
        for i in comparison
    }
    early = sum(d < SPLIT for d in corrected)
    late = len(corrected) - early
    # Adversarial values bound the influence of missing whole nights. The
    # smallest early-minus-late shift assigns bright early and faint late dates.
    early_min = median(corrected.get(d, -math.inf) for d in visit_dates if d < SPLIT)
    late_max = median(corrected.get(d, math.inf) for d in visit_dates if d >= SPLIT)
    missing_bound = early_min - late_max
    delta = shift(corrected)
    passed = (
        min(early, late) >= 20
        and delta is not None
        and delta >= 0.5
        and sum(
            v["paired_target_shift_mag"] is not None
            and v["paired_target_shift_mag"] > 0
            for v in cs.values()
        )
        >= 3
    )
    red = red_reference()
    return {
        "object": "ZTF18abqclig",
        "target": target,
        "comparison_stars": controls,
        "green_response_rows": len(rows),
        "early_dates": early,
        "late_dates": late,
        "comparison_field_dates_early": sum(d < SPLIT for d in visit_dates),
        "comparison_field_dates_late": sum(d >= SPLIT for d in visit_dates),
        "worst_case_missing_night_shift_lower_bound_mag": (
            round(missing_bound, 4) if math.isfinite(missing_bound) else None
        ),
        "raw_target_shift_mag": shift(nightly(target_days)),
        "corrected_shift_mag": delta,
        "matched_month_shift_mag": round(median(available), 4) if available else None,
        "individual_comparisons": cs,
        "monthly_shifts_mag": per_month,
        "annual_corrected_medians": annual_medians(corrected),
        "exploratory_red": {
            "shift_mag": shift(red),
            "annual_medians": annual_medians(red),
        },
        "criterion_passed": passed,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    result = run(parser.parse_args().fetch)
    pilot.json_write(HERE / "results/brightening.json", result)
    print(json.dumps(result, indent=2))
