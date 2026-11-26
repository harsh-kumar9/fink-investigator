"""Check the compass's three alert patterns with archived ZTF DR24 direct photometry."""

import argparse
from collections import defaultdict
import csv
from datetime import date, datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

import compass
import compass_archive
import pilot
from census_photometry import SERVICE


def capture(folder, name, url, required, fetch):
    receipt = folder / (name + ".json")
    archive = folder / (name + ".csv.gz")
    new = fetch and not receipt.exists()
    raw = pilot.curl(url, []) if new else gzip.decompress(archive.read_bytes())
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig").lstrip("\r\n")))
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"Invalid CSV response for {url}")
    rows = list(reader)
    sha = hashlib.sha256(raw).hexdigest()
    if new:
        folder.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256_uncompressed": sha,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    meta = json.loads(receipt.read_text())
    if meta["url"] != url or meta["sha256_uncompressed"] != sha:
        raise ValueError(f"Source request or checksum changed in {receipt}")
    return rows


def direct_data(case, folder, fetch=False):
    ra, dec = case["ra"], case["dec"]
    query = (
        "SELECT TOP 100 oid,ra,dec,field,ccdid,qid,fid,nobs FROM ztf_objects_dr24 "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{1.5/3600:.12f}))=1 "
        "AND fid IN (1,2)"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
        )
    )
    cone = capture(
        folder,
        "cone",
        url,
        {"oid", "ra", "dec", "fid", "nobs", "field", "ccdid", "qid"},
        fetch,
    )
    if len(cone) >= 100:
        raise ValueError("DR24 source cone reached its row cap")
    by_filter = defaultdict(list)
    for r in cone:
        offset = math.hypot(*pilot.offsets(float(r["ra"]), float(r["dec"]), ra, dec))
        if (
            offset > 1.501
            or not r["oid"].isdigit()
            or int(r["fid"]) not in compass.BANDS
        ):
            raise ValueError("Invalid direct source or position")
        by_filter[int(r["fid"])].append({**r, "separation_arcsec": round(offset, 6)})
    selected = {
        fid: min(
            rows, key=lambda r: (-int(r["nobs"]), r["separation_arcsec"], int(r["oid"]))
        )
        for fid, rows in by_filter.items()
    }
    if not selected:
        return selected, [], len(cone)
    params = [("ID", r["oid"]) for _, r in sorted(selected.items())]
    params += [
        ("COLLECTION", "ztf_dr24"),
        ("FORMAT", "csv"),
        ("BAD_CATFLAGS_MASK", "65535"),
    ]
    rows = capture(
        folder,
        "photometry",
        SERVICE + "?" + urlencode(params),
        {"oid", "expid", "mjd", "filtercode", "mag", "magerr", "catflags"},
        fetch,
    )
    expected = {v["oid"]: "z" + compass.BANDS[k] for k, v in selected.items()}
    if any(
        r["oid"] not in expected or r["filtercode"] != expected[r["oid"]] for r in rows
    ):
        raise ValueError("Unexpected object ID or filter in direct series")
    return selected, rows, len(cone)


def nightly(rows, band, anchor_jd):
    grouped = defaultdict(list)
    seen = {}
    after = set()
    for r in rows:
        if r["filtercode"] != "z" + band or int(r["catflags"]) != 0:
            continue
        t, mag, err = (float(r[k]) for k in ("mjd", "mag", "magerr"))
        if not (0 < err <= 0.2 and math.isfinite(t) and math.isfinite(mag)):
            continue
        expid = r["expid"]
        if expid in seen:
            if seen[expid] != (t, mag):
                raise ValueError("Conflicting direct exposure rows")
            continue
        seen[expid] = (t, mag)
        day = pilot.day_of(t + 2400000.5)
        if t + 2400000.5 <= anchor_jd + 1e-7:
            grouped[day].append(mag)
        else:
            after.add(day)
    return (
        {d: median(v) for d, v in grouped.items()},
        len(after),
        max([*grouped, *after], default=None),
    )


def group(data, days):
    values = [data[d] for d in sorted(days) if d in data]
    return {"dates": len(values), "median_mag": median(values) if values else None}


def gap(a, b, min_a, min_b):
    if a["dates"] >= min_a and b["dates"] >= min_b:
        return round(a["median_mag"] - b["median_mag"], 4)
    return None


def compare(case, history, rows, fid):
    band = compass.BANDS[fid]
    data, after, cutoff = nightly(rows, band, case["jd"])
    alerts = [r for r in history if r["i:jd"] <= case["jd"] and compass.quality(r)]
    signs = {
        s: {
            pilot.day_of(r["i:jd"])
            for r in alerts
            if r["i:fid"] == fid and compass.sign(r) == s
        }
        for s in [1, -1]
    }
    result = {
        "before_or_at_anchor_dates": len(data),
        "later_dates": after,
        "last_quality_date": str(cutoff) if cutoff else None,
        "off_alert": group(data, data.keys() - signs[1] - signs[-1]),
    }
    route = case["selected_for"]
    if route == "onset":
        first = date.fromisoformat(case["eligible"][route]["first_positive_date"])
        baseline = group(data, [d for d in data if 15 <= (first - d).days <= 365])
        event = group(data, signs[1])
        result.update(
            baseline=baseline,
            event=event,
            baseline_minus_event_mag=gap(baseline, event, 5, 3),
        )
    elif route == "recurrent":
        a, b, c, d = map(date.fromisoformat, case["eligible"][route]["example_dates"])
        early = group(data, [t for t in data if a <= t <= b])
        late = group(data, [t for t in data if c <= t <= d])
        other = group(data, [t for t in data if not (a <= t <= b or c <= t <= d)])
        result.update(
            early=early,
            late=late,
            outside=other,
            outside_minus_early_mag=gap(other, early, 10, 2),
            outside_minus_late_mag=gap(other, late, 10, 2),
        )
    elif route == "fading":
        bright = group(data, signs[1] - signs[-1])
        faint = group(data, signs[-1] - signs[1])
        result.update(
            positive=bright,
            negative=faint,
            negative_minus_positive_mag=gap(faint, bright, 5, 5),
        )
    return result


def run(snapshot, ledger, sources, fetch=False):
    anchors, histories, _ = compass_archive.load(snapshot)
    selected = [
        r for r in json.loads(ledger.read_text())["objects"] if r["selected_for"]
    ]
    raw_anchors = {pilot.oid(a): a for a in anchors}
    result = {}
    for r in selected:
        name = r["id"]
        a = raw_anchors[name]
        measures = compass.measure(a, histories[name])
        if (
            str(a["i:candid"]) != r["candid"]
            or a["i:jd"] != r["jd"]
            or a["i:ra"] != r["ra"]
            or a["i:dec"] != r["dec"]
            or r["eligible"] != measures["eligible"]
            or r["selected_for"] not in r["eligible"]
        ):
            raise ValueError(f"Ledger and source snapshot disagree for {name}")
        series, rows, count = direct_data(r, sources / name, fetch)
        result[name] = {
            "route": r["selected_for"],
            "cone_series": count,
            "chosen_series": series,
            "bands": {
                band: compare(r, histories[name], rows, fid)
                for fid, band in compass.BANDS.items()
            },
        }
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("snapshot", "ledger", "sources"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = run(args.snapshot, args.ledger, args.sources, args.fetch)
    pilot.json_write(args.out, result)
