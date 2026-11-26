"""Measure the early and late bright dates of J0722 against same-exposure stars."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json
import math
from pathlib import Path
from statistics import median
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
import compass_direct as direct
import pilot

NAME = "ZTF25aacfjde"
JAN = date(2025, 1, 19)
APR = date(2025, 4, 16)
START = date(2024, 1, 1)
NMW = 2460696.4416
BANDS = ("g", "r")


def distance(row, other):
    return math.hypot(
        *pilot.offsets(*[float(x[k]) for x in (row, other) for k in ("ra", "dec")])
    )


def curves(rows):
    sources = defaultdict(dict)
    for row in rows:
        exp = row["expid"]
        t, mag, err = [float(row[k]) for k in ("mjd", "mag", "magerr")]
        if (
            int(row["catflags"])
            or not 0 < err <= 0.2
            or not all(map(math.isfinite, (t, mag)))
        ):
            continue
        if exp in sources[row["oid"]] and sources[row["oid"]][exp] != (t, mag, err):
            raise ValueError("Conflicting exposures")
        sources[row["oid"]][exp] = (t, mag, err)
    return sources


def nights(data):
    values = defaultdict(list)
    for mjd, mag, _ in data.values():
        values[pilot.day_of(mjd + 2400000.5)].append(mag)
    return {d: median(v) for d, v in values.items()}


def window(series, start, stop):
    vals = [v for d, v in series.items() if start <= d <= stop]
    return {"n": len(vals), "median": median(vals) if vals else None}


def deltas(series):
    out = {}
    for label, day, lo, hi in [
        ("january", JAN, START, JAN - timedelta(days=1)),
        ("previous_month", JAN, JAN - timedelta(days=30), JAN - timedelta(days=1)),
        ("april_before", APR, date(2025, 4, 5), date(2025, 4, 11)),
        ("april_after", APR, date(2025, 4, 21), date(2025, 5, 3)),
    ]:
        ref = window(series, lo, hi)
        ref["day_mag"] = series.get(day)
        ref["brightening_mag"] = (
            ref["median"] - series[day] if ref["n"] and day in series else None
        )
        out[label] = ref
    return out


def acquire(target, fid, fetch, bright_panel=False):
    cone_dir = HERE / "onset-field"
    query = (
        "SELECT TOP 2000 oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag FROM ztf_objects_dr24 "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{target['ra']},{target['dec']},{120/3600}))=1 "
        f"AND field={target['field']} AND ccdid={target['ccdid']} AND qid={target['qid']} AND fid={fid} AND nobs>=100"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    )
    fields = direct.capture(
        cone_dir,
        f"cone-{fid}",
        url,
        {"oid", "ra", "dec", "field", "ccdid", "qid", "fid", "nobs", "meanmag"},
        fetch,
    )
    if len(fields) >= 2000 or any(
        any(r[k] != target[k] for k in ("field", "ccdid", "qid", "fid")) for r in fields
    ):
        raise ValueError("Unexpected control cone")
    (center,) = [r for r in fields if r["oid"] == target["oid"]]
    match_mag = (
        (16.4 if fid == 1 else 16.9) if bright_panel else float(center["meanmag"])
    )
    controls = []
    for r in sorted(
        fields,
        key=lambda r: (
            abs(float(r["meanmag"]) - match_mag),
            distance(r, center),
            r["oid"],
        ),
    ):
        if (
            10 <= distance(r, center) <= 120
            and abs(float(r["meanmag"]) - match_mag) <= 1.5
            and all(distance(r, c) >= 5 for c in controls)
        ):
            controls.append(r)
            if len(controls) == 5:
                break
    if len(controls) < 3:
        raise ValueError("Too few field stars")
    label = "bright-stars" if bright_panel else "stars"
    return controls, photometry(
        cone_dir, f"{label}-{fid}", [r["oid"] for r in controls], fetch
    )


def photometry(folder, label, ids, fetch):
    params = [("ID", i) for i in ids] + [
        ("COLLECTION", "ztf_dr24"),
        ("FORMAT", "csv"),
        ("BAD_CATFLAGS_MASK", "65535"),
    ]
    rows = direct.capture(
        folder,
        label,
        direct.SERVICE + "?" + urlencode(params),
        {"oid", "mjd", "mag", "magerr", "catflags", "expid", "filtercode"},
        fetch,
    )
    if set(r["oid"] for r in rows) - set(ids):
        raise ValueError("Photometry outside requested objects")
    return curves(rows)


def field_correction(series, controls):
    offsets = defaultdict(list)
    for star in controls.values():
        zero = window(nights(star), START, JAN - timedelta(days=1))
        if zero["n"] < 10:
            raise ValueError("Insufficient field baseline")
        for exp, (t, mag, err) in star.items():
            offsets[exp].append(mag - zero["median"])
    corrected = {
        e: (t, v - median(offsets[e]), err)
        for e, (t, v, err) in series.items()
        if len(offsets[e]) >= 3
    }
    return corrected, offsets


def run(fetch=False):
    ledger = json.loads((HERE / "results/ledger.json").read_text())
    (case,) = [
        r for r in ledger["objects"] if r["id"] == NAME and r["selected_for"] == "onset"
    ]
    selected, rows, _ = direct.direct_data(case, HERE / "direct" / NAME, False)
    sources = curves(rows)
    primary = {r["oid"] for r in selected.values()}
    import csv, gzip

    with gzip.open(HERE / "direct" / NAME / "cone.csv.gz", "rt") as fh:
        more = [r for r in csv.DictReader(fh) if r["oid"] not in primary]
    alt = photometry(
        HERE / "onset-field", "other-field", sorted(r["oid"] for r in more), fetch
    )
    output = {"other_field": {}, "bands": {}}
    for row in more:
        series = nights(alt.get(row["oid"], {}))
        output["other_field"][row["oid"]] = {
            "fid": int(row["fid"]),
            "field": row["field"],
            "first": str(min(series)) if series else None,
            "last": str(max(series)) if series else None,
            "dates": len(series),
            "january": series.get(JAN),
            "april": series.get(APR),
        }
    all_nights = {}
    for fid, band in enumerate(BANDS, 1):
        target = selected[fid]
        points = sources[target["oid"]]
        ordered, stars = acquire(target, fid, fetch)
        corrected, offsets = field_correction(points, stars)
        bright_order, bright_stars = acquire(target, fid, fetch, bright_panel=True)
        bright_curve, bright_offsets = field_correction(points, bright_stars)
        series = nights(corrected)
        all_nights[band] = series
        summary = deltas(series)
        zero = summary["january"]["median"]
        returns = [d for d, v in series.items() if d > APR and abs(v - zero) <= 0.3]
        bright = sorted(
            d
            for d, v in series.items()
            if returns and d > min(returns) and zero - v >= 2
        )
        cutoff = max(series)
        events = []
        for exp, (t, mag, err) in sorted(points.items(), key=lambda x: x[1][0]):
            d = pilot.day_of(t + 2400000.5)
            if d in (JAN, APR):
                events.append(
                    {
                        "day": str(d),
                        "expid": exp,
                        "mjd": t,
                        "days_from_NMW": t + 2400000.5 - NMW,
                        "mag": mag,
                        "err": err,
                        "control_count": len(offsets[exp]),
                        "field_corrected_mag": (
                            corrected[exp][1] if exp in corrected else None
                        ),
                    }
                )
        output["bands"][band] = {
            "source_oid": target["oid"],
            "raw": deltas(nights(points)),
            "corrected": summary,
            "events": events,
            "controls": {r["oid"]: deltas(nights(stars[r["oid"]])) for r in ordered},
            "bright_panel": {
                "corrected": deltas(nights(bright_curve)),
                "comparison_oids": [r["oid"] for r in bright_order],
                "april_control_counts": [
                    len(bright_offsets[e["expid"]])
                    for e in events
                    if e["day"] == str(APR)
                ],
                "individual": {
                    r["oid"]: deltas(nights(bright_stars[r["oid"]]))
                    for r in bright_order
                },
            },
            "first_return_within_0.3": str(min(returns)) if returns else None,
            "last_observed": str(cutoff),
            "dates_after_first_return": (
                sum(d > min(returns) for d in series) if returns else 0
            ),
            "subsequent_bright_dates": list(map(str, bright)),
            "subsequent_pairs": [
                [str(a), str(b)]
                for a, b in zip(bright, bright[1:])
                if (b - a).days <= 14
            ],
            "event_context": [
                [str(d), v]
                for d, v in sorted(series.items())
                if JAN - timedelta(days=12) <= d <= JAN + timedelta(days=2)
                or date(2025, 4, 1) <= d <= date(2025, 5, 5)
            ],
        }
    both = sorted(d for d in all_nights["g"].keys() & all_nights["r"].keys() if d < JAN)
    as_bright = [
        str(d)
        for d in both
        if all(all_nights[b][d] <= all_nights[b][JAN] for b in BANDS)
    ]
    output["january_empirical"] = {
        "earlier_paired_dates": len(both),
        "as_bright_in_both": as_bright,
        "first": str(both[0]),
        "last": str(both[-1]),
    }
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    print(json.dumps(run(args.fetch), indent=2, allow_nan=False))
