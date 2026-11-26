"""Check repeat detectors, neighboring stars and the calendar of missing target exposures."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import csv, gzip, json, math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

import numpy as np
from measure import HERE, BATCH, NAMES, load, baseline, previous, compass_direct, pilot

SRC = HERE / "sources"


def quality(v):
    return (
        not int(v["catflags"])
        and 0 < float(v["magerr"]) <= 0.2
        and all(math.isfinite(float(v[k])) for k in ("mjd", "hjd", "mag"))
    )


def phot(folder, label, oids, fetch):
    params = [("ID", oid) for oid in sorted(oids, key=int)] + [
        ("COLLECTION", "ztf_dr24"),
        ("FORMAT", "csv"),
        ("BAD_CATFLAGS_MASK", "65535"),
    ]
    rows = compass_direct.capture(
        folder,
        label,
        compass_direct.SERVICE + "?" + urlencode(params),
        {
            "oid",
            "mjd",
            "expid",
            "hjd",
            "catflags",
            "magerr",
            "mag",
            "filtercode",
            "limitmag",
        },
        fetch,
    )
    if set(oids) != {v["oid"] for v in rows}:
        raise ValueError("DR24 series missing from photometry response")
    return rows


def offset(a, case):
    return math.hypot(
        *pilot.offsets(float(a["ra"]), float(a["dec"]), case["ra"], case["dec"])
    )


def splits(d, boundary):
    v = (d - baseline.EPOCH).days - boundary
    return "before" if -180 <= v < 0 else "after" if 0 <= v < 180 else "outside"


def summarize(curve, boundary, anchor, exclude=()):
    clean = {d: m for d, (t, m) in curve.items() if d not in exclude}
    outside = [m for d, m in clean.items() if splits(d, boundary) == "outside"]
    if len(outside) < 5:
        return {"outside_n": len(outside), "status": "insufficient outside dates"}
    bright, mid, faint = np.quantile(outside, [0.2, 0.5, 0.8])
    result = {
        "outside_n": len(outside),
        "outside_q20_q50_q80": [float(bright), float(mid), float(faint)],
        "windows": {},
    }
    for s in ("before", "after"):
        dates = sorted(d for d in clean if splits(d, boundary) == s)
        vals = [clean[d] for d in dates]
        bins = defaultdict(list)
        for d in dates:
            bins[((d - baseline.EPOCH).days - 58000) // 30].append(clean[d])
        counts = []
        for vs in bins.values():
            if len(vs) >= 3:
                counts.append(
                    [np.mean(np.array(vs) <= bright), np.mean(np.array(vs) >= faint)]
                )
        result["windows"][s] = {
            "n": len(dates),
            "first": str(dates[0]) if dates else None,
            "last": str(dates[-1]) if dates else None,
            "max_internal_gap": max(
                ((b - a).days for a, b in zip(dates, dates[1:])), default=None
            ),
            "q20_q50_q80": (
                list(np.round(np.quantile(vals, [0.2, 0.5, 0.8]), 5)) if vals else None
            ),
            "bright": sum(m <= bright for m in vals).item() if vals else 0,
            "faint": sum(m >= faint for m in vals).item() if vals else 0,
            "complete_bins_with_3_nights": len(counts),
            "equal_bin_weight_bright_faint_fractions": (
                list(np.round(np.mean(counts, axis=0), 5)) if counts else None
            ),
        }
    curve_used = {d: v for d, v in curve.items() if d in clean}
    result["bin_contrast"] = baseline.change(
        baseline.bin_curve(curve_used, anchor), boundary
    )
    return result


def field_case(name, fetch):
    case, pre, _ = load(name)
    folder = SRC / name
    anchor = date.fromisoformat(case["date"])
    boundary = json.loads((HERE / "results.json").read_text())[name]["state"][
        "boundary_mjd"
    ]
    selected, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / name)
    with gzip.open(BATCH / "direct" / name / "cone.csv.gz", "rt") as f:
        cone = list(csv.DictReader(f))
    alt = [r for r in cone if r["oid"] not in {v["oid"] for v in selected.values()}]
    alt_rows = phot(folder, "alternates", [r["oid"] for r in alt], fetch) if alt else []
    columns = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag"
    q = (
        f"SELECT TOP 10000 {columns} FROM ztf_objects_dr24 WHERE CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{case['ra']},{case['dec']},{90/3600}))=1 AND fid IN (1,2) AND nobs>=40"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
    )
    field = compass_direct.capture(folder, "field", url, set(columns.split(",")), fetch)
    pool = {r["oid"]: r for r in field}
    if len(field) >= 10000 or len(pool) != len(field):
        raise ValueError("Truncated or duplicated field")
    out = {
        "state_boundary": str(baseline.EPOCH + timedelta(days=boundary)),
        "primary": {},
        "other_series": {},
    }
    for fid, target in selected.items():
        b = {1: "g", 2: "r"}[fid]
        oid = target["oid"]
        typical = float(pool[oid]["meanmag"])
        other = sorted(
            (
                r
                for r in field
                if r["oid"] != oid
                and all(r[k] == target[k] for k in ("fid", "field", "ccdid", "qid"))
                and 10 <= offset(r, case) <= 90
                and r["meanmag"]
                and abs(float(r["meanmag"]) - typical) <= 1
            ),
            key=lambda r: (offset(r, case), int(r["oid"])),
        )[:3]
        summary = summarize(pre[b], boundary, anchor)
        result = {"primary_oid": oid, "uncorrected": summary, "field_objects": other}
        out["primary"][b] = result
        if len(other) < 2:
            continue
        refs = phot(folder, b + "-controls", [r["oid"] for r in other], fetch)
        counts = defaultdict(dict)
        levels = {}
        for c in other:
            vv = [
                v
                for v in refs
                if v["oid"] == c["oid"]
                and quality(v)
                and float(v["mjd"]) + 2400000.5 <= case["jd"] + 1e-7
            ]
            outside = [
                float(v["mag"])
                for v in vv
                if splits(pilot.day_of(float(v["mjd"]) + 2400000.5), boundary)
                == "outside"
            ]
            if not outside:
                continue
            levels[c["oid"]] = median(outside)
            for v in vv:
                counts[v["expid"]][c["oid"]] = v
        target_rows = {v["expid"]: v for v in rows if v["oid"] == oid and quality(v)}
        grouped = defaultdict(list)
        corrected = []
        status = defaultdict(dict)
        for expid, rr in counts.items():
            if len(rr) < 2:
                continue
            tt = float(next(iter(rr.values()))["mjd"])
            day = pilot.day_of(tt + 2400000.5)
            key = splits(day, boundary)
            if key != "outside":
                grouped[day].append((-len(rr), tt, expid))
            if expid not in target_rows:
                continue
            v = target_rows[expid]
            cal = median(float(m["mag"]) - levels[k] for k, m in rr.items())
            corrected.append({**v, "mag": float(v["mag"]) - cal})
        result["field_corrected"] = summarize(
            previous.nightly(corrected, case["jd"] + 1e-7)[b], boundary, anchor
        )
        bright = summary["outside_q20_q50_q80"][0]
        for d, v in sorted(grouped.items()):
            expid = min(v)[2]
            rr = counts[expid]
            point = target_rows.get(expid)
            lim = median(float(r["limitmag"]) for r in rr.values())
            status[splits(d, boundary)][str(d)] = {
                "expid": expid,
                "limitmag": round(lim, 5),
                "target_mag": float(point["mag"]) if point else None,
                "n_controls": len(rr),
                "above_bright_threshold": (
                    float(point["mag"]) <= bright if point else None
                ),
            }
        coverage = {}
        for k, by_day in status.items():
            dd = list(by_day.values())
            n = len(dd)
            seen = sum(r["target_mag"] is not None for r in dd)
            br = sum(r["above_bright_threshold"] is True for r in dd)
            coverage[k] = {
                "field_dates": n,
                "target_dates": seen,
                "bright_dates": br,
                "bright_fraction_range": [
                    round(br / n, 5),
                    round((br + n - seen) / n, 5),
                ],
                "median_limitmag": median(r["limitmag"] for r in dd),
                "field_calendar": by_day,
            }
        result["coverage"] = coverage
        result["control_series"] = {
            c["oid"]: summarize(
                previous.nightly(
                    [v for v in refs if v["oid"] == c["oid"]], case["jd"] + 1e-7
                )[b],
                boundary,
                anchor,
            )
            for c in other
        }
    for v in alt:
        b = {"1": "g", "2": "r"}[v["fid"]]
        arr = previous.nightly(
            [r for r in alt_rows if r["oid"] == v["oid"]], case["jd"] + 1e-7
        )[b]
        out["other_series"][v["oid"]] = {
            "band": b,
            "field": v["field"],
            "all_dates": summarize(arr, boundary, anchor),
            "absent_primary_dates": summarize(arr, boundary, anchor, pre[b]),
        }
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--name", choices=NAMES)
    p.add_argument("--out", type=Path, default=HERE / "field.json")
    args = p.parse_args()
    pilot.json_write(
        args.out,
        {n: field_case(n, args.fetch) for n in ([args.name] if args.name else NAMES)},
    )
