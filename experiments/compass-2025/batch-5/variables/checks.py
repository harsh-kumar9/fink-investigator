"""Check timing shifts, independent flux and nearby-star sampling for MGAB-V960."""

import argparse
from collections import defaultdict
from datetime import date, timedelta
import json, math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

import numpy as np
from scipy.stats import spearmanr

from analyze_vars import HERE, BATCH, prior, previous, compass_direct, load, pilot
from cadence import NAME
from duration import SRC, archive
from alternating import group_peaks
from tess_data import bin6, cube, sectors, SELECTION
from tess_events import corrected


def imagetest():
    case, pre, _ = load(NAME)
    options = sectors(case, False)
    out = {}
    groups = [g for g in group_peaks() if "2019-07-15" <= g["first"] <= "2019-10-31"]
    peaks = [
        pilot.utc_jd(date.fromisoformat(g["r"]["date"] if g["r"] else g["first"])) + 0.5
        for g in groups
    ]
    for s in SELECTION[NAME]:
        m, rt, p = cube(case, s, *options[s], False)
        arr = bin6(rt, p)
        t = arr[:, 0]
        train = np.ones(len(t), bool)
        for peak in peaks:
            train &= abs(t - peak) > 1.5
        passed = []
        for j in range(100):
            v = arr.copy()
            v[:, 1:] = np.roll(v[:, 1:], round(j * len(v) / 100), axis=0)
            imgs = corrected(v, train)[:, 1:].reshape(-1, 9, 9)
            success = True
            for b in "rg":
                pairs = [
                    (k, val)
                    for d, (target_time, val) in pre[b].items()
                    for k in [np.argmin(abs(t - target_time))]
                    if abs(t[k] - target_time) <= 0.25
                ]
                if len(pairs) < 6:
                    success = False
                    break
                k, y = zip(*pairs)
                k = np.array(k)
                y = np.array(y)
                order = np.argsort(y)
                n = len(y) // 3
                diff = np.median(imgs[k[order[:n]]], axis=0) - np.median(
                    imgs[k[order[-n:]]], axis=0
                )
                py, px = np.unravel_index(diff.argmax(), diff.shape)
                x0, y0 = m["target_xy"]
                if math.hypot(px - x0, py - y0) > 1:
                    success = False
                    break
                if spearmanr(imgs[k, 3:6, 3:6].sum(axis=(1, 2)), -y).statistic < 0.65:
                    success = False
                    break
            if success:
                passed.append(j)
        out[str(s)] = {"shifts_passing_both_bands": passed, "trials": 100}
    common = sorted(
        set.intersection(*[set(v["shifts_passing_both_bands"]) for v in out.values()])
    )
    return {"per_sector": out, "same_fraction_passing_all": common, "trials": 100}


def field(fetch):
    case, pre, _ = load(NAME)
    anchor = date.fromisoformat(case["date"])
    selected, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    q = (
        f"SELECT TOP 10000 oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag FROM ztf_objects_dr24 WHERE CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{case['ra']},{case['dec']},{90/3600}))=1 AND fid IN (1,2) AND nobs>=40"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q})
    )
    table = compass_direct.capture(
        SRC,
        "v960-control-field",
        url,
        {"oid", "ra", "dec", "field", "fid", "meanmag"},
        fetch,
    )
    if len(table) >= 10000:
        raise ValueError("Field truncated")
    edge = json.loads((HERE / "results.json").read_text())[NAME]["state"][
        "boundary_mjd"
    ]

    def state(data):
        return {
            k.replace("mag", "flux"): v
            for k, v in prior.change(prior.bin_curve(data, anchor), edge).items()
        }

    def span(d):
        return ((d - prior.EPOCH).days - edge) // 180

    _, independent = archive(case, False)
    out = {
        "boundary": str(prior.EPOCH + timedelta(days=edge)),
        "asas_bin_flux_change_mjy": state(
            {d: (pilot.utc_jd(d), v["flux"]) for d, v in independent.items()}
        ),
        "bands": {},
    }
    for fid, target in selected.items():
        b = {1: "g", 2: "r"}[fid]
        target_row = next(v for v in table if v["oid"] == target["oid"])
        cands = []
        for v in table:
            if (
                any(v[k] != target[k] for k in ("field", "ccdid", "qid", "fid"))
                or not v["meanmag"]
            ):
                continue
            sep = math.hypot(
                *pilot.offsets(float(v["ra"]), float(v["dec"]), case["ra"], case["dec"])
            )
            if (
                10 <= sep <= 90
                and abs(float(v["meanmag"]) - float(target_row["meanmag"])) <= 1
            ):
                cands.append((sep, v["oid"]))
        chosen = sorted(cands)[:3]
        params = [("ID", i) for _, i in chosen] + [
            ("COLLECTION", "ztf_dr24"),
            ("FORMAT", "csv"),
            ("BAD_CATFLAGS_MASK", "65535"),
        ]
        raw = compass_direct.capture(
            SRC,
            "v960-" + b + "-refs",
            compass_direct.SERVICE + "?" + urlencode(params),
            {"oid", "expid", "mag", "catflags", "magerr"},
            fetch,
        )

        def good(v):
            return (
                int(v["catflags"]) == 0
                and 0 < float(v["magerr"]) <= 0.2
                and float(v["mjd"]) + 2400000.5 <= case["jd"] + 1e-7
            )

        group = defaultdict(dict)
        for r in raw:
            if good(r):
                group[r["expid"]][r["oid"]] = r
        tr = {r["expid"]: r for r in rows if r["oid"] == target["oid"] and good(r)}
        before = {}
        after = {}
        for exp, rr in group.items():
            if len(rr) < 2:
                continue
            v = next(iter(rr.values()))
            d = pilot.day_of(float(v["mjd"]) + 2400000.5)
            dest = before if span(d) == -1 else after if span(d) == 0 else None
            if dest is None:
                continue
            candidate = (-len(rr), float(v["mjd"]), exp)
            if d not in dest or candidate < dest[d]:
                dest[d] = candidate
        external = [mag for d, (t, mag) in pre[b].items() if span(d) not in (-1, 0)]
        cut = float(np.quantile(external, 0.2))

        def coverage(data):
            vals = [float(tr[exp]["mag"]) for num, t, exp in data.values() if exp in tr]
            br = sum(u <= cut for u in vals)
            n = len(data)
            return {
                "field_dates": n,
                "matched": len(vals),
                "bright": br,
                "bright_fraction_bounds": (
                    [round(br / n, 4), round((br + n - len(vals)) / n, 4)]
                    if n
                    else None
                ),
            }

        controls = {}
        for sep, oid in chosen:
            a = previous.nightly(
                [r for r in raw if r["oid"] == oid], case["jd"] + 1e-7
            )[b]
            controls[oid] = {
                "separation_arcsec": round(sep, 3),
                "contrast_mag": prior.change(prior.bin_curve(a, anchor), edge)[
                    "after_minus_before_mag"
                ],
            }
        out["bands"][b] = {
            "controls": controls,
            "before": coverage(before),
            "after": coverage(after),
            "outside_bright_threshold_mag": round(cut, 5),
        }
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "checks.json")
    a = p.parse_args()
    pilot.json_write(a.out, {"state": field(a.fetch), "tess_shifts": imagetest()})
