"""Resolve the selected astrometric displacement by ZTF detector and reference."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

import compass
import compass_archive
import compass_astrometry
import compass_direct
import pilot


def field_photometry(anchor, fid, sources, fetch):
    ra, dec = anchor["i:ra"], anchor["i:dec"]
    query = (
        "SELECT TOP 100 oid,ra,dec,field,ccdid,qid,fid,nobs FROM ztf_objects_dr24 "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{1.5/3600:.12f}))=1 AND fid={fid}"
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
        )
    )
    cone = compass_direct.capture(
        sources, "direct-cone", url, {"oid", "ra", "dec", "field", "nobs", "fid"}, fetch
    )
    if len(cone) >= 100:
        raise ValueError("DR24 cone hit cap")
    groups = defaultdict(list)
    for r in cone:
        sep = math.hypot(*pilot.offsets(float(r["ra"]), float(r["dec"]), ra, dec))
        if sep > 1.501 or int(r["fid"]) != fid or not r["oid"].isdigit():
            raise ValueError("Bad cone row")
        groups[int(r["field"])].append((int(r["nobs"]), sep, r["oid"]))
    series = {
        k: min(rows, key=lambda v: (-v[0], v[1], int(v[2])))[2]
        for k, rows in groups.items()
    }
    if not series:
        return {}, []
    url = (
        compass_direct.SERVICE
        + "?"
        + urlencode(
            [("ID", oid) for _, oid in sorted(series.items())]
            + [
                ("COLLECTION", "ztf_dr24"),
                ("FORMAT", "csv"),
                ("BAD_CATFLAGS_MASK", "65535"),
            ]
        )
    )
    raw = compass_direct.capture(
        sources,
        "direct",
        url,
        {
            "oid",
            "expid",
            "filtercode",
            "mjd",
            "mag",
            "magerr",
            "catflags",
            "ccdid",
            "qid",
            "ra",
            "dec",
            "field",
        },
        fetch,
    )
    band = compass.BANDS[fid]
    if any(
        r["oid"] not in series.values() or r["filtercode"] != "z" + band for r in raw
    ):
        raise ValueError("Unexpected recovered filter or object ID")
    return {
        k: (
            oid,
            compass_direct.nightly(
                [r for r in raw if r["oid"] == oid], band, anchor["i:jd"]
            )[0],
        )
        for k, oid in series.items()
    }, raw


def matched_positions(by_candid, packet, photometry):
    """Compare identical pairs of exposures with both position estimators."""
    candids = {
        c
        for pairs in packet["pairs"].values()
        for p in pairs
        for c in p["positive_candids"] + p["negative_candids"]
    }
    edges = []
    for r in photometry:
        if int(r["catflags"]) or not 0 < float(r["magerr"]) <= 0.2:
            continue
        t = float(r["mjd"]) + 2400000.5
        rcid = 4 * (int(r["ccdid"], 0) - 1) + int(r["qid"], 0) - 1
        position = (float(r["ra"]), float(r["dec"]))
        if not all(map(math.isfinite, position)):
            continue
        for c in candids:
            row = by_candid[c]
            dt = abs(t - row["i:jd"]) * 86400
            if dt <= 60 and int(r["field"]) == row["i:field"] and rcid == row["i:rcid"]:
                edges.append((dt, c, r["expid"], position))
    matched, used = {}, set()
    for _, c, e, xy in sorted(edges):
        if c not in matched and e not in used:
            matched[c] = xy
            used.add(e)
    band = packet["case"]["best_band"]
    dx, dy = packet["case"]["bands"][band]["training_east_north_arcsec"]
    scale = math.hypot(dx, dy)
    result = {}
    for split, pairs in packet["pairs"].items():
        direct, difference, dates = [], [], []
        for pair in pairs:
            sides = [pair["positive_candids"], pair["negative_candids"]]
            if any(c not in matched for cs in sides for c in cs):
                continue
            origin = packet["case"]
            points = [
                [pilot.offsets(*matched[c], origin["ra"], origin["dec"]) for c in cs]
                for cs in sides
            ]
            direct.append(
                [
                    median(v[i] for v in points[0]) - median(v[i] for v in points[1])
                    for i in (0, 1)
                ]
            )
            difference.append([pair["east"], pair["north"]])
            dates.append([pair["positive"], pair["negative"]])
        out = {"matched_pairs": len(direct), "original_pairs": len(pairs)}
        out["paired_projections"] = [
            {
                "utc_dates_positive_negative": d,
                "direct_arcsec": (dx * u[0] + dy * u[1]) / scale,
                "difference_arcsec": (dx * v[0] + dy * v[1]) / scale,
            }
            for d, u, v in zip(dates, direct, difference)
        ]
        for name, vecs in [("direct", direct), ("difference", difference)]:
            if vecs:
                x, y = [median(v[i] for v in vecs) for i in (0, 1)]
                out[name] = {
                    "east_north_arcsec": [x, y],
                    "training_direction_projection_arcsec": (dx * x + dy * y) / scale,
                }
        result[split] = out
    return result


def run(snapshot, packet_path, sources, fetch=False):
    packet = json.loads(packet_path.read_text())
    case = packet["case"]
    if case is None:
        return {"status": "No astrometry case admitted"}
    anchors, histories, _ = compass_archive.load(snapshot)
    name = case["id"]
    a = next(r for r in anchors if pilot.oid(r) == name)
    fid = next(k for k, v in compass.BANDS.items() if v == case["best_band"])
    stats, pairs = compass_astrometry.assess(a, histories[name], fid)
    if (
        str(a["i:candid"]) != case["candid"]
        or stats["status"] != "passed"
        or any(a["i:" + k] != case[k] for k in ("ra", "dec", "jd"))
    ):
        raise ValueError("Astrometry packet differs from the frozen source")
    compass.check_replay(pairs, packet["pairs"], "astrometry pairs")
    compass.check_replay(
        stats, case["bands"][case["best_band"]], "astrometry measurements"
    )
    cols = (
        pilot.COLUMNS
        + ",i:field,i:rcid,i:jdstartref,i:jdendref,i:fwhm,i:ranr,i:decnr,i:magnr"
    )
    params = {"objectId": name, "columns": cols, "output-format": "json"}
    manifest = sources / "manifest.json"
    if fetch:
        sources.mkdir(parents=True, exist_ok=True)
    meta = json.loads(manifest.read_text()) if manifest.exists() else {"sources": {}}
    if not fetch and "history.json.gz" not in meta["sources"]:
        raise ValueError("Reference history not archived")
    archive = pilot.store_raw(
        sources, meta, "history.json.gz", pilot.FINK + "/objects", params
    )
    by_candid = {str(r["i:candid"]): r for r in pilot.histories(archive)[name]}
    groups = defaultdict(list)
    key_by_candid = {}
    for old in histories[name]:
        if old["i:fid"] != fid or old["i:jd"] > a["i:jd"] or not compass.quality(old):
            continue
        cand = str(old["i:candid"])
        row = by_candid[cand]
        if any(
            row[k] != old[k]
            for k in ("i:ra", "i:dec", "i:jd", "i:fid", "i:isdiffpos", "i:objectId")
        ):
            raise ValueError("Measurement changed between archive queries")
        key = tuple(row[k] for k in ("i:field", "i:rcid", "i:jdstartref", "i:jdendref"))
        if any(v is None for v in key):
            raise ValueError("Incomplete detector/reference identification")
        key_by_candid[cand] = key
        groups[key].append(row)
    phot, raw_direct = field_photometry(a, fid, sources, fetch)
    pg = packet["pairs"]
    pair_stats = {}
    for split, pairs in pg.items():
        counts = defaultdict(int)
        for p in pairs:
            pk = {key_by_candid[i] for i in p["positive_candids"]}
            nk = {key_by_candid[i] for i in p["negative_candids"]}
            if len(pk) == len(nk) == 1 and pk == nk:
                counts["same_reference"] += 1
            elif pk.isdisjoint(nk):
                counts["disjoint_references"] += 1
            else:
                counts["mixed"] += 1
        pair_stats[split] = dict(counts)
    summaries = []
    for k, rows in sorted(groups.items()):
        dates = compass_astrometry.daily(a, rows, fid)
        by_sign = {}
        oid, direct = phot.get(k[0], (None, {}))
        for s in [1, -1]:
            positions = [v for (d, sgn), v in dates.items() if sgn == s]
            ds = {d for d, sgn in dates if sgn == s} - {
                d for d, sgn in dates if sgn == -s
            }
            vals = [r for r in rows if compass.sign(r) == s]
            by_sign[str(s)] = {
                "dates": len(positions),
                "median_offset_arcsec": (
                    [median(x[i] for x in positions) for i in [0, 1]]
                    if positions
                    else None
                ),
                "median_psf_fwhm_pixels": (
                    median(r["i:fwhm"] for r in vals) if vals else None
                ),
                "direct": compass_direct.group(direct, ds),
            }
        matched, _ = compass_astrometry.assess(a, rows, fid)
        summaries.append(
            {
                "field": k[0],
                "rcid": k[1],
                "ref_start_jd": k[2],
                "ref_end_jd": k[3],
                "series_id": oid,
                "signs": by_sign,
                "within_reference": matched,
                "direct_negative_minus_positive_mag": compass_direct.gap(
                    by_sign["-1"]["direct"], by_sign["1"]["direct"], 5, 5
                ),
            }
        )
    return {
        "id": name,
        "original_pairs_by_reference": pair_stats,
        "reference_groups": summaries,
        "same_exposure_positions": matched_positions(by_candid, packet, raw_direct),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("snapshot", "packet", "sources", "out"):
        p.add_argument("--" + key, required=True, type=Path)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    pilot.json_write(
        args.out, run(args.snapshot, args.packet, args.sources, args.fetch)
    )
