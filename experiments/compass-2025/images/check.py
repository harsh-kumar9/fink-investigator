"""Replay the compass case image triplets, Gaia cones and signed residual positions."""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parents[1]))
import compass_archive
import cutouts
import pilot


def quality(r):
    return (
        r["i:fid"] in (1, 2)
        and (r["i:drb"] or 0) >= 0.8
        and 0 < (r["i:sigmapsf"] or 99) <= 0.2
    )


def positive(r):
    sign = str(r["i:isdiffpos"]).lower()
    if sign in {"t", "true", "1"}:
        return True
    if sign in {"f", "false", "0"}:
        return False
    return None


def image_case(anchor):
    return {
        "id": pilot.oid(anchor),
        "jd": anchor["i:jd"],
        "filter": anchor["i:fid"],
        "candid": str(anchor["i:candid"]),
    }


def stamps(rows, anchor, folder, fetch):
    case = image_case(anchor)
    if fetch and not folder.exists():
        cutouts.capture_images(case, rows, folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    expected = cutouts.epochs(rows, case)
    if manifest["id"] != case["id"] or len(expected) != len(manifest["epochs"]):
        raise ValueError("Invalid cutout manifest")
    for r, e in zip(expected, manifest["epochs"], strict=True):
        if (
            e["candid"] != str(r["i:candid"])
            or e["jd"] != r["i:jd"]
            or e["isdiffpos"] != r["i:isdiffpos"]
            or [v["kind"] for v in e["images"]] != ["Science", "Template", "Difference"]
        ):
            raise ValueError("Cutout selection changed")
        for item in e["images"]:
            name = f"{e['candid']}-{item['kind'].lower()}.png"
            url = pilot.image_url(
                {"id": case["id"], "candid": e["candid"]}, item["kind"]
            )
            if (
                name != item["file"]
                or url != item["url"]
                or hashlib.sha256((folder / name).read_bytes()).hexdigest()
                != item["sha256"]
            ):
                raise ValueError("Unexpected cutout or checksum")
    return [e["candid"] for e in manifest["epochs"]]


def gaia(anchor, fetch):
    query = (
        "SELECT TOP 2000 source_id,ra,dec,ref_epoch,pmra,pmdec,pmra_error,pmdec_error,phot_g_mean_mag "
        "FROM gaia_dr3_source WHERE " + pilot.cone_clause([anchor], 60)
    )
    url = (
        pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "CSV", "QUERY": query}
        )
    )
    folder = HERE / pilot.oid(anchor)
    receipt = folder / "gaia.json"
    data = folder / "gaia.csv.gz"
    if fetch and not receipt.exists():
        raw = pilot.curl(url, [])
        pilot.parse_csv(raw)
        data.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            },
        )
    raw = gzip.decompress(data.read_bytes())
    m = json.loads(receipt.read_text())
    if m["url"] != url or m["sha256_uncompressed"] != hashlib.sha256(raw).hexdigest():
        raise ValueError("Gaia archive mismatch")
    stars = pilot.parse_csv(raw)
    if len({r["source_id"] for r in stars}) != len(stars) or any(
        sep(anchor, r, False) > 60.1 for r in stars
    ):
        raise ValueError("Gaia cone mismatch")
    return stars


def sep(row, star, propagate):
    return math.hypot(*relative_position(row, star, propagate))


def relative_position(row, star, propagate):
    x, y = pilot.offsets(
        row["i:ra"], row["i:dec"], float(star["ra"]), float(star["dec"])
    )
    if propagate and star["pmra"] and star["pmdec"]:
        years = (row["i:jd"] - 2451545) / 365.25 + 2000 - float(star["ref_epoch"])
        x -= years * float(star["pmra"]) / 1000
        y -= years * float(star["pmdec"]) / 1000
    return x, y


def case_summary(anchor, history, fetch):
    oid = pilot.oid(anchor)
    folder = HERE / oid
    rows = [
        r
        for r in history
        if r["i:jd"] <= anchor["i:jd"] and quality(r) and positive(r) is not None
    ]
    image_ids = stamps(history, anchor, folder / "epochs", fetch)
    negatives = [
        r for r in rows if positive(r) is False and r["i:fid"] == anchor["i:fid"]
    ]
    if negatives and not any(str(r["i:candid"]) in image_ids for r in negatives):
        last = max(negatives, key=lambda r: (r["i:jd"], r["i:candid"]))
        image_ids += stamps([last], last, folder / "negative", fetch)
    stars = gaia(anchor, fetch)
    nearby = []
    for s in stars:
        if min(sep(anchor, s, False), sep(anchor, s, True)) > 5:
            continue
        obj = {
            "id": s["source_id"],
            "g_mag": float(s["phot_g_mean_mag"]) if s["phot_g_mean_mag"] else None,
            "ref_epoch": float(s["ref_epoch"]),
            "ref_ra": float(s["ra"]),
            "ref_dec": float(s["dec"]),
            "proper_motion_applied": bool(s["pmra"] and s["pmdec"]),
            "ref_separation_arcsec": round(sep(anchor, s, False), 3),
            "anchor_separation_arcsec": round(sep(anchor, s, True), 3),
        }
        for sign, label in [(True, "positive"), (False, "negative")]:
            groups = defaultdict(list)
            for row in rows:
                if positive(row) == sign:
                    groups[pilot.day_of(row["i:jd"])].append(
                        relative_position(row, s, True)
                    )
            vals = sorted(
                median(math.hypot(*point) for point in night)
                for night in groups.values()
            )
            obj[label] = {
                "date_count": len(vals),
                "median_xy_arcsec": (
                    [
                        round(
                            median(
                                median(point[i] for point in night)
                                for night in groups.values()
                            ),
                            3,
                        )
                        for i in (0, 1)
                    ]
                    if vals
                    else None
                ),
                "median_distance_arcsec": round(median(vals), 3) if vals else None,
                "p90_distance_arcsec": (
                    round(vals[math.ceil(0.9 * len(vals)) - 1], 3) if vals else None
                ),
            }
        nearby.append(obj)
    moving = pilot.stellar_matches(
        anchor, {s["source_id"]: s for s in stars if s["pmra"] and s["pmdec"]}
    )
    return {
        "id": oid,
        "checked_image_candids": image_ids,
        "gaia_rows_60arcsec": len(stars),
        "nearby_5arcsec": sorted(nearby, key=lambda r: r["anchor_separation_arcsec"]),
        "moving_star_tracks": moving,
    }


def run(fetch):
    packet = json.loads((RUN / "results/packet.json").read_text())
    anchors, histories, _ = compass_archive.load(RUN / "snapshot")
    if not packet["gate_passed"]:
        raise ValueError("Review packet closed")
    by_id = {pilot.oid(r): r for r in anchors}
    return [
        case_summary(by_id[c["id"]], histories[c["id"]], fetch) for c in packet["cases"]
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    results = run(args.fetch)
    if args.out:
        pilot.json_write(args.out, results)
    else:
        print(json.dumps(results, indent=2))
