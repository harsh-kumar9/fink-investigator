"""Replay ASAS-SN Sky Patrol forced flux on dates of the long red-selected group."""

import argparse
from base64 import encodebytes
from collections import defaultdict
from datetime import date, datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import median

import numpy as np
import pyarrow.parquet as pq
import requests
from scipy.stats import spearmanr

from analyze import HERE, curves, pilot
from interval import NAME, START, END, SRC

RADIUS = 6 / 3600
RA, DEC = 99.526908, 45.2491446
HOST = "http://asassn-lb01.ifa.hawaii.edu:9006"
COLS = ["asas_sn_id", "ra_deg", "dec_deg", "catalog_sources"]
POST = {"catalog": "master_list", "cols": COLS, "format": "arrow", "download": True}


def obtain(fetch, local, neighbors=False, ra=RA, dec=DEC, stem=None, sources=SRC):
    stem = stem or ("asas-field" if neighbors else "asas")
    radius = 60 / 3600 if neighbors else RADIUS
    cone_url = f"{HOST}/lookup_cone/radius{radius}_ra{ra}_dec{dec}"
    key = encodebytes(
        f"conectr-{radius}_conera-{ra}_conedec-{dec}|catalog-master_list|cols-{'/'.join(COLS)}".encode()
    ).decode()
    receipt = sources / (stem + ".json")
    cone = sources / (stem + "-cone.parquet")
    phot = sources / (stem + "-forced.parquet")
    if (fetch or local) and not receipt.exists():
        if local:
            servers = json.loads((local / "fink-mgab-asas-servers.json").read_text())
            cone_bytes = (local / "fink-mgab-asas-cone.msg").read_bytes()
        else:
            response = requests.get(HOST + "/get_block_servers", timeout=30)
            response.raise_for_status()
            servers = response.json()
            response = requests.post(cone_url, json=POST, timeout=45)
            response.raise_for_status()
            cone_bytes = response.content
        cone_table = pq.ParquetFile(io.BytesIO(cone_bytes)).read(use_threads=False)
        if set(COLS) - set(cone_table.schema.names):
            raise ValueError("Unexpected Sky Patrol catalog schema")
        if not len(cone_table):
            cone.write_bytes(cone_bytes)
            pilot.json_write(
                receipt,
                {
                    "cone_url": cone_url,
                    "cone_post": POST,
                    "empty": True,
                    "cone_sha256": hashlib.sha256(cone_bytes).hexdigest(),
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(
                        timespec="seconds"
                    ),
                },
            )
            return [], []
        server = servers[0]
        if not server.endswith(".ifa.hawaii.edu") or not all(
            c.isalnum() or c in ".-" for c in server
        ):
            raise ValueError("Unknown ASAS-SN block domain")
        url = f"http://{server}:9006/get_block/query_hash-{key}-block_idx-0-catalog-master_list"
        if local:
            data = (local / "fink-mgab-asas-data.parquet").read_bytes()
        else:
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            data = response.content
        pq.ParquetFile(io.BytesIO(cone_bytes)).read(use_threads=False)
        pq.ParquetFile(io.BytesIO(data)).read(use_threads=False)
        cone.write_bytes(cone_bytes)
        phot.write_bytes(data)
        pilot.json_write(
            receipt,
            {
                "cone_url": cone_url,
                "cone_post": POST,
                "servers_url": HOST + "/get_block_servers",
                "chosen_block_server": server,
                "block_url": requests.Request("GET", url).prepare().url,
                "cone_sha256": hashlib.sha256(cone_bytes).hexdigest(),
                "photometry_sha256": hashlib.sha256(data).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    meta = json.loads(receipt.read_text())
    if meta.get("empty"):
        if (
            meta["cone_url"] != cone_url
            or meta["cone_post"] != POST
            or hashlib.sha256(cone.read_bytes()).hexdigest() != meta["cone_sha256"]
            or len(pq.ParquetFile(cone).read(use_threads=False))
        ):
            raise ValueError("Empty Sky Patrol result failed replay")
        return [], []
    server = meta["chosen_block_server"]
    block = (
        requests.Request(
            "GET",
            f"http://{server}:9006/get_block/query_hash-{key}-block_idx-0-catalog-master_list",
        )
        .prepare()
        .url
    )
    if (
        meta["cone_url"] != cone_url
        or meta["cone_post"] != POST
        or meta["block_url"] != block
    ):
        raise ValueError("Sky Patrol request changed")
    for path, label in ((cone, "cone_sha256"), (phot, "photometry_sha256")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta[label]:
            raise ValueError("Sky Patrol source archive checksum changed")
    cone_rows = pq.ParquetFile(cone).read(use_threads=False).to_pylist()
    for r in cone_rows:
        r["separation_arcsec"] = round(
            math.hypot(*pilot.offsets(r["ra_deg"], r["dec_deg"], ra, dec)), 5
        )
    if len(cone_rows) > 100 or any(
        r["separation_arcsec"] > radius * 3600 + 0.001 for r in cone_rows
    ):
        raise ValueError("Sky Patrol cone is incomplete or has unexpected positions")
    target = min(cone_rows, key=lambda r: (r["separation_arcsec"], r["asas_sn_id"]))
    if target["separation_arcsec"] > 3:
        raise ValueError("No Sky Patrol counterpart within 3 arcsec")
    full = pq.ParquetFile(phot).read(use_threads=False).to_pylist()
    if not full:
        return cone_rows, []
    ids = {r["asas_sn_id"] for r in full}
    if (
        not ids.issubset({r["asas_sn_id"] for r in cone_rows})
        or target["asas_sn_id"] not in ids
    ):
        raise ValueError("Unexpected Sky Patrol IDs in the block")
    return cone_rows, (
        full
        if neighbors
        else [r for r in full if r["asas_sn_id"] == target["asas_sn_id"]]
    )


def run(fetch, local):
    cone, raw = obtain(fetch, local)
    case, pre = curves(NAME)
    anchor = date.fromisoformat(case["date"])
    thresholds = json.loads((HERE / "eruptions.json").read_text())[
        "training_thresholds_mag"
    ]["g"]
    groups = defaultdict(list)
    for r in raw:
        if (
            r["quality"] != "G"
            or not all(math.isfinite(r[k]) for k in ("jd", "flux", "flux_err"))
            or r["flux_err"] <= 0
        ):
            continue
        d = pilot.day_of(r["jd"])
        if date(2019, 1, 1) <= d <= anchor:
            groups[d].append(r)
    flux = {d: median(v["flux"] for v in vv) for d, vv in groups.items()}
    green = pre["g"]
    overlap = sorted(set(green) & set(flux))
    common = [d for d in overlap if START <= d <= END]
    training = [d for d in overlap if not START <= d <= END]
    bright = [d for d in training if green[d][1] <= thresholds["bright_q25"]]
    faint = [d for d in training if green[d][1] >= thresholds["faint_q75"]]
    basis = {
        k: round(median(flux[d] for d in dd), 6)
        for k, dd in (("bright", bright), ("faint", faint))
    }

    def rank(dd):
        return (
            round(
                float(
                    spearmanr(
                        [flux[d] for d in dd], [-green[d][1] for d in dd]
                    ).statistic
                ),
                5,
            )
            if len(dd) > 2
            else None
        )

    within = [
        {
            "date": str(d),
            "flux_mjy": round(flux[d], 6),
            "median_single_epoch_error_mjy": round(
                median(r["flux_err"] for r in groups[d]), 6
            ),
            "cameras": sorted({r["camera"] for r in groups[d]}),
            "ztf_green_mag": round(green[d][1], 5) if d in green else None,
            "ztf_red_mag": round(pre["r"][d][1], 5) if d in pre["r"] else None,
        }
        for d in sorted(flux)
        if START <= d <= END
    ]
    return {
        "cone": cone,
        "raw_measurements": len(raw),
        "good_nights_2019_to_anchor": len(flux),
        "outside_interval_reference": {
            "bright_ztf_quartile_dates": len(bright),
            "faint_ztf_quartile_dates": len(faint),
            "median_asassn_flux_mjy": basis,
        },
        "all_overlapping_dates": len(overlap),
        "all_rank": rank(overlap),
        "interval_overlapping_dates": len(common),
        "interval_rank": rank(common),
        "interval_nights": within,
        "at_or_below_faint_state_median": [
            r["date"] for r in within if r["flux_mjy"] <= basis["faint"]
        ],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--fetch", action="store_true")
    g.add_argument("--import-from", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "asas-results.json")
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch, a.import_from))
