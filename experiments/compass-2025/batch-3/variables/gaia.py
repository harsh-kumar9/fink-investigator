"""Rank earlier Gaia DR3 visits using red ZTF fits for IRAS 19245+0727."""

import argparse
from collections import defaultdict
from datetime import date, datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode
from xml.etree import ElementTree

import numpy as np
from scipy.stats import spearmanr

from analyze import HERE, curves, previous, pilot, arrays

NAME = "ZTF18aazutgw"
SOURCE = "4295582339667975168"
QUERY = "SELECT * FROM gaiadr3.epoch_photometry WHERE source_id=" + SOURCE
URL = "https://gaia.aip.de/tap/sync?" + urlencode(
    {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "QUERY": QUERY}
)
DATA = HERE / "sources/gaia.xml.gz"
RECEIPT = HERE / "sources/gaia.json"


def parse(raw):
    root = ElementTree.fromstring(raw)
    info = {el.get("name"): el.get("value") for el in root.findall(".//{*}INFO")}
    if info.get("QUERY_STATUS") != "OK" or info.get("QUERY") != QUERY:
        raise ValueError("Unexpected Gaia query status")
    rows = root.findall(".//{*}TR")
    if not rows:
        return {}
    (row,) = rows
    cols = {
        f.get("name"): v.text.split()
        for f, v in zip(root.findall(".//{*}FIELD"), row, strict=True)
    }
    assert cols["source_id"] == [SOURCE]
    n = int(cols["n_transits"][0])
    out = {}
    for b in ("g", "rp", "bp"):
        tag = "g_transit" if b == "g" else b
        flags = ["variability_flag_" + b + "_reject"]
        if b != "g":
            flags += [
                "photometry_flag_" + b + "_reject",
                "photometry_flag_" + b + "_unavailable",
            ]
        time = "g_transit_time" if b == "g" else b + "_obs_time"
        bins = defaultdict(list)
        for i in range(n):
            if b == "g" and int(cols["g_transit_n_obs"][i]) < 5:
                continue
            if any(cols[k][i].lower() != "false" for k in flags):
                continue
            t, m, snr = [
                float(cols[k][i])
                for k in (time, tag + "_mag", tag + "_flux_over_error")
            ]
            if snr > 5 and math.isfinite(m) and math.isfinite(t):
                bins[math.floor(t / 10)].append((2455197.5 + t, m))
        out[b] = [
            (median(t for t, m in v), median(m for t, m in v), len(v))
            for _, v in sorted(bins.items())
        ]
    return out


def run(fetch):
    if fetch and not RECEIPT.exists():
        raw = pilot.curl(URL, [])
        parse(raw)
        DATA.parent.mkdir(exist_ok=True)
        DATA.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            RECEIPT,
            {
                "url": URL,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            },
        )
    meta = json.loads(RECEIPT.read_text())
    raw = gzip.decompress(DATA.read_bytes())
    assert (
        meta["url"] == URL
        and meta["sha256_uncompressed"] == hashlib.sha256(raw).hexdigest()
    )
    points = parse(raw)
    base = json.loads((HERE / "results.json").read_text())[NAME]
    case, pre = curves(NAME)
    tr = [d for d in pre["r"] if d <= date.fromisoformat(base["training_last"])]
    t, y = arrays(pre["r"], tr)
    scores = {}
    for key, model in base["clocks"].items():
        f = 1 / model["period_days"]
        _, c = previous.fit(t - t[0], y, f)
        score = {}
        for b, v in points.items():
            if len(v) < 3:
                continue
            tt, m, nn = np.array(v).T
            yy = previous.design(tt - t[0], f) @ c
            rho = float(spearmanr(yy, m).statistic)
            score[b] = {
                "transits": int(sum(nn)),
                "bins": len(m),
                "spearman": round(rho, 5),
                "passes": bool(len(m) >= 10 and rho >= 0.5),
            }
        scores[key] = score
    return {"source_id": SOURCE, "gaia_10_day_bins": points, "clock_scores": scores}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "gaia-results.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
