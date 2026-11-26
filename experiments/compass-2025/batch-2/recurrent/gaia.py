"""Compare fixed ZTF waveforms with earlier Gaia DR3 epoch measurements."""

from collections import defaultdict
from datetime import date, datetime, timezone
import argparse
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

from check import BATCH, HERE, NAME, compass_direct, design, fit, nightly, pilot

SOURCE = "4239802740358274688"
QUERY = "SELECT * FROM gaiadr3.epoch_photometry WHERE source_id=" + SOURCE
URL = "https://gaia.aip.de/tap/sync?" + urlencode(
    {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": QUERY}
)
SOURCE_PATH = HERE / "sources/gaia.xml.gz"
RECEIPT = HERE / "sources/gaia.json"
EPOCH = 2455197.5


def gaia_curves(archive):
    root = ElementTree.fromstring(archive)
    info = {i.get("name"): i.get("value") for i in root.findall(".//{*}INFO")}
    rows = root.findall(".//{*}TR")
    if info.get("QUERY_STATUS") != "OK" or info.get("QUERY") != QUERY or len(rows) != 1:
        raise ValueError("Unexpected Gaia result")
    columns = {
        f.get("name"): v.text.split()
        for f, v in zip(root.findall(".//{*}FIELD"), rows[0], strict=True)
    }
    if columns["source_id"] != [SOURCE]:
        raise ValueError("Wrong Gaia source")
    count = int(columns["n_transits"][0])
    curves = {}
    for key in ("g", "bp", "rp", "g_average"):
        band = "g" if key == "g_average" else key
        tag = "g_transit" if band == "g" else band
        time = tag + "_time" if band == "g" else band + "_obs_time"
        reject = [f"variability_flag_{band}_reject"]
        if key == "g":
            reject.append("photometry_flag_noisy_data")
            reject.extend(f"photometry_flag_af{i}_reject" for i in range(1, 10))
        elif band != "g":
            reject.extend(
                [
                    f"photometry_flag_{band}_reject",
                    f"photometry_flag_{band}_unavailable",
                ]
            )
        fields = [time, tag + "_mag", tag + "_flux_over_error", *reject]
        if any(len(columns[field]) != count for field in fields):
            raise ValueError("Misaligned Gaia epoch arrays")
        bn = defaultdict(list)
        for i in range(count):
            if key == "g_average" and int(columns["g_transit_n_obs"][i]) < 5:
                continue
            flags = [columns[field][i].lower() for field in reject]
            if any(f not in ("true", "false") for f in flags):
                raise ValueError("Unexpected Gaia quality flag")
            if "true" in flags:
                continue
            t, mag, snr = [float(columns[field][i]) for field in fields[:3]]
            if math.isfinite(t) and math.isfinite(mag) and snr > 5:
                bn[math.floor(t / 10)].append((t + EPOCH, mag))
        curves[key] = [
            (median(t for t, m in pts), median(m for t, m in pts), len(pts))
            for k, pts in sorted(bn.items())
        ]
    return curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--out", type=Path, default=HERE / "gaia-results.json")
    args = parser.parse_args()
    if args.fetch and not RECEIPT.exists():
        archive = pilot.curl(URL, [])
        gaia_curves(archive)
        SOURCE_PATH.write_bytes(gzip.compress(archive, mtime=0))
        pilot.json_write(
            RECEIPT,
            {
                "url": URL,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(archive).hexdigest(),
            },
        )
    receipt = json.loads(RECEIPT.read_text())
    archive = gzip.decompress(SOURCE_PATH.read_bytes())
    if (
        receipt["url"] != URL
        or receipt["sha256_uncompressed"] != hashlib.sha256(archive).hexdigest()
    ):
        raise ValueError("Gaia source receipt mismatch")
    curves = gaia_curves(archive)
    case = next(
        r
        for r in json.loads((BATCH / "results/ledger.json").read_text())["objects"]
        if r["id"] == NAME
    )
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    timing = json.loads((HERE / "results.json").read_text())["timing"]
    cutoff = date.fromisoformat(timing["training_last"])
    red = nightly(rows, case["jd"] + 1e-7)["r"]
    t, y = np.array([v for d, v in red.items() if d <= cutoff]).T
    origin = t[0]
    scores = {}
    for name, model in timing["models"].items():
        f = 1 / model["period_days"]
        _, coeff = fit(t - origin, y, f)
        scores[name] = {}
        for b, samples in curves.items():
            if len(samples) < 3:
                scores[name][b] = {
                    "bins": len(samples),
                    "spearman": None,
                    "passes": False,
                }
                continue
            tt, gg, _ = np.array(samples).T
            yy = design(tt - origin, f) @ coeff
            rho = float(spearmanr(gg, yy).statistic)
            item = {
                "bins": len(gg),
                "spearman": round(rho, 5),
                "passes": bool(len(gg) >= 10 and rho >= 0.5),
            }
            if b in ("g", "g_average"):
                loo = [
                    float(spearmanr(np.delete(gg, i), np.delete(yy, i)).statistic)
                    for i in range(len(gg))
                ]
                item["leave_one_bin_out_range"] = list(
                    np.round((min(loo), max(loo)), 5)
                )
            scores[name][b] = item
    pilot.json_write(
        args.out,
        {"source_id": SOURCE, "bins": curves, "scores": scores},
    )


if __name__ == "__main__":
    main()
