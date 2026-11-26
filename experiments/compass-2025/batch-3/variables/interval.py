"""Look for intervening faint dates inside MGAB-V243's red-selected 68-day group."""

import argparse
from collections import Counter
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

import numpy as np

from analyze import BATCH, HERE, compass_direct, curves, previous, pilot

NAME = "ZTF17aabshbq"
START, END = date(2020, 11, 10), date(2021, 1, 17)
SRC = HERE / "sources"


def aid(fetch):
    lo, hi = pilot.utc_jd(START) - 30, pilot.utc_jd(END) + 31 - 1e-6
    url = "https://vsx.aavso.org/index.php?" + urlencode(
        {
            "view": "api.delim",
            "ident": "MGAB-V243",
            "fromjd": f"{lo:.7f}",
            "tojd": f"{hi:.7f}",
            "delimiter": "|",
        }
    )
    file, receipt = SRC / "observers.aid.gz", SRC / "observers.json"
    if fetch and not receipt.exists():
        raw = pilot.curl(url, [])
        if not raw.decode("utf-8-sig").lstrip().startswith("JD|"):
            raise ValueError("Unexpected observer archive")
        file.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            },
        )
    raw = gzip.decompress(file.read_bytes())
    meta = json.loads(receipt.read_text())
    assert (
        meta["url"] == url
        and hashlib.sha256(raw).hexdigest() == meta["sha256_uncompressed"]
    )
    dit = csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig")), delimiter="|", quoting=csv.QUOTE_NONE
    )
    rows = [r for r in dit if lo <= float(r["JD"]) <= hi]
    if any(r["starName"].upper() != "MGAB-V243" for r in rows) or len(rows) >= 10000:
        raise ValueError("Observer record mismatch")
    types = Counter((r["obsType"], r["val"], r["fainterThan"]) for r in rows)
    group = []
    for r in rows:
        if (
            r["obsType"] not in ("CCD", "PEP", "DSLR")
            or r["val"] not in ("V", "Z")
            or r["fainterThan"] != "0"
        ):
            continue
        if (
            not 0 < float(r["uncert"] or 0) <= 0.2
            or "ZTF"
            in " ".join(r[k] for k in ("credit", "comment", "software")).upper()
        ):
            continue
        group.append(
            {
                "observer": r["by"],
                "band": r["band"],
                "date": str(pilot.day_of(float(r["JD"]))),
                "mag": float(r["mag"]),
            }
        )
    return {
        "all_records": len(rows),
        "record_types": {" / ".join(k): v for k, v in types.items()},
        "usable": group,
    }


def summarize(curve, cutoff, red):
    tr = [v[1] for d, v in curve.items() if d <= cutoff]
    dates = [d for d in curve if START <= d <= END]
    out = {
        "training_dates": len(tr),
        "interval_dates": len(dates),
        "dates_absent_primary_red": len(set(dates) - red.keys()),
    }
    if len(tr) < 30 or not dates:
        return out
    bright, faint = np.quantile(tr, [0.25, 0.75])
    lows = [d for d in dates if curve[d][1] >= faint]
    vals = [curve[d][1] for d in dates]
    return {
        **out,
        "training_bright_faint_mag": [float(bright), float(faint)],
        "bright_dates": sum(m <= bright for m in vals).item(),
        "range_mag": [min(vals), max(vals)],
        "maximum_gap_days": max(
            (b - a).days for a, b in zip([START, *dates], [*dates, END])
        ),
        "faint": [
            {"date": str(d), "mag": curve[d][1], "new_red_date": d not in red}
            for d in lows
        ],
        "nightly": [[str(d), round(curve[d][1], 5)] for d in dates],
    }


def run(fetch):
    case, pre = curves(NAME)
    cutoff = date.fromisoformat(
        json.loads((HERE / "results.json").read_text())[NAME]["training_last"]
    )
    main, _, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    with gzip.open(BATCH / "direct" / NAME / "cone.csv.gz", "rt") as stream:
        cone = list(csv.DictReader(stream))
    alt = [r for r in cone if r["oid"] not in {v["oid"] for v in main.values()}]
    req = [("ID", r["oid"]) for r in sorted(alt, key=lambda r: int(r["oid"]))] + [
        ("COLLECTION", "ztf_dr24"),
        ("FORMAT", "csv"),
        ("BAD_CATFLAGS_MASK", "65535"),
    ]
    rows = compass_direct.capture(
        SRC,
        "alternate",
        compass_direct.SERVICE + "?" + urlencode(req),
        {"oid", "mag", "filtercode", "expid", "mjd", "hjd", "magerr", "catflags"},
        fetch,
    )
    assert {r["oid"] for r in alt} == {r["oid"] for r in rows}
    out = {}
    for fid, row in main.items():
        b = {1: "g", 2: "r"}[fid]
        out[row["oid"]] = {
            "field": row["field"],
            "band": b,
            **summarize(pre[b], cutoff, pre["r"]),
        }
    for r in alt:
        b = {"1": "g", "2": "r"}[r["fid"]]
        norm = previous.nightly(
            [x for x in rows if x["oid"] == r["oid"]], case["jd"] + 1e-7
        )[b]
        out[r["oid"]] = {
            "field": r["field"],
            "band": b,
            **summarize(norm, cutoff, pre["r"]),
        }
    return {"series": out, "aavso": aid(fetch)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "interval.json")
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch))
