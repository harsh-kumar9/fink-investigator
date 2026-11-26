"""Check the independent AAVSO observer windows registered for FT Cam and V1060 Cyg."""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import date, timedelta
import io
import json
import math
from statistics import median
from urllib.parse import urlencode

from check import HERE, pilot, source
import compass
import compass_archive


def retrieve(star, key, lo, hi, fetch):
    params = {
        "view": "api.delim",
        "ident": star,
        "fromjd": f"{lo:.7f}",
        "tojd": f"{hi:.7f}",
        "delimiter": "|",
    }
    url = "https://vsx.aavso.org/index.php?" + urlencode(params)
    raw = source(key + ".aid", url, fetch)
    reader = csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig")), delimiter="|", quoting=csv.QUOTE_NONE
    )
    if not {
        "JD",
        "mag",
        "uncert",
        "band",
        "by",
        "starName",
        "obsID",
        "fainterThan",
        "obsType",
        "val",
    }.issubset(reader.fieldnames or []):
        raise ValueError(f"Unexpected AAVSO response for {star}: {raw[:200]!r}")
    rows = list(reader)
    if len(rows) >= 10000 or any(
        r["starName"].upper() != star.upper() or None in r for r in rows
    ):
        raise ValueError("Truncated, mismatched or malformed AAVSO table")
    # Some date APIs round their endpoints to a date. Cut precisely here.
    rows = [r for r in rows if lo <= float(r["JD"]) <= hi]
    if len({r["obsID"] for r in rows}) != len(rows):
        raise ValueError("Duplicate AAVSO observation")
    series = defaultdict(lambda: defaultdict(list))
    types = Counter()
    codes = set()
    for r in rows:
        types[(r["obsType"], r["fainterThan"], r["val"])] += 1
        if (
            r["obsType"] not in {"CCD", "DSLR", "PEP"}
            or r["val"] not in {"V", "Z"}
            or r["fainterThan"] != "0"
        ):
            continue
        if not (0 < float(r["uncert"] or 0) <= 0.2 and math.isfinite(float(r["mag"]))):
            continue
        if "ZTF" in (r["credit"] + " " + r["comment"] + " " + r["software"]).upper():
            continue
        key = r["by"] + "/" + r["band"]
        codes.add(r["by"])
        series[key][pilot.day_of(float(r["JD"]))].append(float(r["mag"]))
    collapsed = {
        k: {d: median(v) for d, v in group.items()} for k, group in series.items()
    }
    return collapsed, {
        "received_in_bounds": len(rows),
        "observer_codes_with_usable_data": sorted(codes),
        "counts_by_type_limit_validation": {
            " / ".join(k): v for k, v in sorted(types.items())
        },
    }


def summarize(data, condition):
    vals = [v for d, v in data.items() if condition(d)]
    return {"nights": len(vals), "median_mag": round(median(vals), 4) if vals else None}


def run(fetch):
    root = HERE.parent
    _, histories, _ = compass_archive.load(root / "snapshot")
    ledger = json.loads((root / "results/ledger.json").read_text())
    cases = {r["selected_for"]: r for r in ledger["objects"] if r["selected_for"]}
    results = {}
    ft = cases["recurrent"]
    start1, stop1, start2, stop2 = map(
        date.fromisoformat, ft["eligible"]["recurrent"]["example_dates"]
    )
    for i, (a, b) in enumerate([(start1, stop1), (start2, stop2)], 1):
        lo, hi = (
            pilot.utc_jd(a - timedelta(days=30)),
            pilot.utc_jd(b + timedelta(days=31)) - 1e-6,
        )
        series, info = retrieve("FT Cam", f"FT-Cam-{i}", lo, hi, fetch)
        info.update(
            dates=[str(a), str(b)],
            series={
                key: {
                    "in_window": summarize(data, lambda d: a <= d <= b),
                    "surrounding": summarize(data, lambda d: d < a or d > b),
                }
                for key, data in sorted(series.items())
            },
        )
        results[f"FT Cam window {i}"] = info
    cyg = cases["fading"]
    end = cyg["jd"]
    series, info = retrieve("V1060 Cyg", "V1060-Cyg", end - 365, end, fetch)
    signs = defaultdict(set)
    for r in histories[cyg["id"]]:
        if end - 365 <= r["i:jd"] <= end and r["i:fid"] == 1 and compass.quality(r):
            signs[pilot.day_of(r["i:jd"])].add(compass.sign(r))
    info["series"] = {}
    for key, data in sorted(series.items()):
        bright = summarize(data, lambda d: signs[d] == {1})
        faint = summarize(data, lambda d: signs[d] == {-1})
        info["series"][key] = {
            "positive": bright,
            "negative": faint,
            "negative_minus_positive_mag": (
                round(faint["median_mag"] - bright["median_mag"], 4)
                if min(bright["nights"], faint["nights"]) >= 5
                else None
            ),
        }
    results["V1060 Cyg"] = info
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    print(json.dumps(run(p.parse_args().fetch), indent=2, allow_nan=False))
