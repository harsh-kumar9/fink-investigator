"""Compare equal windows on opposite halves of the frozen 11.62-day cycle."""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

TIMING = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TIMING))
from analyze import CENSUS, nights
from green import check as green


def medians(data, fit):
    phase = (data[:, 1] - fit["center_hjd"]) / fit["period_days"]
    half_width = fit["duration_days"] / (2 * fit["period_days"])
    deep = np.abs(phase - np.round(phase)) < half_width
    offset = phase - 0.5
    intervening = np.abs(offset - np.round(offset)) < half_width
    result = {}
    for key, mask in (("deep",deep), ("opposite",intervening), ("outside_both", ~(deep | intervening))):
        values = data[mask, 2]
        result[key] = {"n":len(values), "median_mag":round(float(np.median(values)),4) if len(values) else None}
    return result


def atlas_fields():
    folder = TIMING / "catalog"
    raw = gzip.decompress((folder / "cone.xml.gz").read_bytes())
    digest = json.loads((folder / "source.json").read_text())["sha256_uncompressed"]
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError("ATLAS cone checksum mismatch")
    table = ET.fromstring(raw).find(".//{*}TABLE[@name='J/AJ/156/241/table4']")
    fields = [f.attrib["name"] for f in table.findall("{*}FIELD")]
    tr, = table.findall(".//{*}TR")
    row = dict(zip(fields, [(td.text or "").strip() for td in tr], strict=True))
    selected = ("Class","fp-LSper","fp-lngfitper","fp-domper-c","fp-domper-o","fp-fit",
                "fp-period","fp-domperiod","ls-Pday","P(CBF)","P(CBH)","P(DBF)","P(DBH)","P(dubious)")
    return {k:row[k] for k in selected}


def run():
    dips = json.loads((TIMING / "dips.json").read_text())
    fit = dips["fits"][dips["training_winner"]]
    selected = json.loads((CENSUS / "results/direct.json").read_text())["exploration"]
    case, = [r for r in selected if r["id"] == green.NAME]
    with gzip.open(CENSUS / "sources/direct" / (green.NAME + ".csv.gz"), "rt") as f:
        ids = {r["oid"] for r in csv.DictReader(f) if r["filtercode"] == "zr"}
    red = {oid:nights({**case, "direct_source_id":oid}) for oid in sorted(ids)}
    red_dates = np.concatenate([np.floor(v[:,0]) for v in red.values()])
    primary = red[case["direct_source_id"]]
    _, stars = green.select(False)
    green_series = green.green_data(stars, False)
    result = {"period_days":fit["period_days"], "center_hjd":fit["center_hjd"],
              "full_window_days":fit["duration_days"], "atlas_catalog":atlas_fields(), "green":{}, "red":{}}
    for r in stars:
        oid = r["oid"]
        data = green_series[oid]
        only_green = data[~np.isin(np.floor(data[:,0]), red_dates)]
        result["green"][oid] = {"all":medians(data, fit), "absent_from_both_red_fields":medians(only_green, fit)}
    n = len(primary)*2//3
    result["red"][case["direct_source_id"]] = {"all":medians(primary, fit), "later":medians(primary[n:], fit)}
    for oid in sorted(ids - {case["direct_source_id"]}):
        data = red[oid]
        separate = data[~np.isin(np.floor(data[:,0]), np.floor(primary[:,0]))]
        result["red"][oid] = {"all":medians(data, fit), "disjoint_from_primary":medians(separate, fit)}
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
