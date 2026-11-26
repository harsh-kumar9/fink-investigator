"""Repeat the checked catalog cones on the third batch, including the BTS registry."""

import argparse
import csv
import importlib.util
import io
import json
from pathlib import Path
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
source_path = HERE.parents[1] / "batch-2/catalog/check.py"
spec = importlib.util.spec_from_file_location("batch_two_catalog", source_path)
catalog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(catalog)
catalog.HERE = HERE
catalog.GAIA_VARIABLES = {
    "ZTF18aazutgw": ("I/358/vlpv", "4295582339667975168", ["Freq", "e_Freq", "Amp"]),
}


def run(fetch=False):
    result = catalog.run(fetch)
    cases = json.loads((HERE.parent / "results/ledger.json").read_text())["objects"]
    (onset,) = [r for r in cases if r["selected_for"] == "onset"]
    params = {
        "f": "s", "startra": onset["ra"] - 0.01, "endra": onset["ra"] + 0.01,
        "startdec": onset["dec"] - 0.01, "enddec": onset["dec"] + 0.01,
        "showdiscdate": "y", "showpeakdate": "y", "format": "csv",
    }
    url = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?" + urlencode(params)
    raw = catalog.source("bts.csv", url, fetch)
    reader = csv.DictReader(io.StringIO(raw.decode().lstrip()))
    if not {"ZTFID", "IAUID", "RA", "Dec", "discdate", "type"}.issubset(reader.fieldnames or []):
        raise ValueError("BTS did not return its source table")
    result[onset["id"]]["bts_in_box"] = list(reader)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    print(json.dumps(run(args.fetch), indent=2))
