"""Run the same coordinate, variable-catalog and BTS queries used in batch 3."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
previous = HERE.parents[1] / "batch-3/catalog/check.py"
spec = importlib.util.spec_from_file_location("batch_three_catalog", previous)
catalog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(catalog)
catalog.HERE = catalog.catalog.HERE = HERE
catalog.catalog.GAIA_VARIABLES = {}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    print(json.dumps(catalog.run(p.parse_args().fetch), indent=2))
