"""Replay the existing SIMBAD, variable-catalog and BTS cones on batch 5."""

import argparse
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
catalog = runpy.run_path(str(HERE.parents[1] / 'batch-4/catalog/check.py'))['catalog']
catalog.HERE = catalog.catalog.HERE = HERE
archive = catalog.catalog

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fetch', action='store_true')
    print(json.dumps(catalog.run(p.parse_args().fetch), indent=2))
