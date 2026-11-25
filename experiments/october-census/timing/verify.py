"""Recompute the timing dossier offline and check the archived conclusions."""

import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CHECKS = [
    ("analyze.py", "results.json"),
    ("catalog/check.py", "catalog/results.json"),
    ("green/check.py", "green/results.json"),
    ("catalog/ps1.py", "catalog/ps1-results.json"),
    ("alternatives.py", "alternatives.json"),
    ("epochs.py", "epochs.json"),
    ("space.py", "tess.json"),
    ("dips.py", "dips.json"),
    ("dip-check/compare.py", "dip-check/results.json"),
    ("centroids.py", "centroids.json"),
    ("space.py --windows", "tess-windows.json"),
]


def check(a, b, path):
    if isinstance(a, dict) and isinstance(b, dict) and a.keys() == b.keys():
        for k in a:
            check(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            check(x, y, f"{path}[{i}]")
    elif type(a) is float and type(b) is float:
        if not math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-7):
            raise ValueError(f"{path}: {a} != {b}")
    elif type(a) is not type(b) or a != b:
        raise ValueError(f"Mismatch at {path}")


if __name__ == "__main__":
    for command, result in CHECKS:
        output = subprocess.check_output([sys.executable, *command.split()], cwd=ROOT)
        check(json.loads(output), json.loads((ROOT / result).read_text()), result)
        print(f"Verified {result}")
