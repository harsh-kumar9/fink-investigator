"""Rebuild the fourth 2025 review packet and test the archived scientific measurements."""

import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from batch_replay import compare, execute, queue_and_direct


def verify(science):
    with TemporaryDirectory() as tmp:
        out = Path(tmp)
        queue_and_direct(
            HERE,
            out,
            (
                "results/packet.json",
                "astrometry/packet.json",
                "batch-2/results/packet.json",
                "batch-3/results/packet.json",
            ),
        )
        printed = [
            ("catalog/check.py", "catalog/results.json"),
            ("catalog/priors.py", "catalog/priors.json"),
        ]
        saved = [
            ("onset.py", [], "results/onset.json"),
            ("comparison.py", [], "results/comparison.json"),
            ("sensitivity.py", [], "results/sensitivity.json"),
        ]
        if science:
            for f in ("images/host/check.py", "images/host/legacy.py"):
                execute(HERE / f)
            printed += [
                ("catalog/ps1.py", "catalog/ps1.json"),
                ("catalog/field.py", "catalog/field.json"),
            ]
            saved += [
                ("late.py", [], "results/late.json"),
                ("late.py", ["--field", "633"], "results/late-other.json"),
            ]
            for script, filename in (
                ("measure", "results"),
                ("crosscheck", "crosscheck"),
                ("sky", "sky"),
                ("field", "field"),
                ("pulsation", "pulsation"),
                ("template", "template"),
                ("pixels", "pixels"),
            ):
                saved.append(
                    (f"variables/{script}.py", [], f"variables/{filename}.json")
                )
        for script, expected in printed:
            compare(execute(HERE / script), HERE / expected)
        for script, args, expected in saved:
            output = out / "actual.json"
            execute(HERE / script, *args, "--out", output)
            compare(output.read_bytes(), HERE / expected)
    print(
        "Verified batch 4 queue, source identities, image triplets, forced fluxes and classified comparison samples."
        + (
            " Original-image measurements, host analysis, 2013 eruption and long-history tests reproduced."
            if science
            else ""
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--science", action="store_true")
    verify(p.parse_args().science)
