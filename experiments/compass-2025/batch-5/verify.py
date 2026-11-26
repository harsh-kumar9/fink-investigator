"""Rebuild the fifth 2025 queue, catalogs and source measurements offline."""

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
        prior = [
            "results/packet.json",
            "astrometry/packet.json",
            *[f"batch-{n}/results/packet.json" for n in (2, 3, 4)],
        ]
        queue_and_direct(HERE, out, prior)
        for stem, name in [
            ("check", "results"),
            ("priors", "priors"),
            ("ps1", "ps1"),
            ("early_state", "early-state"),
        ]:
            compare(execute(HERE / f"catalog/{stem}.py"), HERE / f"catalog/{name}.json")
        jobs = [
            ("light_curve.py", [], "results/light-curve.json"),
            ("type_ii.py", [], "results/type-ii.json"),
        ]
        if science:
            for script in ("images/host/field.py", "images/host/hst/hubble.py"):
                execute(HERE / script)
            jobs += [
                (f"variables/{script}.py", [], f"variables/{expected}.json")
                for script, expected in (
                    ("analyze_vars", "results"),
                    ("duration", "duration"),
                    ("cadence", "cadence"),
                    ("alternating", "alternating"),
                    ("tess_data", "tess"),
                    ("tess_events", "tess-events"),
                    ("checks", "checks"),
                )
            ]
            jobs += [
                ("season.py", ["--band", b], f"results/season-{b}.json") for b in "ri"
            ]
            jobs += [
                ("standards.py", [], "results/ps1-calibrated.json"),
                ("standards.py", ["--wide"], "results/field-calibrated.json"),
            ]
        for script, args, result in jobs:
            dst = out / "actual.json"
            execute(HERE / script, *args, "--out", dst)
            compare(dst.read_bytes(), HERE / result)
    print(
        "Verified the fifth 2025 queue, catalog identities, alert images, forced fluxes and Type II controls."
        + (
            " Original images, Hubble registration, TESS and cross-season tests reproduced."
            if science
            else ""
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--science", action="store_true")
    verify(p.parse_args().science)
