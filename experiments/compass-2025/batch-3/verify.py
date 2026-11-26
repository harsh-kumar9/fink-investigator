"""Replay the third queue and compare the measured evidence with the frozen reports."""

import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
YEAR = HERE.parent
sys.path.insert(0, str(YEAR))
from batch_replay import compare as compare_result, execute, queue_and_direct


def compare(data, expected):
    compare_result(data, HERE / expected)


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
            ),
        )
        result = out / "replay.json"
        for script, expected in (
            ("catalog/check.py", "catalog/results.json"),
            ("catalog/recurrent_prior.py", "catalog/recurrent-prior.json"),
            ("catalog/sources.py", "catalog/sources.json"),
            ("images/ps1/check.py", "images/ps1/results.json"),
            ("images/ps1/field.py", "images/ps1/field.json"),
        ):
            compare(execute(HERE / script), expected)
        outputs = [("timeline.py", "results/timeline.json")]
        if science:
            execute(HERE / "fallback-validation/verify.py")
            compare(
                execute(HERE / "catalog/sources.py", "--verify-parquet"),
                "catalog/sources.json",
            )
            compare(execute(HERE / "catalog/timing.py"), "catalog/timing.json")
            outputs += [("onset.py", "results/apertures.json")]
            for script, filename in (
                ("analyze", "results"),
                ("gaia", "gaia-results"),
                ("interval", "interval"),
                ("eruptions", "eruptions"),
                ("tess", "tess-results"),
                ("asas", "asas-results"),
                ("recurrence", "recurrence"),
                ("controls", "controls"),
                ("mira_survey", "mira-survey"),
            ):
                outputs.append((f"variables/{script}.py", f"variables/{filename}.json"))
        for script, expected in outputs:
            execute(HERE / script, "--out", result)
            compare(result.read_bytes(), expected)
    print(
        "Verified batch 3 queue, catalog crossmatches, images, PS1 controls and ZTF timeline."
        + (
            " Original-image photometry and multi-survey analyses reproduced."
            if science
            else ""
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--science",
        action="store_true",
        help="Also reduce original images and multi-survey light curves",
    )
    verify(p.parse_args().science)
