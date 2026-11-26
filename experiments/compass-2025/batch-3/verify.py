"""Replay the third queue and compare the measured evidence with the frozen reports."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
YEAR = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from compass import check_replay


def execute(script, *args):
    return subprocess.check_output([sys.executable, script, *args], cwd=ROOT)


def compare(data, expected):
    check_replay(json.loads(data), json.loads((HERE / expected).read_text()), expected)


def verify(science):
    with TemporaryDirectory() as tmp:
        out = Path(tmp)
        packets = [
            YEAR / p
            for p in (
                "results/packet.json",
                "astrometry/packet.json",
                "batch-2/results/packet.json",
            )
        ]
        execute(
            ROOT / "compass.py",
            "--snapshot",
            YEAR / "snapshot",
            "--out",
            out,
            *[x for p in packets for x in ("--reviewed-packet", p)],
        )
        for name in ("packet.json", "ledger.json", "calibration.json"):
            compare((out / name).read_bytes(), "results/" + name)
        if (out / "report.md").read_bytes() != (
            HERE / "results/report.md"
        ).read_bytes():
            raise ValueError("Queue report changed")
        result = out / "replay.json"
        execute(
            ROOT / "compass_direct.py",
            "--snapshot",
            YEAR / "snapshot",
            "--ledger",
            out / "ledger.json",
            "--sources",
            HERE / "direct",
            "--out",
            result,
        )
        compare(result.read_bytes(), "results/direct.json")
        execute(YEAR / "batch-2/images/replay.py", "--batch", HERE)
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
