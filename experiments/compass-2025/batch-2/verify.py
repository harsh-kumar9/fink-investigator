"""Replay the second 2025 queue and compare each measurement to its archived result."""

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


def compare(actual, path):
    check_replay(actual, json.loads(path.read_text()), str(path.relative_to(HERE)))


def verify(science=False):
    with TemporaryDirectory() as tmp:
        out = Path(tmp)
        execute(
            ROOT / "compass.py",
            "--snapshot",
            YEAR / "snapshot",
            "--out",
            out,
            "--reviewed-packet",
            YEAR / "results/packet.json",
            "--reviewed-packet",
            YEAR / "astrometry/packet.json",
        )
        for name in ("ledger.json", "calibration.json", "packet.json"):
            compare(json.loads((out / name).read_text()), HERE / "results" / name)
        if (out / "report.md").read_bytes() != (
            HERE / "results/report.md"
        ).read_bytes():
            raise ValueError("Queue changed")
        execute(
            ROOT / "compass_direct.py",
            "--snapshot",
            YEAR / "snapshot",
            "--ledger",
            out / "ledger.json",
            "--sources",
            HERE / "direct",
            "--out",
            out / "direct.json",
        )
        compare(
            json.loads((out / "direct.json").read_text()), HERE / "results/direct.json"
        )
        for script, result in [
            ("onset.py", "results/onset.json"),
            ("broker.py", "results/broker.json"),
            ("catalog/check.py", "catalog/results.json"),
            ("catalog/observer_priority.py", "catalog/observers.json"),
        ]:
            compare(json.loads(execute(HERE / script)), HERE / result)
        for script in ("images/replay.py", "images/rebrightenings.py"):
            execute(HERE / script)
        if science:
            compare(
                json.loads(execute(HERE / "timing.py")), HERE / "results/timing.json"
            )
            for script, args, result in [
                ("check.py", [], "results.json"),
                ("field.py", [], "field-results.json"),
                ("field.py", ["--per-band"], "field-per-band.json"),
                ("color.py", [], "color.json"),
                ("gaia.py", [], "gaia-results.json"),
            ]:
                execute(HERE / "recurrent" / script, *args, "--out", out / result)
                compare(
                    json.loads((out / result).read_text()), HERE / "recurrent" / result
                )
    print(
        "Verified the second queue, direct and forced photometry, field checks, catalogs and images."
        + (
            " Stellar timing, injections and Gaia predictions also reproduced."
            if science
            else ""
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--science",
        action="store_true",
        help="Include NumPy, SciPy and Astropy timing analyses",
    )
    verify(parser.parse_args().science)
