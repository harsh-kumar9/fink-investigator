"""Replay the registered wide search and its source-level checks without a network."""

import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


sys.path.insert(0, str(ROOT))
from compass import check_replay as check


def execute(script, *args):
    return subprocess.check_output([sys.executable, script, *args], cwd=ROOT)


def verify():
    snapshot = HERE / "snapshot"
    with TemporaryDirectory() as tmp:
        out = Path(tmp)
        execute("compass.py", "--snapshot", snapshot, "--out", out)
        for name in ("ledger.json", "calibration.json", "packet.json"):
            check(
                json.loads((out / name).read_text()),
                json.loads((HERE / "results" / name).read_text()),
                name,
            )
        if (out / "report.md").read_text() != (HERE / "results/report.md").read_text():
            raise ValueError("Changed ranked queue")
        execute(
            "compass_direct.py",
            "--snapshot",
            snapshot,
            "--ledger",
            out / "ledger.json",
            "--sources",
            HERE / "direct",
            "--out",
            out / "direct.json",
        )
        check(
            json.loads((out / "direct.json").read_text()),
            json.loads((HERE / "results/direct.json").read_text()),
            "direct",
        )
        execute(
            "compass_astrometry.py",
            "--snapshot",
            snapshot,
            "--ledger",
            out / "ledger.json",
            "--out",
            HERE / "astrometry",
            "--verify",
        )
        execute(
            "compass_reference.py",
            "--snapshot",
            snapshot,
            "--packet",
            HERE / "astrometry/packet.json",
            "--sources",
            HERE / "astrometry/reference-sources",
            "--out",
            out / "reference.json",
        )
        check(
            json.loads((out / "reference.json").read_text()),
            json.loads((HERE / "astrometry/reference.json").read_text()),
            "reference",
        )
    for script, args, result in [
        ("behavior.py", [], "results/behavior.json"),
        ("images/check.py", [], "images/results.json"),
        ("images/centroid_check.py", [], "images/centroid_check.json"),
        ("images/offset_case.py", [], "images/offset-results.json"),
        ("identifications/check.py", [], "identifications/results.json"),
        (
            "identifications/check.py",
            ["--astrometry"],
            "identifications/astrometry.json",
        ),
        ("identifications/observers.py", [], "identifications/observers.json"),
    ]:
        check(
            json.loads(execute(HERE / script, *args)),
            json.loads((HERE / result).read_text()),
            result,
        )
    print(
        "Verified the 2025 archive, calibration, queue, photometry, astrometry, catalogs and case measurements."
    )


if __name__ == "__main__":
    verify()
