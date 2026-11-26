"""Replay the frozen queue and matching direct photometry for a later review packet."""

import json
from pathlib import Path
import subprocess
import sys

YEAR = Path(__file__).resolve().parent
ROOT = YEAR.parents[1]
sys.path.insert(0, str(ROOT))
from compass import check_replay


def execute(script, *args):
    return subprocess.check_output([sys.executable, script, *args], cwd=ROOT)


def compare(data, expected):
    check_replay(
        json.loads(data),
        json.loads(expected.read_text()),
        str(expected.relative_to(ROOT)),
    )


def queue_and_direct(batch, out, previous):
    packets = [YEAR / p for p in previous]
    execute(
        ROOT / "compass.py",
        "--snapshot",
        YEAR / "snapshot",
        "--out",
        out,
        *[value for p in packets for value in ("--reviewed-packet", p)]
    )
    for name in ("packet.json", "ledger.json", "calibration.json"):
        compare((out / name).read_bytes(), batch / "results" / name)
    if (out / "report.md").read_bytes() != (batch / "results/report.md").read_bytes():
        raise ValueError("Queue report changed")
    result = out / "replay.json"
    execute(
        ROOT / "compass_direct.py",
        "--snapshot",
        YEAR / "snapshot",
        "--ledger",
        out / "ledger.json",
        "--sources",
        batch / "direct",
        "--out",
        result,
    )
    compare(result.read_bytes(), batch / "results/direct.json")
    execute(YEAR / "batch-2/images/replay.py", "--batch", batch)
