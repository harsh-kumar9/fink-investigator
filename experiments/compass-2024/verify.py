"""Rebuild the independent-year cohort, queue and direct-photometry checks offline."""

import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
YEAR = HERE.parent / "compass-2025"
sys.path.insert(0, str(YEAR))
from batch_replay import ROOT, compare, execute


def verify(science):
    with TemporaryDirectory() as tmp:
        out = Path(tmp)
        execute(
            ROOT / "compass.py",
            "--snapshot",
            HERE / "snapshot",
            "--out",
            out,
            "--exclude-ledger",
            YEAR / "results/ledger.json",
        )
        for name in ("calibration.json", "ledger.json", "packet.json"):
            compare((out / name).read_bytes(), HERE / "results" / name)
        if (out / "report.md").read_bytes() != (
            HERE / "results/report.md"
        ).read_bytes():
            raise ValueError("2024 queue changed")
        execute(
            ROOT / "compass_direct.py",
            "--snapshot",
            HERE / "snapshot",
            "--ledger",
            out / "ledger.json",
            "--sources",
            HERE / "direct",
            "--out",
            out / "direct.json",
        )
        compare((out / "direct.json").read_bytes(), HERE / "results/direct.json")
        for script, name in [("check", "results"), ("priors", "priors")]:
            compare(
                execute(HERE / f"catalog/{script}.py"), HERE / f"catalog/{name}.json"
            )
        if science:
            execute(HERE / "waveform.py", "--out", out / "waveform.json")
            compare(
                (out / "waveform.json").read_bytes(), HERE / "results/waveform.json"
            )
    print(
        "Verified the independent 2024 archive, gates, ranked sources and direct photometry."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--science", action="store_true")
    verify(p.parse_args().science)
