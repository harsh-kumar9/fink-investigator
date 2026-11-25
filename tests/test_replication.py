import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pilot
import replication
import vsx


class ReplicationTests(unittest.TestCase):
    def test_frozen_replication_counts_and_replay(self):
        experiment = pilot.ROOT / "experiments/replication/2026-10"
        original = pilot.reconcile_snapshot(pilot.ROOT / "data/pilot")
        original_review = [r for r in original["objects"] if r["status"] == "review"]
        self.assertEqual(sum(vsx.closed_alias(r) for r in original_review), 3)
        with TemporaryDirectory() as tmp:
            dest = Path(tmp)
            packet = pilot.outputs(
                pilot.reconcile_snapshot(experiment / "snapshot"), dest
            )
            replication.summarize(experiment, dest)
            for output in dest.iterdir():
                self.assertEqual(
                    output.read_bytes(),
                    (experiment / "results" / output.name).read_bytes(),
                )
            self.assertTrue(packet["gate_passed"])
            counts = json.loads((dest / "summary.json").read_text())
            self.assertEqual(
                counts["first_pass"],
                {"broker crossmatch failed": 5, "review": 4, "catalog context": 11},
            )
            self.assertEqual(counts["vsx_closed"], 4)
            self.assertEqual(counts["after_vsx"], [])
            self.assertEqual(counts["baseline_object_id_overlap"], [])
            failed = counts["failed_crossmatch_audit"]
            self.assertEqual(len(failed), 5)
            self.assertEqual(sum(bool(r["simbad"]) for r in failed), 3)
            self.assertEqual(sum(bool(r["vsx"]) for r in failed), 2)
