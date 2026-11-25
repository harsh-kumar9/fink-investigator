import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest


DATA = Path(__file__).resolve().parents[1] / "experiments/october-census/infrared"


class InfraredComparisonTest(unittest.TestCase):
    def test_frozen_irsa_sources_and_visit_medians(self):
        for record in json.loads((DATA / "sources.json").read_text())["sources"]:
            self.assertEqual(hashlib.sha256((DATA / record["path"]).read_bytes()).hexdigest(),
                             record["sha256"])
        replay = subprocess.check_output([sys.executable, str(DATA / "compare.py")])
        self.assertEqual((DATA / "results.json").read_bytes(), replay)
