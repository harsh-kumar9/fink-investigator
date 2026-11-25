import json
from pathlib import Path
import subprocess
import sys
import unittest


class EpochPositionCheck(unittest.TestCase):
    def test_borderline_vsx_match_replays_from_the_catalog_archive(self):
        root = Path(__file__).resolve().parents[1]
        program = root / "experiments/replication/2026-10/astrometry/check.py"
        result = json.loads(subprocess.check_output([sys.executable, str(program)]))
        self.assertEqual(result["gaia_id"], "4499667881166202624")
        self.assertEqual(result["quality_positive_nights"], 269)
        self.assertEqual(result["median_nightly_separation_gaia_arcsec"], 0.134)
        self.assertEqual(result["median_nightly_separation_vsx_arcsec"], 0.896)
