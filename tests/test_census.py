import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import census
import census_photometry
import census_report


class CensusTests(unittest.TestCase):
    def test_all_registered_ids_and_followups_replay(self):
        with TemporaryDirectory() as tmp:
            dest = Path(tmp)
            counts = census_report.summarize(dest)
            self.assertEqual(counts["archive_packets"], 190)
            self.assertEqual(
                counts["cohort_sizes"], {"initial sample": 20, "remainder": 110}
            )
            with (dest / "outcomes.csv").open() as f:
                rows = {r["id"]: r for r in csv.DictReader(f)}
            self.assertEqual(len(rows), 130)
            self.assertEqual(counts["first_pass"]["review"], 20)
            self.assertEqual(counts["by_first_pass"]["review"]["VSX counterpart"], 14)
            self.assertEqual(
                rows["ZTF26acbcwpz"]["disposition"], "registered type Ia supernova"
            )
            self.assertEqual(rows["ZTF26acbcwpz"]["packet_tns"], "not supplied")
            self.assertEqual(
                rows["ZTF18ablrppd"]["first_pass"], "broker crossmatch failed"
            )
            self.assertEqual(rows["ZTF18aaaduso"]["vsx_oid"], "838077")
            self.assertEqual(
                rows["ZTF22abgqdcj"]["disposition"], "moving-star explanation supported"
            )
            files = list(dest.glob("review_packet_*.json"))
            ids = set()
            for name in files:
                packet = json.loads(name.read_text())
                self.assertTrue(packet["gate_passed"])
                self.assertEqual(len(packet["cases"]), 3)
                ids.update(c["id"] for c in packet["cases"])
            self.assertEqual(ids, set(census_report.IMAGE_READINGS))
            for name in dest.iterdir():
                self.assertEqual(
                    name.read_bytes(),
                    (census.CENSUS / "results" / name.name).read_bytes(),
                )

    def test_direct_control_and_registered_sign_coverage(self):
        data = census_photometry.run(False)
        self.assertTrue(data["calibration_passed"])
        self.assertEqual(data["control"]["matched_dates"], 158)
        target = next(r for r in data["exploration"] if r["id"] == "ZTF17aacydyv")
        self.assertEqual((target["matched_dates"], target["negative_dates"]), (130, 6))
        self.assertEqual(target["spearman"], 0.995)
        for r in data["exploration"]:
            if min(r["positive_dates"], r["negative_dates"]) < 5:
                self.assertFalse(r["threshold_supported"])
                self.assertIsNone(r["negative_minus_positive_median_direct_mag"])


if __name__ == "__main__":
    unittest.main()
