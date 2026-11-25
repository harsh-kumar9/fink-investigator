import json
import shutil
import motion
from datetime import date
import gzip
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pilot as p


class ScreeningTests(unittest.TestCase):
    def test_deduplicated_temporal_selection_excludes_stop(self):
        start = date(2024, 10, 1)
        stop = date(2024, 10, 3)

        def alert(name, delta, score):
            return {
                "i:objectId": name,
                "i:jd": p.utc_jd(start) + delta,
                "d:anomaly_score": score,
            }

        a, b, c, d = [f"ZTF20aaaaaa{letter}" for letter in "abcd"]
        rows = [
            alert(a, 0, -0.2),
            alert(a, 0.1, -0.3),
            alert(b, 1, -0.4),
            alert(c, 1.5, -0.1),
            alert(d, 2, -0.6),
        ]
        sampled = p.select_alerts(rows, 2, start, stop)
        self.assertEqual([r["i:objectId"] for r in sampled], [a, c])
        self.assertEqual(sampled[0]["d:anomaly_score"], -0.3)

    def test_motion_track_with_ra_wrap_and_pmra_cos_dec(self):
        # 0.6 arcsec/yr east for ten years. Source starts just west of RA=0.
        star = {
            "source_id": "1",
            "ra": "359.999",
            "dec": "60",
            "ref_epoch": "2016",
            "pmra": "600",
            "pmdec": "0",
            "pmra_error": "1",
            "pmdec_error": "1",
        }
        anchor = {
            "i:jd": 2451545.0 + 26 * 365.25,
            "i:ra": (359.999 + 5 / (3600 * 0.5)) % 360,
            "i:dec": 60.0001,
        }
        self.assertAlmostEqual(
            p.stellar_matches(anchor, {"1": star})[0]["distance_to_track_arcsec"],
            0.36,
            places=2,
        )
        anchor["i:dec"] += 3 / 3600
        self.assertEqual(p.stellar_matches(anchor, {"1": star}), [])

    def test_low_significance_and_slow_pm_excluded(self):
        star = {
            "source_id": "1",
            "ra": "10",
            "dec": "0",
            "ref_epoch": "2016",
            "pmra": "800",
            "pmdec": "0",
            "pmra_error": "200",
            "pmdec_error": "1",
        }
        anchor = {"i:jd": 2451545.0 + 26 * 365.25, "i:ra": 10, "i:dec": 0}
        self.assertFalse(p.stellar_matches(anchor, {"1": star}))
        self.assertFalse(
            p.stellar_matches(anchor, {"1": {**star, "pmra_error": "1", "pmra": "149"}})
        )

    def test_metadata_only_history_records_have_no_detection(self):
        row = {
            "i:objectId": "nan",
            "i:candid": None,
            "i:jd": None,
            "d:cdsxmatch": "EclBin",
        }
        self.assertEqual(p.histories([row]), {})
        with self.assertRaisesRegex(ValueError, "unexpected Fink object ID"):
            p.histories([{**row, "i:jd": 2459945.5}])

    def test_checksum_rejects_mutated_source(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "x.gz").write_bytes(gzip.compress(b"[]"))
            (folder / "manifest.json").write_text(
                json.dumps({"sources": {"x.gz": {"sha256_uncompressed": "0" * 64}}})
            )
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                p.snapshot(folder)

    def test_reference_date_check_from_archived_responses(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp) / "moving-star"
            checked = p.ROOT / "results/moving-star"
            shutil.copytree(checked, folder)
            self.assertEqual(
                motion.investigate(
                    p.ROOT / "data/pilot", "ZTF18abgobyb", folder, False
                ),
                "inside",
            )
            self.assertEqual(
                (folder / "report.md").read_bytes(),
                (checked / "report.md").read_bytes(),
            )

    def test_replay_assertions_and_gate(self):
        with TemporaryDirectory() as tmp:
            data = p.reconcile_snapshot(p.ROOT / "data/pilot")
            self.assertEqual(len(data["objects"]), 20)
            self.assertEqual(len(data["calibration"]), 7)
            self.assertTrue(data["gate_passed"])
            initial = [r for r in data["objects"] if r["status"] == "review"]
            self.assertEqual(
                {r["vsx_matches"][0]["oid"] for r in initial},
                {"2214357", "2222911", "2223122"},
            )
            self.assertTrue(
                all(r["adjudication"].startswith("closed") for r in initial)
            )
            ids = {r["id"] for r in data["calibration"]}
            for r in data["objects"]:
                self.assertNotIn(r["id"], ids)
            packet = p.outputs(data, Path(tmp))
            self.assertGreater(len(packet["cases"]), 0)
            self.assertEqual(
                packet, json.loads((p.ROOT / "results/review_packet.json").read_text())
            )
            data["gate_passed"] = False
            self.assertEqual(p.outputs(data, Path(tmp))["cases"], [])


if __name__ == "__main__":
    unittest.main()
