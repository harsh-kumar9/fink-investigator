import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import compass
import compass_archive


class RouteTests(unittest.TestCase):
    def test_history_cutoff_unique_nights_and_missing_first_detection(self):
        def row(day, band=1, positive="t"):
            return {
                "i:jd": 2459500.5 + day,
                "i:fid": band,
                "i:isdiffpos": positive,
                "i:drb": 0.9,
                "i:sigmapsf": 0.05,
            }

        a = {**row(3), "i:jdstarthist": 2459500.5}
        history = [row(0), row(0), row(1), row(3), row(4, 2)]
        self.assertEqual(compass.measure(a, history)["eligible"], {})
        history.append(row(2, 2))
        self.assertEqual(
            compass.measure(a, history)["eligible"]["onset"]["score"], [0, 4]
        )
        self.assertNotIn(
            "onset", compass.measure({**a, "i:jdstarthist": None}, history)["eligible"]
        )

    def test_next_batch_preserves_ranks_and_excludes_both_packets(self):
        base = Path(__file__).resolve().parents[1] / "experiments/compass-2025"
        anchors, histories, cover = compass_archive.load(base / "snapshot")
        ledger = json.loads((base / "results/ledger.json").read_text())
        old = [
            (r["id"], r["ra"], r["dec"])
            for r in ledger["objects"]
            if r["prior_overlap"]
        ]
        reviewed = compass.already_reviewed(
            [base / "results/packet.json", base / "astrometry/packet.json"], anchors
        )
        self.assertEqual(len(reviewed), 4)
        rows = compass.select(
            anchors,
            histories,
            old + reviewed,
            ledger["summary"]["route_passed"],
            cover["unavailable_history_ids"],
        )
        self.assertEqual(
            {r["id"] for r in rows if r["selected_for"]},
            {"ZTF25aacfjde", "ZTF18abmfrjj", "ZTF18aagstwh"},
        )
        original = {r["id"]: r for r in ledger["objects"]}
        for r in rows:
            for route, d in r["route_decisions"].items():
                self.assertEqual(
                    d["rank"], original[r["id"]]["route_decisions"][route]["rank"]
                )
            if not r["history_available"]:
                self.assertFalse(r["eligible"])

    def test_unknown_reviewed_ids_do_not_silently_change_sample(self):
        with TemporaryDirectory() as tmp:
            p = Path(tmp) / "packet.json"
            p.write_text('{"cases":[{"id":"ZTF18abcdefg"}]}')
            with self.assertRaisesRegex(ValueError, "different sample"):
                compass.already_reviewed([p], [])


if __name__ == "__main__":
    unittest.main()
