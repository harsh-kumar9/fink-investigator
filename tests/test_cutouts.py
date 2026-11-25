import hashlib
import json
from pathlib import Path
import unittest

import cutouts
import pilot


class CutoutVerification(unittest.TestCase):
    def test_frozen_stamps_and_epochs(self):
        root = Path(__file__).resolve().parents[1]
        folder = root / "data/pilot"
        _, raw = pilot.snapshot(folder)
        rows = [row for name, pack in raw.items() if name.startswith("history-") for row in pack]
        histories = pilot.histories(rows)
        objects = {r["id"]: r for r in pilot.evaluate(folder)["objects"]}
        for case in ("ZTF18abgobyb", "ZTF17aaagxae"):
            with self.subTest(case=case):
                gallery = root / "review" / case
                manifest = json.loads((gallery / "manifest.json").read_text())
                expected = cutouts.epochs(histories[case], objects[case])
                self.assertEqual([str(r["i:candid"]) for r in expected],
                                 [r["candid"] for r in manifest["epochs"]])
                for epoch in manifest["epochs"]:
                    self.assertEqual([s["kind"] for s in epoch["images"]], ["Science", "Template", "Difference"])
                    for stamp in epoch["images"]:
                        png = (gallery / stamp["file"]).read_bytes()
                        self.assertEqual(hashlib.sha256(png).hexdigest(), stamp["sha256"])
