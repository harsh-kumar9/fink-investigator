"""Replay this batch's Fink thumbnails, Gaia neighbors and signed offsets."""

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
YEAR = BATCH.parent
sys.path.insert(0, str(YEAR / "images"))
import check as images


def run(fetch=False):
    # The first batch's capture and validation routines use its archive directory.
    images.HERE = HERE
    packet = json.loads((BATCH / "results/packet.json").read_text())
    if not packet["gate_passed"]:
        raise ValueError("Closed review packet")
    anchors, histories, _ = images.compass_archive.load(YEAR / "snapshot")
    by_id = {images.pilot.oid(a): a for a in anchors}
    results = []
    for case in packet["cases"]:
        name = case["id"]
        anchor = by_id[name]
        record = images.case_summary(anchor, histories[name], fetch)
        stars = images.gaia(anchor, False)
        ordered = sorted(stars, key=lambda s: images.sep(anchor, s, True))
        record["nearest_two_at_anchor"] = [
            {"id": s["source_id"], "separation_arcsec": round(images.sep(anchor, s, True), 3)}
            for s in ordered[:2]
        ]
        results.append(record)
    return results


def plot(results):
    from figure import panels
    import figure

    figure.HERE = HERE
    for case in results:
        folder = HERE / case["id"]
        rows = []
        for receipt in (folder / "epochs/manifest.json", folder / "negative/manifest.json"):
            if not receipt.exists():
                continue
            rows.extend((receipt.parent, e) for e in json.loads(receipt.read_text())["epochs"])
        panels(case["id"], sorted(rows, key=lambda pair: pair[1]["jd"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    results = run(args.fetch)
    expected = HERE / "results.json"
    if args.fetch:
        images.pilot.json_write(expected, results)
    elif json.loads(expected.read_text()) != results:
        raise ValueError("Archived image or Gaia results changed")
    if args.plot:
        plot(results)
    print("Verified", sum(len(c["checked_image_candids"]) for c in results), "cutout triplets across", len(results), "sources")
