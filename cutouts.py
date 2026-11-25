"""Archive contemporaneous ZTF stamps for a shortlisted case or a flagged moving star."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import pilot


def epochs(history, anchor):
    """Earliest and middle usable prior dates in the alert band, plus the alert."""
    prior = {}
    anchor_day = pilot.day_of(anchor["jd"])
    for row in sorted(history, key=lambda r: r["i:jd"]):
        day = pilot.day_of(row["i:jd"])
        if (day < anchor_day and row["i:fid"] == anchor["filter"]
                and (row["i:drb"] or 0) >= 0.8 and row["i:sigmapsf"] is not None
                and 0 < row["i:sigmapsf"] <= 0.2):
            prior.setdefault(day, row)
    rows = list(prior.values())
    selected = [rows[i] for i in sorted({0, len(rows) // 2})] if rows else []
    selected.append(next(r for r in history if str(r["i:candid"]) == str(anchor["candid"])))
    return selected


def capture(snapshot, case_id, destination):
    analysis = pilot.evaluate(snapshot)
    eligible = {"review", "moving-star alternative"}
    case = next((r for r in analysis["objects"] if r["id"] == case_id and r["status"] in eligible), None)
    if not analysis["gate_passed"] or case is None:
        raise ValueError("case was not admitted for review")
    if destination.exists():
        raise ValueError("choose a new output directory to preserve the previous capture")
    _, raw = pilot.snapshot(snapshot)
    history = pilot.histories([r for key, rows in raw.items() if key.startswith("history-") for r in rows])
    capture_images(case, history[case_id], destination)


def capture_images(case, history, destination):
    case_id = case["id"]
    rows = epochs(history, case)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=destination.parent) as tmp:
        stage = Path(tmp) / "capture"
        stage.mkdir()
        manifest = {"id": case_id, "selection": epochs.__doc__, "epochs": []}
        gallery = [f"# ZTF stamps for {case_id}", "", epochs.__doc__, "",
                   "Dates are UTC. The archive stretches each PNG for display. Pixel brightness between panels is not calibrated.",
                   "", "| Date and candid | Science | Template | Difference |", "|---|---|---|---|"]
        for row in rows:
            day = str(pilot.day_of(row["i:jd"]))
            candid = str(row["i:candid"])
            epoch = {"date_utc": day, "jd": row["i:jd"], "candid": candid, "filter": row["i:fid"],
                     "isdiffpos": row["i:isdiffpos"], "diff_mag": row["i:magpsf"], "images": []}
            cells = [f"{day} · {candid}"]
            for kind in ("Science", "Template", "Difference"):
                url = pilot.image_url({"id": case_id, "candid": candid}, kind)
                png = subprocess.check_output(["curl", "--fail", "--silent", "--show-error",
                                               "--max-time", "30", "--max-filesize", "500000", url])
                if not png.startswith(b"\x89PNG\r\n\x1a\n"):
                    raise ValueError(f"archive returned a non-PNG for {url}")
                name = f"{candid}-{kind.lower()}.png"
                (stage / name).write_bytes(png)
                epoch["images"].append({"kind": kind, "url": url, "file": name,
                                        "sha256": hashlib.sha256(png).hexdigest()})
                cells.append(f"[![{kind} on {day}]({name})]({url})")
            gallery.append("| " + " | ".join(cells) + " |")
            manifest["epochs"].append(epoch)
        (stage / "gallery.md").write_text("\n".join(gallery) + "\n")
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id")
    parser.add_argument("--snapshot", type=Path, default=pilot.ROOT / "data/pilot")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    capture(args.snapshot, args.case_id, args.out)
