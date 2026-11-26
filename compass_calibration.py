"""Replay the moving-star gate and the independently labeled compass controls."""

import argparse
from datetime import date
import json

import pilot

HERE = pilot.ROOT / "experiments/compass-2025/calibration"
ROUTES = ("onset", "recurrent", "fading")
COLUMNS = (
    ",".join(k for k in pilot.COLUMNS.split(",") if k.startswith("i:"))
    + ",i:jdstarthist"
)


def calibration_inputs(folder=HERE, fetch=False):
    labels = json.loads((folder / "expectations.json").read_text())
    names = sorted(r["id"] for r in labels)
    if len(set(names)) != len(names):
        raise ValueError("Duplicate calibration ID")
    params = {"objectId": ",".join(names), "columns": COLUMNS, "output-format": "json"}
    path = folder / "manifest.json"
    manifest = json.loads(path.read_text()) if path.exists() else {"sources": {}}
    filename = "histories.json.gz"
    if not fetch and filename not in manifest["sources"]:
        raise ValueError("Missing frozen calibration histories")
    rows = pilot.store_raw(folder, manifest, filename, pilot.FINK + "/objects", params)
    history = pilot.histories(rows)
    if set(history) - set(names):
        raise ValueError("Unexpected object ID in control response")
    return labels, history


def within_bounds(control, history):
    start = pilot.utc_jd(date.fromisoformat(control["start"]))
    stop = pilot.utc_jd(date.fromisoformat(control["stop_exclusive"]))
    return [r for r in history if start <= r["i:jd"] < stop]


def select_anchor(control, history):
    from compass import quality

    eligible = (r for r in within_bounds(control, history) if quality(r))
    return max(eligible, key=lambda r: (r["i:jd"], -r["i:candid"]), default=None)


def original_controls():
    snapshot = pilot.ROOT / "data/pilot"
    calibrated = pilot.evaluate(snapshot)
    _, raw = pilot.snapshot(snapshot)
    history = pilot.histories(
        [r for f, rows in raw.items() if f.startswith("history-") for r in rows]
    )
    counts = {"tp": 0, "tn": 0, "fp": 0, "fn": 0}
    positions = []
    cases = []
    for c in calibrated["calibration"]:
        name = c["id"]
        anchor = next(r for r in history[name] if str(r["i:candid"]) == c["candid"])
        observed = bool(c["matches"])
        expected = c["expected"]
        counts[("t" if expected == observed else "f") + ("p" if observed else "n")] += 1
        positions.append({"id": name, "ra": anchor["i:ra"], "dec": anchor["i:dec"]})
        cases.append(
            {
                "id": name,
                "name": c["name"],
                "candid": c["candid"],
                "matched": observed,
                "expected": expected,
                "correct_identity_or_rejection": c["pass"],
            }
        )
    positions.append({"id": "ZTF23aacafew", "ra": None, "dec": None})
    return {"passed": calibrated["gate_passed"], **counts, "cases": cases}, positions


def control_positions():
    _, positions = original_controls()
    for folder in (HERE, HERE / "replacement"):
        labels, history = calibration_inputs(folder)
        for c in labels:
            rows = history.get(c["id"], [])
            anchor = select_anchor(c, rows)
            # A position outside the test interval is still relevant to the alias check.
            location = anchor or (rows[0] if rows else None)
            positions.append(
                {
                    "id": c["id"],
                    "ra": location["i:ra"] if location else None,
                    "dec": location["i:dec"] if location else None,
                }
            )
    return positions


def score_panel(labels, history):
    from compass import measure

    counts = {
        route: {"tp": 0, "tn": 0, "fp": 0, "fn": 0, "unevaluable": 0}
        for route in ROUTES
    }
    cases = []
    for c in labels:
        rows = history.get(c["id"], [])
        anchor = select_anchor(c, rows)
        observed = {}
        if anchor:
            # Only numeric measurements, signs and IDs enter the routing function.
            bounded = [r for r in rows if r["i:jd"] <= anchor["i:jd"]]
            observed = measure(anchor, bounded)["eligible"]
        for route, expected in c["expected"].items():
            if anchor is None:
                counts[route]["unevaluable"] += 1
                continue
            selected = route in observed
            key = ("t" if expected == selected else "f") + ("p" if selected else "n")
            counts[route][key] += 1
        cases.append(
            {
                "id": c["id"],
                "name": c["name"],
                "expected": c["expected"],
                "rows_in_anchor_interval": len(within_bounds(c, rows)),
                "candid": str(anchor["i:candid"]) if anchor else None,
                "date": str(pilot.day_of(anchor["i:jd"])) if anchor else None,
                "first_detection_jd": anchor.get("i:jdstarthist") if anchor else None,
                "history_rows_through_anchor": len(bounded) if anchor else 0,
                "observed_eligible": observed,
            }
        )
    for row in counts.values():
        row["passed"] = bool(
            row["tp"] and row["tn"] and not (row["fp"] + row["fn"] + row["unevaluable"])
        )
    return {"cases": cases, "routes": counts}


def run():
    original, _ = original_controls()
    labels, history = calibration_inputs()
    v1 = score_panel(labels, history)
    missing = "ZTF18abjedws"
    if any(
        r["i:jd"] < pilot.utc_jd(date(2025, 1, 1)) for r in history.get(missing, [])
    ):
        raise ValueError("The V2 amendment applies only to absence of pre-2025 data")
    pool, extras = calibration_inputs(HERE / "replacement")
    coverage = {c["id"]: len(within_bounds(c, extras.get(c["id"], []))) for c in pool}
    chosen = next((c for c in pool if coverage[c["id"]]), None)
    if chosen:
        panel = [chosen if c["id"] == missing else c for c in labels]
        v2 = score_panel(panel, {**history, chosen["id"]: extras[chosen["id"]]})
    else:
        v2 = v1
    return {
        "original": original,
        "v1": v1,
        "effective_panel": "V2",
        "coverage_amendment": {
            "unobservable_v1_id": missing,
            "replacement_pool_interval_counts": coverage,
            "chosen": chosen["id"] if chosen else None,
        },
        **v2,
        "passed_routes": [
            k for k in ROUTES if original["passed"] and v2["routes"][k]["passed"]
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fetch", action="store_true", help="Fill a missing historical source response"
    )
    args = parser.parse_args()
    if args.fetch:
        with pilot.locked_directory(HERE):
            calibration_inputs(fetch=True)
            calibration_inputs(HERE / "replacement", fetch=True)
    print(json.dumps(run(), indent=2, allow_nan=False))
