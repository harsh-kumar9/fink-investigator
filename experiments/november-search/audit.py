"""Replay the registered November sample and independently audit broker failures."""

import argparse
from collections import Counter
import csv
from datetime import date
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import pilot
import vsx


def identifiers(row):
    return (
        {("Gaia", x["source_id"]) for x in row["gaia_matches"]}
        | {("SIMBAD", x["name"]) for x in row["simbad_neighbors"]}
        | {("VSX", x["oid"]) for x in row.get("vsx_matches", [])}
    )


def past_positions():
    prior = []
    old = json.loads((ROOT / "results/analysis.json").read_text())
    for r in old["objects"]:
        prior.append(("2024-10", r["id"], r["ra"], r["dec"], identifiers(r)))
    path = ROOT / "experiments/october-census/results/outcomes.csv"
    with path.open() as f:
        for r in csv.DictReader(f):
            ids = set()
            for field, label in (
                ("motion_source_id", "Gaia"),
                ("gaia_extended_source", "Gaia"),
                ("vsx_oid", "VSX"),
                ("simbad_name", "SIMBAD"),
            ):
                ids.update((label, i) for i in r[field].split("|") if i)
            prior.append(
                (
                    "2026-10",
                    r["id"],
                    float(r["anchor_ra_deg"]),
                    float(r["anchor_dec_deg"]),
                    ids,
                )
            )
    return prior


def overlaps(r, others):
    result = []
    for cohort, oid, ra, dec, ids in others:
        reasons = []
        sep = math.hypot(*pilot.offsets(r["ra"], r["dec"], ra, dec))
        if oid == r["id"]:
            reasons.append("exact object ID")
        if sep <= 5:
            reasons.append("within 5 arcsec")
        reasons.extend(f"{cat} {key}" for cat, key in sorted(identifiers(r) & ids))
        if reasons:
            result.append(
                {
                    "cohort": cohort,
                    "id": oid,
                    "arcsec": round(sep, 3),
                    "reasons": reasons,
                }
            )
    return result


def audit(fetch):
    snapshot = HERE / "snapshot"
    result = pilot.reconcile_snapshot(snapshot)
    rows = result["objects"]
    failures = [r for r in rows if r["status"] == "broker crossmatch failed"]
    # Use the existing October audit exactly: VSX after a broker failure and empty SIMBAD cone.
    vsx.annotate_cases(
        HERE / "failure-vsx",
        [r for r in failures if not r["simbad_neighbors"]],
        (lambda url: pilot.curl(url, [])) if fetch else None,
    )
    prior = past_positions()
    previous_in_sample = []
    review, audits, details = [], [], []
    for r in rows:
        past = overlaps(r, prior)
        within = overlaps(r, previous_in_sample)
        previous_in_sample.append(
            ("2026-11", r["id"], r["ra"], r["dec"], identifiers(r))
        )
        quality = (
            r["drb"] is not None
            and r["drb"] >= 0.8
            and r["usable_nights"] >= 3
            and r["positive_nights"] >= 2
        )
        current = {
            "id": r["id"],
            "primary": r["status"],
            "simbad": r["simbad_neighbors"],
            "vsx": r.get("vsx_matches"),
            "meets_photometric_cuts": quality,
            "tns": r["tns"],
            "previous_cohort_matches": past,
            "within_sample_matches": within,
        }
        details.append(current)
        if (
            not result["gate_passed"]
            or past
            or within
            or "adjudication" not in r
            or vsx.closed_alias(r)
        ):
            continue
        if r["status"] == "review":
            review.append(r)
        elif quality and not r["tns"]:
            audits.append(r)
    _, raw = pilot.snapshot(snapshot)
    start, end = (
        pilot.utc_jd(date.fromisoformat(result[k])) for k in ("start", "stop_exclusive")
    )
    archive = raw["anomaly.json.gz"]
    in_window = [r for r in archive if start <= r["i:jd"] < end]
    dates = sorted({str(pilot.day_of(r["i:jd"])) for r in in_window})
    counts = dict(Counter(r["status"] for r in rows))
    out = {
        "utc_start": result["start"],
        "utc_stop_exclusive": result["stop_exclusive"],
        "fetched_anomaly_packets": len(archive),
        "in_window_anomaly_packets": result["anomaly_alerts"],
        "archive_broker_labels": dict(Counter(r["d:cdsxmatch"] for r in in_window)),
        "utc_dates": dates,
        "distinct_archive_ids": result["distinct_anomaly_objects"],
        "sample_size": len(rows),
        "first_pass": counts,
        "control_agreements": sum(c["pass"] for c in result["calibration"]),
        "catalog_closed_primary": sum(
            r["status"] == "review" and vsx.closed_alias(r) for r in rows
        ),
        "eligible_primary": [r["id"] for r in review],
        "assigned_primary": [r["id"] for r in review[:3]],
        "audit_quality_survivors": [r["id"] for r in audits],
        "audit_simbad_neighbors": sum(bool(r["simbad_neighbors"]) for r in failures),
        "audit_unique_vsx_within_1_arcsec": sum(
            vsx.closed_alias(r) for r in failures if "adjudication" in r
        ),
        "objects": details,
    }
    pilot.json_write(HERE / "results/audit.json", out)
    for name, cases in (
        ("new_source_packet", review[:3]),
        ("outage_followup_packet", audits[:3]),
    ):
        pilot.json_write(
            HERE / f"results/{name}.json",
            {
                "gate_passed": result["gate_passed"],
                "cases": [pilot.case_packet(r) for r in cases],
            },
        )
    print(json.dumps({k: v for k, v in out.items() if k != "objects"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    audit(parser.parse_args().fetch)
