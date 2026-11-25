"""Summarize the preregistered cohort from frozen primary and catalog responses."""

import argparse
from collections import Counter
import json
from pathlib import Path
import pilot
import vsx


def summarize(folder, out, fetch_audit=False):
    out.mkdir(parents=True, exist_ok=True)
    rows = json.loads((folder.parent / "coverage/coverage.json").read_text())
    coverage = next(r for r in rows if r["eligible"])
    result = pilot.reconcile_snapshot(folder / "snapshot")
    if (coverage["start"], coverage["anomaly_packets"], coverage["unique_objects"]) != (
        result["start"],
        result["anomaly_alerts"],
        result["distinct_anomaly_objects"],
    ):
        raise ValueError("coverage and analysis snapshots disagree")
    original = {r["id"] for r in pilot.evaluate(pilot.ROOT / "data/pilot")["objects"]}
    _, raw = pilot.snapshot(folder / "snapshot")
    histories = [
        r for key, rows in raw.items() if key.startswith("history-") for r in rows
    ]
    incomplete = sum(pilot.empty_history_record(r) for r in histories)
    objects = result["objects"]
    overlap = sorted(r["id"] for r in objects if r["id"] in original)
    stages = Counter(r["status"] for r in objects)
    eligible = (
        [r for r in objects if r["status"] == "review"] if result["gate_passed"] else []
    )
    survivors = [r for r in eligible if not vsx.closed_alias(r)]
    closed = [r for r in eligible if r not in survivors]
    failed = [dict(r) for r in objects if r["status"] == "broker crossmatch failed"]
    pending = [r for r in failed if not r["simbad_neighbors"]]
    audit_dir = folder / "crossmatch-audit"
    vsx.annotate_cases(
        audit_dir, pending, (lambda url: pilot.curl(url, [])) if fetch_audit else None
    )
    sample_dates = {str(pilot.day_of(r["jd"])) for r in objects}
    metrics = {
        "window": [result["start"], result["stop_exclusive"]],
        "archive_utc_dates": coverage["utc_dates_with_anomalies"],
        "baseline_object_id_overlap": overlap,
        "sample_size": len(objects),
        "sample_utc_dates": len(sample_dates),
        "science_alerts": coverage["science_alerts_in_statistics"],
        "observing_nights": coverage["observing_nights_in_statistics"],
        "alerts_in_anomaly_archive": coverage["anomaly_packets"],
        "distinct_archive_ids": coverage["unique_objects"],
        "history_response_rows": len(histories),
        "metadata_only_history_rows": incomplete,
        "control_agreements": sum(r["pass"] for r in result["calibration"]),
        "control_count": len(result["calibration"]),
        "first_pass": dict(stages),
        "vsx_closed": len(closed),
        "after_vsx": [r["id"] for r in survivors],
    }
    metrics["failed_crossmatch_audit"] = [
        {
            "id": r["id"],
            "broker": r["cdsxmatch"],
            "simbad": r["simbad_neighbors"],
            "vsx": r.get("vsx_matches"),
            "vsx_query": r.get("vsx_query"),
        }
        for r in failed
    ]
    pilot.json_write(out / "summary.json", metrics)
    lines = [
        "# Fresh-month replication",
        "",
        f"The preregistered UTC interval runs from {result['start']} through {result['stop_exclusive']} (exclusive). Fink lists {coverage['science_alerts_in_statistics']:,} processed alerts across {coverage['observing_nights_in_statistics']} observing-night rows. The anomaly archive contains {len(set(coverage['anomalies_per_utc_date']))} UTC dates, {result['anomaly_alerts']} packets and {result['distinct_anomaly_objects']} distinct ZTF IDs. Both coverage thresholds passed. These two services count different products. [Archived coverage](../../coverage/coverage.json) · [Fink statistics](https://doc.ztf.fink-broker.org/services/search/statistics/)",
        "",
        f"The frozen sampling rule selected {len(objects)} IDs on {len(sample_dates)} UTC dates. {len(overlap)} IDs overlap the October 2024 pilot. The original {metrics['control_agreements']}/{metrics['control_count']} control outcomes reproduced on the current snapshot. The sample is new, while the seven controls are reused. [Registered protocol](../../protocol.md) · [Control and object measurements](report.md)",
        "",
        "| First-pass disposition | Objects |",
        "|---|---:|",
    ]
    lines.extend(f"| {reason} | {n} |" for reason, n in sorted(stages.items()))
    lines.extend(
        [
            "",
            f"Among {len(eligible)} first-pass objects, the 5″ VSX check closed {len(closed)} as unique variable-star associations within 1″. **{len(survivors)} remained after the frozen rules.** [Machine-readable stage counts](summary.json) · [Bounded review packet](review_packet.json)",
            "",
        ]
    )
    if eligible:
        lines += [
            "| First-pass case | Anchor residual | Usable UTC dates | VSX | Downstream disposition |",
            "|---|---|---:|---|---|",
        ]
        for r in eligible:
            m = (
                ", ".join(
                    f"[{n['name']}]({n['source']}) ({n['type']}, {n['separation_arcsec']:.2f}″)"
                    for n in r["vsx_matches"]
                )
                or "No entry within 5″"
            )
            sign = "positive" if r["isdiffpos"] in ("t", "1", 1) else "negative"
            lines.append(
                f"| [{r['id']}]({pilot.PORTAL}{r['id']}) | {sign} | {r['usable_nights']} | {m} | {r['adjudication']} |"
            )
        lines += [""]
    if failed:
        lines.extend(
            [
                "## Broker service failures",
                "",
                "A separate coverage audit examined the archived SIMBAD cone for every historical broker crossmatch failure. Entries with no SIMBAD neighbor were queried against VSX. This exploratory audit leaves the frozen sample dispositions above intact.",
                "",
                "| Object | Broker status | Independent lookup |",
                "|---|---|---|",
            ]
        )
        for r in failed:
            if r["simbad_neighbors"]:
                info = ", ".join(
                    f"{n['name']} ({n['type']}, {n['separation_arcsec']:.2f}″)"
                    for n in r["simbad_neighbors"]
                )
                note = f"SIMBAD · {info}"
            else:
                info = ", ".join(
                    f"[{n['name']}]({n['source']}) ({n['type']}, {n['separation_arcsec']:.2f}″)"
                    for n in r["vsx_matches"]
                )
                note = (
                    f"VSX · {info}"
                    if info
                    else "No SIMBAD neighbor within 3″ or VSX entry within 5″"
                )
            lines.append(
                f"| [{r['id']}]({pilot.PORTAL}{r['id']}) | {r['cdsxmatch']} | {note} |"
            )
        lines.append("")
    lines.extend(
        [
            f"The archive returned {len(histories):,} history rows across controls and exploration objects. {incomplete} contained broker metadata with null object ID, time, position, candidate ID and all photometric fields. Those rows remain in the raw archive and contribute zero observations to the analysis. No alert-selection or classification threshold changed.",
            "",
            "Source ID overlap alone cannot exclude all astrophysical aliases. Weather and observing cadence affect coverage. Broker crossmatch failures are listed separately and may hide additional sources. These counts measure a fixed 20-object sample, so multiple unreviewed objects and observing dates remain in the month.",
            "",
        ]
    )
    (out / "summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experiment",
        type=Path,
        default=pilot.ROOT / "experiments/replication/2026-10",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--fetch-audit",
        action="store_true",
        help="capture missing VSX checks for ambiguous broker failures",
    )
    args = parser.parse_args()
    out = args.out or args.experiment / "results"
    pilot.outputs(pilot.reconcile_snapshot(args.experiment / "snapshot"), out)
    summarize(args.experiment, out, args.fetch_audit)
