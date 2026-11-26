"""Route Fink alerts into three observational questions; replay a registered snapshot."""

import argparse
import csv
import json
import math
from pathlib import Path

import pilot

ROUTES = ("onset", "recurrent", "fading")
BANDS = {1: "g", 2: "r"}
QUESTIONS = {
    "onset": "Did a source newly brighten here, and was the event already reported?",
    "recurrent": "Do separate bright detections trace distinct episodes or ongoing variability?",
    "fading": "Does a source actually fade at these coordinates, or can motion or blending explain the residuals?",
}


def check_replay(a, b, path):
    if type(a) is not type(b):
        raise ValueError(f"Type mismatch at {path}")
    if isinstance(a, dict) and a.keys() == b.keys():
        for k in a:
            check_replay(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            check_replay(x, y, f"{path}[{i}]")
    elif isinstance(a, float):
        if not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-8):
            raise ValueError(f"{path}: {a} differs from {b}")
    elif a != b:
        raise ValueError(f"Mismatch at {path}")


def sign(row):
    value = str(row.get("i:isdiffpos")).lower()
    if value in {"t", "true", "1"}:
        return 1
    if value in {"f", "false", "0"}:
        return -1
    return 0


def quality(row):
    return (
        row.get("i:fid") in BANDS
        and (row.get("i:drb") or 0) >= 0.8
        and 0 < (row.get("i:sigmapsf") or 99) <= 0.2
        and sign(row) != 0
    )


def recent_pairs(dates):
    """Adjacent pairs suffice to find two non-overlapping groups separated by 90 days."""
    ordered = sorted(dates)
    pairs = [(a, b) for a, b in zip(ordered, ordered[1:]) if (b - a).days <= 30]
    if not pairs or (pairs[-1][0] - pairs[0][1]).days < 90:
        return None
    years = sorted({a.year for a, b in pairs if a.year == b.year})
    return {
        "score": [len(years), len(ordered)],
        "years_with_pair": years,
        "example_dates": [
            d.isoformat() for pair in (pairs[0], pairs[-1]) for d in pair
        ],
    }


def measure(anchor, history):
    """Use measurements only. Catalog labels, future points and expected classes are excluded."""
    end = anchor["i:jd"]
    good = [r for r in history if r["i:jd"] <= end and quality(r)]
    dates = {b: {1: set(), -1: set()} for b in BANDS}
    for r in good:
        dates[r["i:fid"]][sign(r)].add(pilot.day_of(r["i:jd"]))
    start = anchor.get("i:jdstarthist")
    age = (
        end - start
        if isinstance(start, (float, int)) and math.isfinite(start) and start > 0
        else None
    )
    pos = dates[1][1] | dates[2][1]
    neg = dates[1][-1] | dates[2][-1]
    span = max(r["i:jd"] for r in good) - min(r["i:jd"] for r in good) if good else None
    eligible = {}
    result = {
        "anchor_quality": quality(anchor),
        "reported_age_days": age,
        "quality_span_days": span,
        "dates": len(pos | neg),
        "positive_dates": len(pos),
        "negative_dates": len(neg),
        "by_filter": {
            BANDS[b]: {"positive": len(v[1]), "negative": len(v[-1])}
            for b, v in dates.items()
        },
        "eligible": eligible,
    }
    if not quality(anchor):
        return result
    if (
        sign(anchor) == 1
        and age is not None
        and 0 <= age <= 90
        and span is not None
        and span <= 90
        and len(pos) >= 3
        and dates[1][1]
        and dates[2][1]
    ):
        eligible["onset"] = {
            "score": [len(dates[1][1] & dates[2][1]), len(pos)],
            "first_positive_date": min(pos).isoformat(),
            "last_positive_date": max(pos).isoformat(),
        }
    for band, signs in dates.items():
        bright, faint = signs[1], signs[-1]
        recurrent = recent_pairs(bright)
        fading = None
        if min(len(bright), len(faint)) >= 5 and (max(faint) - min(faint)).days >= 90:
            fading = {
                "score": [len({d.year for d in faint}), min(len(bright), len(faint))],
                "first_negative_date": min(faint).isoformat(),
                "last_negative_date": max(faint).isoformat(),
            }
        for route, evidence in (("recurrent", recurrent), ("fading", fading)):
            if evidence and (
                route not in eligible or evidence["score"] > eligible[route]["score"]
            ):
                eligible[route] = {"filter": BANDS[band], **evidence}
    return result


def project_positions(extra_ledgers=()):
    old = []
    for source in [
        pilot.ROOT / "results/analysis.json",
        pilot.ROOT / "experiments/november-search/results/analysis.json",
    ]:
        for r in json.loads(source.read_text())["objects"]:
            old.append((r["id"], r["ra"], r["dec"]))
    with (pilot.ROOT / "experiments/october-census/results/outcomes.csv").open() as f:
        for r in csv.DictReader(f):
            old.append((r["id"], float(r["anchor_ra_deg"]), float(r["anchor_dec_deg"])))
    for ledger in extra_ledgers:
        for r in json.loads(ledger.read_text())["objects"]:
            old.append((r["id"], r["ra"], r["dec"]))
    return old


def overlap(row, positions):
    for name, ra, dec in positions:
        if row["id"] == name:
            return f"same ID as {name}"
        if (
            ra is not None
            and dec is not None
            and math.hypot(*pilot.offsets(row["ra"], row["dec"], ra, dec)) <= 5
        ):
            return f"within 5 arcsec of {name}"
    return None


def already_reviewed(paths, anchors):
    current = {pilot.oid(a): a for a in anchors}
    names = set()
    for path in paths:
        previous = json.loads(path.read_text())
        if "cases" in previous:
            names.update(c["id"] for c in previous["cases"])
        elif "case" in previous and "pairs" in previous:
            if previous["case"]:
                names.add(previous["case"]["id"])
        else:
            raise ValueError(f"Unrecognized review packet {path}")
    if names - current.keys():
        raise ValueError(
            "For packets from a different sample use --exclude-ledger instead"
        )
    return [(n, current[n]["i:ra"], current[n]["i:dec"]) for n in sorted(names)]


def select(anchors, histories, positions, gates, unavailable=()):
    objects = []
    for a in anchors:
        name = pilot.oid(a)
        available = name not in unavailable
        row = {
            "id": name,
            "candid": str(a["i:candid"]),
            "jd": a["i:jd"],
            "date": str(pilot.day_of(a["i:jd"])),
            "ra": a["i:ra"],
            "dec": a["i:dec"],
            "history_available": available,
            **measure(a, histories[name] if available else []),
        }
        row.update(
            prior_overlap=overlap(row, positions), route_decisions={}, selected_for=None
        )
        objects.append(row)
    assigned = []
    for route in ROUTES:
        ranked = sorted(
            (r for r in objects if route in r["eligible"]),
            key=lambda r: (*[-n for n in r["eligible"][route]["score"]], r["id"]),
        )
        filled = False
        for rank, r in enumerate(ranked, 1):
            repeat = r["prior_overlap"] or overlap(r, assigned)
            if not gates[route]:
                status = "calibration failed"
            elif repeat:
                status = repeat
            elif filled:
                status = "backlog"
            else:
                status = "admitted"
                filled = True
                r["selected_for"] = route
                assigned.append((r["id"], r["ra"], r["dec"]))
            r["route_decisions"][route] = {"rank": rank, "status": status}
    for r in objects:
        repeat = overlap(r, assigned)
        if repeat:
            for decision in r["route_decisions"].values():
                if decision["status"] == "backlog":
                    decision["status"] = repeat
    return objects


def packet(objects, anchors, gates):
    cases = []
    by_id = {pilot.oid(a): a for a in anchors}
    for r in objects:
        route = r["selected_for"]
        if not route:
            continue
        a = by_id[r["id"]]
        url = pilot.PORTAL + r["id"]
        cases.append(
            {
                "id": r["id"],
                "question": QUESTIONS[route],
                "evidence": [
                    {
                        "id": "position",
                        "source": url,
                        "detail": f"ICRS {r['ra']}, {r['dec']} deg at JD {r['jd']}, UTC {r['date']}; candid {r['candid']}.",
                    },
                    {
                        "id": "anchor",
                        "source": url,
                        "detail": f"{BANDS[a['i:fid']]} difference magnitude {a['i:magpsf']} ± {a['i:sigmapsf']}; sign {sign(a)}, DRB {a['i:drb']}. Reported first-detection age {r['reported_age_days']} days.",
                    },
                    {
                        "id": "history",
                        "source": url,
                        "detail": f"Through the anchor: {r['dates']} distinct quality dates; positive/negative dates by filter {r['by_filter']}. Selected on {route}: {r['eligible'][route]}.",
                    },
                    {
                        "id": "context",
                        "source": url,
                        "detail": f"At the anchor Fink reported CDS={a.get('d:cdsxmatch')!r}, TNS={a.get('d:tns')!r}. Check catalog IDs, astrometry and the registry separately.",
                    },
                ],
            }
        )
    return {"gate_passed": any(gates.values()), "cases": cases}


def render(summary, objects):
    coverage = summary["coverage"]
    lines = [
        f"# Fink discovery queue for {coverage['year']}",
        "",
        f"The anomaly archive contained {coverage['in_window_packets']:,} in-window packets on "
        f"{coverage['dates']} UTC dates from {coverage['distinct_ids']:,} distinct ZTF IDs. The "
        f"sample has {len(objects)} IDs. Verified anchor histories were available for "
        f"{len(objects)-summary['history_unavailable']}. {summary['eligible_ids']} satisfy at least one route. "
        f"{summary['prior_overlaps']} sampled positions overlap earlier investigations or controls.",
        "",
        "| Question | Gate | Eligible IDs | Selected ID | Next available IDs in rank order |",
        "|---|---|---:|---|---|",
    ]
    for route in ROUTES:
        queued = sorted(
            (r["route_decisions"][route]["rank"], r["id"])
            for r in objects
            if r["route_decisions"].get(route, {}).get("status") == "backlog"
            and not r["selected_for"]
        )
        top = (
            ", ".join(f"[{name}]({pilot.PORTAL}{name})" for _, name in queued[:3])
            or "—"
        )
        name = summary["admitted"].get(route)
        chosen = f"[{name}]({pilot.PORTAL}{name})" if name else "—"
        lines.append(
            f"| {QUESTIONS[route]} | {'pass' if summary['route_passed'][route] else 'closed'} | "
            f"{summary['eligible_by_route'][route]} | {chosen} | {top} |"
        )
    lines += [
        "",
        "Route counts can overlap. Each selected ID occupies one review slot. The complete "
        "ledger gives the scores, exclusions and backlog for every sampled ID. Coordinates use "
        "ICRS degrees at the alert observation epoch. The archive ranks only the nightly anomaly "
        "selections. This sample does not estimate rates in the full ZTF alert stream. "
        "[Object ledger](ledger.json) · [Controls](calibration.json) · [Reviewer packet](packet.json)",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--reviewed-packet",
        action="append",
        default=[],
        type=Path,
        help="Advance the same sample past a previously reviewed packet; repeat for each batch",
    )
    parser.add_argument(
        "--exclude-ledger",
        action="append",
        default=[],
        type=Path,
        help="Previous compass ledger, repeatable",
    )
    args = parser.parse_args()
    if any(
        p.resolve().parent == args.out.resolve()
        for p in args.reviewed_packet + args.exclude_ledger
    ):
        parser.error("Use a new output directory to preserve earlier decisions")
    import compass_archive
    import compass_calibration

    anchors, histories, coverage = compass_archive.load(args.snapshot)
    calibration = compass_calibration.run()
    gates = {
        k: calibration["original"]["passed"] and calibration["routes"][k]["passed"]
        for k in ROUTES
    }
    prior = project_positions(args.exclude_ledger) + [
        (r["id"], r.get("ra"), r.get("dec"))
        for r in compass_calibration.control_positions()
    ]
    reviewed = already_reviewed(args.reviewed_packet, anchors)
    prior.extend(reviewed)
    rows = select(anchors, histories, prior, gates, coverage["unavailable_history_ids"])
    summary = {
        "coverage": coverage,
        "sampled_ids": len(rows),
        "quality_anchor": sum(r["anchor_quality"] for r in rows),
        "prior_overlaps": sum(bool(r["prior_overlap"]) for r in rows),
        "eligible_ids": sum(bool(r["eligible"]) for r in rows),
        "history_unavailable": sum(not r["history_available"] for r in rows),
        "no_route_with_history": sum(
            r["history_available"] and not r["eligible"] for r in rows
        ),
        "eligible_by_route": {k: sum(k in r["eligible"] for r in rows) for k in ROUTES},
        "admitted": {r["selected_for"]: r["id"] for r in rows if r["selected_for"]},
        "route_passed": gates,
    }
    if reviewed:
        summary["previously_reviewed_packet_ids"] = [v[0] for v in reviewed]
    args.out.mkdir(parents=True, exist_ok=True)
    pilot.json_write(args.out / "ledger.json", {"summary": summary, "objects": rows})
    pilot.json_write(args.out / "calibration.json", calibration)
    pilot.json_write(args.out / "packet.json", packet(rows, anchors, gates))
    (args.out / "report.md").write_text(render(summary, rows))
    print(json.dumps({k: v for k, v in summary.items() if k != "coverage"}, indent=2))


if __name__ == "__main__":
    main()
