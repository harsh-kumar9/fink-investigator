"""Compare positive and negative difference centroids on temporally disjoint date pairs."""

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from datetime import timedelta
import json
import math
from pathlib import Path
from statistics import median

import compass
import compass_archive
import pilot


def coordinate(row):
    ra, dec = row.get("i:ra"), row.get("i:dec")
    return (
        isinstance(ra, (int, float))
        and 0 <= ra < 360
        and isinstance(dec, (int, float))
        and -90 <= dec <= 90
    )


def daily(anchor, history, band):
    groups = defaultdict(list)
    for row in history:
        if (
            row["i:jd"] <= anchor["i:jd"]
            and row["i:fid"] == band
            and compass.quality(row)
            and coordinate(row)
        ):
            date = pilot.day_of(row["i:jd"])
            x, y = pilot.offsets(
                row["i:ra"], row["i:dec"], anchor["i:ra"], anchor["i:dec"]
            )
            groups[date, compass.sign(row)].append((x, y, str(row["i:candid"])))
    return {
        (date, sign): (
            median(x for x, _, _ in vals),
            median(y for _, y, _ in vals),
            sorted(c for _, _, c in vals),
        )
        for (date, sign), vals in groups.items()
    }


def date_pairs(positions):
    positive = sorted(d for d, sign in positions if sign > 0)
    negative = sorted(d for d, sign in positions if sign < 0)
    edges = []
    for p in positive:
        lo = bisect_left(negative, p - timedelta(days=30))
        hi = bisect_right(negative, p + timedelta(days=30))
        edges.extend((abs((p - n).days), p, n) for n in negative[lo:hi])
    used = set()
    pairs = []
    for gap, p, n in sorted(edges):
        if p in used or n in used:
            continue
        used.update((p, n))
        px, py, pc = positions[p, 1]
        nx, ny, nc = positions[n, -1]
        pairs.append(
            {
                "positive": p.isoformat(),
                "negative": n.isoformat(),
                "lag_days": gap,
                "east": px - nx,
                "north": py - ny,
                "positive_candids": pc,
                "negative_candids": nc,
            }
        )
    return sorted(
        pairs,
        key=lambda r: (
            max(r["positive"], r["negative"]),
            min(r["positive"], r["negative"]),
        ),
    )


def assess(anchor, history, band):
    positions = daily(anchor, history, band)
    pairs = date_pairs(positions)
    n = 2 * len(pairs) // 3
    training = pairs[:n]
    end_training = (
        max(training[-1]["positive"], training[-1]["negative"]) if training else ""
    )
    test = [p for p in pairs[n:] if min(p["positive"], p["negative"]) > end_training]
    stats = {
        "positive_dates": sum(s == 1 for _, s in positions),
        "negative_dates": sum(s == -1 for _, s in positions),
        "pairs": len(pairs),
        "training_pairs": n,
        "test_pairs": len(test),
        "boundary_pairs_dropped": len(pairs) - n - len(test),
        "training_end": end_training or None,
        "test_start": min(
            (min(p["positive"], p["negative"]) for p in test), default=None
        ),
    }
    evidence = {"training": training, "held_out": test}
    if n < 8 or len(test) < 4:
        return {**stats, "status": "insufficient"}, evidence
    tx, ty = (median(p[k] for p in training) for k in ("east", "north"))
    vx, vy = (median(p[k] for p in test) for k in ("east", "north"))
    length = math.hypot(tx, ty)
    dot, cross = tx * vx + ty * vy, tx * vy - ty * vx
    angle = (
        math.degrees(math.atan2(abs(cross), dot))
        if length and math.hypot(vx, vy)
        else None
    )
    projection = dot / length if length else None
    passed = (
        length >= 0.45
        and projection is not None
        and projection >= 0.35
        and angle is not None
        and angle <= 45
    )
    return {
        **stats,
        "training_east_north_arcsec": [tx, ty],
        "test_east_north_arcsec": [vx, vy],
        "training_length_arcsec": length,
        "test_projection_arcsec": projection,
        "angle_deg": angle,
        "score_arcsec": min(length, projection) if projection is not None else None,
        "status": "passed" if passed else "below threshold",
    }, evidence


def controls():
    values = {}
    for case, count, expected in [
        ("stable", 24, "passed"),
        ("common_drift", 24, "below threshold"),
        ("reversed", 24, "below threshold"),
        ("sparse", 11, "insufficient"),
    ]:
        anchor = {"i:ra": 359.99985, "i:dec": 60.0, "i:jd": 2459000 + count * 60}
        rows = []
        for i in range(count):
            for sign in (1, -1):
                time = 60 * i + int((sign == 1) == (i % 2 == 0))
                if case == "common_drift":
                    x, y = 0.08 * time, 0.02 * time
                else:
                    dx = -1 if case == "reversed" and i >= 16 else 1
                    x, y = (0.6 * dx, 0.45 * dx) if sign == 1 else (0.0, 0.0)
                rows.append(
                    {
                        "i:jd": 2459000 + time,
                        "i:fid": 1,
                        "i:candid": len(rows) + 1,
                        "i:drb": 1.0,
                        "i:sigmapsf": 0.05,
                        "i:isdiffpos": "t" if sign == 1 else "f",
                        "i:ra": (anchor["i:ra"] + x / (3600 * 0.5)) % 360,
                        "i:dec": anchor["i:dec"] + y / 3600,
                    }
                )
        result, grouped = assess(anchor, rows, 1)
        dates = [
            {p[k] for p in g for k in ("positive", "negative")}
            for g in grouped.values()
        ]
        good = (
            result["status"] == expected
            and len(dates[0] & dates[1]) == 0
            and sum(len(d) for d in dates) == 2 * count
        )
        if case == "stable":
            good = good and all(
                math.isclose(a, b, abs_tol=1e-8)
                for a, b in zip(result["training_east_north_arcsec"], [0.60, 0.45])
            )
        values[case] = {"expected": expected, "agreement": good, "result": result}
    return {
        "kind": "synthetic geometry controls",
        "passed": all(r["agreement"] for r in values.values()),
        "cases": values,
    }


def run(snapshot, ledger):
    gate = controls()
    anchors, histories, coverage = compass_archive.load(snapshot)
    objects = json.loads(ledger.read_text())["objects"]
    original = {r["id"]: r for r in objects}
    if len(original) != len(objects) or len(anchors) != len(original):
        raise ValueError("Snapshot and ledger membership differ")
    prior_slots = [(r["id"], r["ra"], r["dec"]) for r in objects if r["selected_for"]]
    unavailable = set(coverage["unavailable_history_ids"])
    ranked = []
    records = []
    details = {}
    for a in anchors:
        name = pilot.oid(a)
        base = original.get(name)
        if (
            base is None
            or str(a["i:candid"]) != base["candid"]
            or (name not in unavailable) != base["history_available"]
        ):
            raise ValueError(f"Ledger does not match the snapshot at {name}")
        record = {
            "id": name,
            "candid": str(a["i:candid"]),
            "jd": a["i:jd"],
            "date": str(pilot.day_of(a["i:jd"])),
            "ra": a["i:ra"],
            "dec": a["i:dec"],
            "bands": {},
            "overlap": base["prior_overlap"] or compass.overlap(base, prior_slots),
        }
        records.append(record)
        if name in unavailable:
            record["status"] = "anchor missing"
            continue
        if not coordinate(a):
            raise ValueError(f"Invalid astrometric origin at {name}")
        best = None
        for band, label in compass.BANDS.items():
            stats, pairs = assess(a, histories[name], band)
            record["bands"][label] = stats
            if stats["status"] == "passed" and (
                best is None or stats["score_arcsec"] > best
            ):
                best = stats["score_arcsec"]
                record["best_band"] = label
                details[name] = pairs
        if best is not None:
            ranked.append((-best, name, record))
        else:
            record["status"] = (
                "below threshold"
                if any(
                    b["status"] == "below threshold" for b in record["bands"].values()
                )
                else "insufficient"
            )
    winner = None
    for rank, (_, name, record) in enumerate(sorted(ranked), 1):
        record["rank"] = rank
        record["status"] = (
            "calibration failed"
            if not gate["passed"]
            else "overlap" if record["overlap"] else "backlog" if winner else "admitted"
        )
        if record["status"] == "admitted":
            winner = record
    measurable = sum(
        any(b["status"] != "insufficient" for b in r["bands"].values()) for r in records
    )
    summary = {
        "year": coverage["year"],
        "sample": len(records),
        "verified_histories": len(records) - len(unavailable),
        "measurable_ids": measurable,
        "passed_ids": len(ranked),
        "dispositions": dict(Counter(r["status"] for r in records)),
        "admitted": winner["id"] if winner else None,
    }
    packet = {"case": winner, "pairs": details[winner["id"]] if winner else {}}
    return gate, {"summary": summary, "objects": records}, packet


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Recompute and compare with all three frozen outputs without rewriting",
    )
    args = parser.parse_args()
    control, result, packet = run(args.snapshot, args.ledger)
    if not args.verify:
        args.out.mkdir(exist_ok=True, parents=True)
    for name, data in [("controls", control), ("results", result), ("packet", packet)]:
        path = args.out / (name + ".json")
        if args.verify:
            compass.check_replay(data, json.loads(path.read_text()), str(path))
        else:
            pilot.json_write(path, data)
    if args.verify and not control["passed"]:
        raise ValueError("Synthetic geometry checks failed")
    print(json.dumps(result["summary"], indent=2))
