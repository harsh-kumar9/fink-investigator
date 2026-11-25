"""Replay or refresh the Fink anomaly pilot. Python 3.10+; refresh uses curl."""

import argparse
import csv
import gzip
import fcntl
import hashlib
import io
import json
import math
import os
from collections import defaultdict
from contextlib import contextmanager
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
import re
from statistics import median
import subprocess
from urllib.parse import urlencode

import vsx

ROOT = Path(__file__).resolve().parent
FINK = "https://api.ztf.fink-portal.org/api/v1"
IRSA = "https://irsa.ipac.caltech.edu/TAP/sync"
SIMBAD = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"
PORTAL = "https://ztf.fink-portal.org/"
COLUMNS = "i:objectId,i:candid,i:ra,i:dec,i:jd,i:drb,i:fid,i:magpsf,i:sigmapsf,i:isdiffpos,d:cdsxmatch,d:tns"
CONTROL_END = date(2026, 2, 1)
PM_MIN = 150  # mas / year
MATCH_RADIUS = 2.5  # arcsec, distance to the star's path
CONE = 60  # arcsec around each Fink anchor at Gaia reference epoch
REVIEW_LIMIT = 3


def utc_jd(day):
    return (day - date(1970, 1, 1)).days + 2440587.5


def day_of(jd):
    return date(1970, 1, 1) + timedelta(days=math.floor(jd - 2440587.5))


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def control_labels():
    labels = json.loads((ROOT / "data/controls.json").read_text())
    positives = [r for r in labels if r["expected"] is True]
    negatives = [r for r in labels if r["expected"] is False]
    if (
        not positives
        or not negatives
        or len(positives) + len(negatives) != len(labels)
        or len({r["id"] for r in labels}) != len(labels)
        or len({r["expected_gaia_dr3_source_id"] for r in positives}) != len(positives)
    ):
        raise ValueError(
            "controls require independent, unique positive and negative sources"
        )
    return labels


def oid(row):
    name = row["i:objectId"]
    if not re.fullmatch(r"ZTF\d{2}[a-z]{7}", name):
        raise ValueError(f"unexpected Fink object ID {name!r}")
    return name


def select_alerts(rows, n, start, stop):
    # Choose each object's most anomalous packet. Spread investigations through the month.
    unique = {}
    for row in rows:
        if not (utc_jd(start) <= row["i:jd"] < utc_jd(stop)):
            continue
        name = oid(row)
        score = row["d:anomaly_score"]
        if not isinstance(score, (float, int)) or not math.isfinite(score):
            raise ValueError(f"invalid anomaly score on {name}")
        if name not in unique or (score, row["i:jd"]) < (
            unique[name]["d:anomaly_score"],
            unique[name]["i:jd"],
        ):
            unique[name] = row
    ordered = sorted(unique.values(), key=lambda r: (r["i:jd"], oid(r)))
    count = min(n, len(ordered))
    if count <= 1:
        return ordered[:count]
    return [ordered[round(i * (len(ordered) - 1) / (count - 1))] for i in range(count)]


def curl(url, extra_args):
    args = [
        "curl",
        "--fail",
        "--silent",
        "--show-error",
        "--max-time",
        "120",
        "--retry",
        "2",
        "--retry-delay",
        "1",
        "--max-filesize",
        "25000000",
    ]
    return subprocess.run(
        args + extra_args + [url], check=True, stdout=subprocess.PIPE
    ).stdout


def fetch_response(url, params, tap=False):
    if tap:
        args = [
            "--get",
            "--data-urlencode",
            "REQUEST=doQuery",
            "--data-urlencode",
            "LANG=ADQL",
            "--data-urlencode",
            "FORMAT=csv",
            "--data-urlencode",
            f"QUERY={params['QUERY']}",
        ]
    else:
        args = [
            "--header",
            "Content-Type: application/json",
            "--data-binary",
            json.dumps(params),
        ]
    data = curl(url, args)
    if tap:
        parse_csv(data)
    elif not isinstance(json.loads(data), list):
        raise ValueError(f"invalid Fink response from {url}")
    return data


def parse_csv(data):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8")))
    if reader.fieldnames is None or not (
        {"source_id", "main_id"} & set(reader.fieldnames)
    ):
        raise ValueError(f"unexpected TAP response: {data[:900]!r}")
    rows = list(reader)
    if len(rows) >= 2000:
        raise ValueError("TAP query reached its row cap")
    return rows


def store_raw(folder, manifest, name, url, params, tap=False):
    if name in manifest["sources"]:
        prior = manifest["sources"][name]
        if prior["params"] != params or prior["url"] != url:
            raise ValueError(
                f"request changed for partial snapshot {name}; select a fresh destination"
            )
        raw = gzip.decompress((folder / name).read_bytes())
        if hashlib.sha256(raw).hexdigest() != prior["sha256_uncompressed"]:
            raise ValueError(f"partial snapshot checksum mismatch: {name}")
        return parse_csv(raw) if tap else json.loads(raw)
    raw = fetch_response(url, params, tap)
    (folder / name).write_bytes(gzip.compress(raw, mtime=0))
    manifest["sources"][name] = {
        "url": url,
        "params": params,
        "method": "GET" if tap else "POST",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
    }
    checkpoint = folder / "manifest.tmp"
    json_write(checkpoint, manifest)
    os.replace(checkpoint, folder / "manifest.json")
    return parse_csv(raw) if tap else json.loads(raw)


@contextmanager
def locked_directory(folder):
    folder.mkdir(parents=True, exist_ok=True)
    fd = os.open(folder, os.O_RDONLY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(fd)


def histories(rows):
    groups = defaultdict(dict)
    for row in rows:
        name = oid(row)
        candid = row["i:candid"]
        if not isinstance(candid, int) or candid <= 0:
            raise ValueError(f"missing numeric candidate id on {name}")
        if candid in groups[name] and groups[name][candid] != row:
            raise ValueError(f"conflicting duplicate {candid}")
        groups[name][candid] = row
    return {
        name: sorted(vals.values(), key=lambda r: r["i:jd"])
        for name, vals in groups.items()
    }


def cone_clause(anchors, radius):
    clauses = []
    for row in anchors:
        ra, dec = row["i:ra"], row["i:dec"]
        if not (0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("invalid sky position")
        clauses.append(
            f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{radius/3600:.10f}))=1"
        )
    return " OR ".join(clauses)


def refresh(folder, start, stop, count):
    if folder.exists():
        raise ValueError(
            f"snapshot directory exists: {folder}. Supply a new destination."
        )
    if not 1 <= (stop - start).days <= 31 or not 1 <= count <= 20:
        raise ValueError("choose 1–31 days and 1–20 objects")
    folder.parent.mkdir(parents=True, exist_ok=True)
    work = folder.with_name(folder.name + ".partial")
    with locked_directory(work):
        manifest = {
            "version": 1,
            "start": str(start),
            "stop_exclusive": str(stop),
            "max_objects": count,
            "sources": {},
        }
        if (work / "manifest.json").exists():
            previous = json.loads((work / "manifest.json").read_text())
            if any(
                manifest[k] != previous.get(k)
                for k in ("version", "start", "stop_exclusive", "max_objects")
            ):
                raise ValueError(
                    "partial snapshot has different bounds; select a fresh destination"
                )
            manifest = previous
        # Ask for headroom above the documented ten per night. Fail on a truncated response.
        n = 10 * ((stop - start).days + 2) + 1
        rows = store_raw(
            work,
            manifest,
            "anomaly.json.gz",
            FINK + "/anomaly",
            {
                "start_date": str(start),
                "stop_date": str(stop),
                "n": n,
                "columns": COLUMNS + ",d:anomaly_score",
                "output-format": "json",
            },
        )
        if len(rows) >= n:
            raise ValueError("anomaly response may be truncated")
        targets = select_alerts(rows, count, start, stop)
        controls = [r["id"] for r in control_labels()]
        ids = sorted(set(controls).union(oid(r) for r in targets))
        for i in range(0, len(ids), 4):
            store_raw(
                work,
                manifest,
                f"history-{i:02}.json.gz",
                FINK + "/objects",
                {
                    "objectId": ",".join(ids[i : i + 4]),
                    "columns": COLUMNS,
                    "output-format": "json",
                },
            )
        history = histories(
            [
                row
                for file in work.glob("history-*.json.gz")
                for row in json.loads(gzip.decompress(file.read_bytes()))
            ]
        )
        anchors = [
            max(
                (r for r in history[name] if r["i:jd"] < utc_jd(CONTROL_END)),
                key=lambda r: r["i:jd"],
            )
            for name in controls
        ]
        anchors += targets
        # IRSA supports one indexed cone per request. SIMBAD accepts a small union of cones.
        for i, anchor in enumerate(anchors):
            subset = [anchor]
            gaia = (
                "SELECT TOP 2000 source_id,ra,dec,ref_epoch,pmra,pmdec,pmra_error,pmdec_error "
                "FROM gaia_dr3_source WHERE "
                + cone_clause(subset, CONE)
                + " AND (pmra*pmra+pmdec*pmdec)>=22500"
            )
            store_raw(
                work, manifest, f"gaia-{i:02}.csv.gz", IRSA, {"QUERY": gaia}, tap=True
            )
        for i in range(0, len(targets), 8):
            simbad = (
                "SELECT TOP 2000 main_id,otype,ra,dec FROM basic WHERE "
                + cone_clause(targets[i : i + 8], 3)
            )
            store_raw(
                work,
                manifest,
                f"simbad-{i:02}.csv.gz",
                SIMBAD,
                {"QUERY": simbad},
                tap=True,
            )
        json_write(work / "manifest.json", manifest)
        # Verify against the source bytes and run analysis before publishing a complete snapshot.
        evaluate(work)
        os.rename(work, folder)


def snapshot(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    raw = {}
    for name, meta in manifest["sources"].items():
        if Path(name).name != name:
            raise ValueError("invalid snapshot filename")
        data = gzip.decompress((folder / name).read_bytes())
        if hashlib.sha256(data).hexdigest() != meta["sha256_uncompressed"]:
            raise ValueError(f"checksum mismatch: {name}")
        raw[name] = parse_csv(data) if ".csv." in name else json.loads(data)
    return manifest, raw


def offsets(ra, dec, ra0, dec0):
    # Signed arcseconds in the local tangent plane, with the RA wrap handled explicitly.
    return ((ra - ra0 + 180) % 360 - 180) * 3600 * math.cos(math.radians(dec0)), (
        dec - dec0
    ) * 3600


def stellar_matches(anchor, stars):
    matches = []
    epoch = 2000 + (anchor["i:jd"] - 2451545.0) / 365.25
    for s in stars.values():
        ra0, dec0 = float(s["ra"]), float(s["dec"])
        x, y = offsets(anchor["i:ra"], anchor["i:dec"], ra0, dec0)
        if math.hypot(x, y) > CONE:
            continue
        vx, vy = float(s["pmra"]), float(s["pmdec"])
        pm = math.hypot(vx, vy)
        if pm < PM_MIN or not s["pmra_error"] or not s["pmdec_error"]:
            continue
        if pm < 5 * math.hypot(float(s["pmra_error"]), float(s["pmdec_error"])):
            continue
        dt = epoch - float(s["ref_epoch"])
        dx, dy = dt * vx / 1000, dt * vy / 1000
        # Both the moving source and its former position can leave subtraction residuals.
        along = min(1, max(0, (x * dx + y * dy) / (dx * dx + dy * dy))) if dt else 0
        dist = math.hypot(x - along * dx, y - along * dy)
        if dist <= MATCH_RADIUS:
            matches.append(
                {
                    "source_id": s["source_id"],
                    "pm_mas_yr": round(pm, 1),
                    "distance_to_track_arcsec": round(dist, 3),
                }
            )
    return sorted(matches, key=lambda r: r["distance_to_track_arcsec"])


def displacement_evidence(star, rows):
    """Compare residual locations to the Gaia track using only pre-alert measurements."""
    vx = float(star["pmra"]) / 1000
    vy = float(star["pmdec"]) / 1000
    ref_epoch = float(star["ref_epoch"])
    negative = [r for r in rows if str(r["i:isdiffpos"]).lower() in {"f", "false", "0"}]
    years = []
    distances = []
    for r in negative:
        x, y = offsets(r["i:ra"], r["i:dec"], float(star["ra"]), float(star["dec"]))
        elapsed = 2000 + (r["i:jd"] - 2451545.0) / 365.25 - ref_epoch
        years.append(ref_epoch + (vx * x + vy * y) / (vx * vx + vy * vy))
        distances.append(math.hypot(x - vx * elapsed, y - vy * elapsed))
    return {
        "negative_rows": len(negative),
        "negative_utc_dates": len({day_of(r["i:jd"]) for r in negative}),
        "median_projected_epoch": round(median(years), 2) if years else None,
        "median_distance_to_contemporaneous_star_arcsec": (
            round(median(distances), 2) if distances else None
        ),
    }


def matched_simbad(anchor, simbad):
    found = []
    for star in simbad.values():
        x, y = offsets(
            float(star["ra"]), float(star["dec"]), anchor["i:ra"], anchor["i:dec"]
        )
        sep = math.hypot(x, y)
        if sep <= 3:
            found.append(
                {
                    "name": star["main_id"],
                    "type": star["otype"],
                    "separation_arcsec": round(sep, 3),
                }
            )
    return sorted(found, key=lambda r: r["separation_arcsec"])


def evaluate(folder):
    manifest, raw = snapshot(folder)
    start, stop = date.fromisoformat(manifest["start"]), date.fromisoformat(
        manifest["stop_exclusive"]
    )
    history = histories(
        [r for name, rows in raw.items() if name.startswith("history-") for r in rows]
    )
    stars = {
        s["source_id"]: s
        for name, rows in raw.items()
        if name.startswith("gaia-")
        for s in rows
    }
    simbad = {
        s["main_id"]: s
        for name, rows in raw.items()
        if name.startswith("simbad-")
        for s in rows
    }
    labels = control_labels()
    calibration = []
    for c in labels:
        anchor = max(
            (r for r in history[c["id"]] if r["i:jd"] < utc_jd(CONTROL_END)),
            key=lambda r: r["i:jd"],
        )
        matches = stellar_matches(anchor, stars)
        observed = {m["source_id"] for m in matches}
        agrees = (
            (c["expected_gaia_dr3_source_id"] in observed)
            if c["expected"]
            else not observed
        )
        calibration.append(
            {**c, "candid": str(anchor["i:candid"]), "matches": matches, "pass": agrees}
        )
    gate = all(r["pass"] for r in calibration)
    anomalies = [
        r for r in raw["anomaly.json.gz"] if utc_jd(start) <= r["i:jd"] < utc_jd(stop)
    ]
    targets = select_alerts(anomalies, manifest["max_objects"], start, stop)
    control_ids = {c["id"] for c in labels} | {
        "ZTF23aacafew"
    }  # G 260-1 alias, paper Table 4
    control_gaia_sources = {m["source_id"] for c in calibration for m in c["matches"]}
    anchor_count = len(labels) + len(targets)
    history_count = len({c["id"] for c in labels}.union(oid(r) for r in targets))
    required = {"anomaly.json.gz"} | {
        f"history-{i:02}.json.gz" for i in range(0, history_count, 4)
    }
    required |= {f"gaia-{i:02}.csv.gz" for i in range(anchor_count)}
    required |= {f"simbad-{i:02}.csv.gz" for i in range(0, len(targets), 8)}
    if set(raw) != required:
        raise ValueError(
            f"incomplete snapshot inventory: {sorted(set(raw) ^ required)}"
        )
    review = []
    for a in targets:
        name = oid(a)
        # Never leak data after the alert into the time-series assessment.
        hist = [r for r in history[name] if r["i:jd"] <= a["i:jd"] + 1e-7]
        history_anchor = next((r for r in hist if r["i:candid"] == a["i:candid"]), None)
        if history_anchor is None:
            raise ValueError(f"missing anchor in history for {name}")
        good = [
            r
            for r in hist
            if isinstance(r["i:drb"], (int, float))
            and r["i:drb"] >= 0.8
            and isinstance(r["i:sigmapsf"], (int, float))
            and 0 < r["i:sigmapsf"] <= 0.2
        ]
        nights = len({day_of(r["i:jd"]) for r in good})
        positive_nights = len(
            {
                day_of(r["i:jd"])
                for r in good
                if str(r["i:isdiffpos"]).lower() in {"1", "t", "true"}
            }
        )
        matches = stellar_matches(a, stars)
        for m in matches:
            m.update(displacement_evidence(stars[m["source_id"]], good))
        neighbors = matched_simbad(a, simbad)
        label = str(a["d:cdsxmatch"])
        tns = str(history_anchor.get("d:tns") or "").strip()
        no_tns = tns.lower() in {"nan", "none", "null", "unknown", "", "-99"}
        if name in control_ids or any(
            m["source_id"] in control_gaia_sources for m in matches
        ):
            status = "calibration overlap"
        elif matches:
            status = "moving-star alternative"
        elif label.lower().startswith("fail"):
            status = "broker crossmatch failed"
        elif label != "Unknown" or neighbors or not no_tns:
            status = "catalog context"
        elif (
            nights < 3 or positive_nights < 2 or a["i:drb"] is None or a["i:drb"] < 0.8
        ):
            status = "limited positive support"
        else:
            status = "review"
        review.append(
            {
                "id": name,
                "candid": str(a["i:candid"]),
                "jd": a["i:jd"],
                "ra": a["i:ra"],
                "dec": a["i:dec"],
                "score": a["d:anomaly_score"],
                "drb": a["i:drb"],
                "filter": a["i:fid"],
                "magpsf": a["i:magpsf"],
                "sigmapsf": a["i:sigmapsf"],
                "history_rows": len(hist),
                "usable_rows": len(good),
                "usable_nights": nights,
                "positive_nights": positive_nights,
                "cdsxmatch": label,
                "tns": None if no_tns else tns,
                "isdiffpos": a["i:isdiffpos"],
                "gaia_matches": matches,
                "simbad_neighbors": neighbors,
                "status": status,
            }
        )
    review.sort(key=lambda r: (r["score"], r["id"]))
    return {
        "gate_passed": gate,
        "start": str(start),
        "stop_exclusive": str(stop),
        "anomaly_alerts": len(anomalies),
        "distinct_anomaly_objects": len({oid(r) for r in anomalies}),
        "calibration": calibration,
        "objects": review,
    }


def image_url(r, kind):
    return (
        FINK
        + "/cutouts?"
        + urlencode({"objectId": r["id"], "candid": r["candid"], "kind": kind})
    )


def reconcile_snapshot(folder, fetch=False):
    result = evaluate(folder)
    if fetch:
        with locked_directory(folder):
            return vsx.reconcile(folder, result, lambda url: curl(url, []))
    return vsx.reconcile(folder, result)


def case_packet(r):
    portal = PORTAL + r["id"]
    sim_url = "https://simbad.cfa.harvard.edu/simbad/sim-coo?" + urlencode(
        {
            "Coord": f"{r['ra']} {r['dec']:+}",
            "CooFrame": "ICRS",
            "CooEpoch": "2000",
            "Radius": 3,
            "Radius.unit": "arcsec",
        }
    )
    tns_context = (
        f"stored TNS label {r['tns']}"
        if r["tns"]
        else "TNS label not supplied for this alert; no independent TNS cone search was made"
    )
    polarity = {"t": "positive", "f": "negative", "1": "positive", "0": "negative"}.get(
        str(r["isdiffpos"]).lower(), "unknown-sign"
    )
    packet = {
        "id": r["id"],
        "question": "Does the evidence warrant inspecting this historical anomaly for an ordinary variable, catalog counterpart, or subtraction artifact?",
        "evidence": [
            {
                "id": "alert",
                "source": portal,
                "detail": f"Fink score {r['score']:.4f}, DRB {r['drb']:.3f}, JD {r['jd']}, candid {r['candid']}; RA {r['ra']}, Dec {r['dec']} deg. The {polarity} difference residual has magnitude {r['magpsf']:.3f} ± {r['sigmapsf']:.3f} in ZTF filter {r['filter']}.",
            },
            {
                "id": "history",
                "source": portal,
                "detail": f"{r['history_rows']} detections through the selected alert. {r['usable_rows']} pass DRB >= 0.8 and sigma <= 0.2 mag on {r['usable_nights']} UTC dates, with positive residuals on {r['positive_nights']} dates.",
            },
            {
                "id": "catalog",
                "source": sim_url,
                "detail": f"Broker CDS label {r['cdsxmatch']!r}; {tns_context}; current 3 arcsec SIMBAD query yielded {len(r['simbad_neighbors'])} neighbors. Coverage and positional differences limit this check.",
            },
            {
                "id": "motion",
                "source": "https://gaia.aip.de/metadata/gaiadr3/gaia_source/",
                "detail": f"{len(r['gaia_matches'])} significant Gaia DR3 stars >=150 mas/yr passed the 2.5 arcsec track screen. Slower, uncataloged or poorly measured stars and image artifacts remain plausible.",
            },
            {
                "id": "image",
                "source": image_url(r, "Difference"),
                "detail": "Epoch-specific difference cutout available for inspection alongside science and template. This link alone supplies no visual verdict.",
            },
        ],
    }
    if r["vsx_matches"]:
        packet["question"] = (
            "Which observations support the published variable-star counterpart, and which observations could distinguish a chance superposition or subtraction residual? Keep the new-source lead closed for a unique match within one arcsecond."
        )
        for i, m in enumerate(r["vsx_matches"]):
            packet["evidence"].append(
                {
                    "id": f"vsx-{i+1}",
                    "source": m["source"],
                    "detail": f"VSX OID {m['oid']}, {m['name']}, variability type {m['type']!r}, separation {m['separation_arcsec']:.2f} arcsec. Triage: {r['adjudication']}. A colon in a VSX type indicates uncertainty.",
                }
            )
    else:
        packet["evidence"].append(
            {
                "id": "vsx-cone",
                "source": r["vsx_query"],
                "detail": "The current five-arcsecond VSX cone returned no variable-star entries; other catalogs and image review remain outstanding.",
            }
        )
    return packet


def outputs(result, folder):
    folder.mkdir(parents=True, exist_ok=True)
    json_write(folder / "analysis.json", result)
    eligible = (
        [r for r in result["objects"] if r["status"] == "review"]
        if result["gate_passed"]
        else []
    )
    open_cases = [r for r in eligible if not r["adjudication"].startswith("close")]
    audits = [r for r in eligible if r["adjudication"].startswith("close")]
    batch = (open_cases + audits)[:REVIEW_LIMIT]
    packet = {
        "gate_passed": result["gate_passed"],
        "cases": [case_packet(r) for r in batch],
    }
    json_write(folder / "review_packet.json", packet)
    cal = result["calibration"]
    lines = [
        "# Historical Fink anomaly pilot",
        "",
        f"{result['anomaly_alerts']} alerts from {result['start']} through {result['stop_exclusive']} (exclusive) contain {result['distinct_anomaly_objects']} unique ZTF IDs. The sample covers {len(result['objects'])} objects.",
        "",
        f"The Gaia moving-star screen agreed with {sum(c['pass'] for c in cal)}/{len(cal)} independently sourced control labels. The review gate {'passed' if result['gate_passed'] else 'failed'}.",
        "",
        "| Control | Literature identification | Gaia track offset (arcsec) | Agreement |",
        "|---|---|---:|---|",
    ]
    for c in cal:
        d = (
            ", ".join(str(m["distance_to_track_arcsec"]) for m in c["matches"])
            or "no match"
        )
        lines.append(
            f"| [{c['id']}]({PORTAL}{c['id']}) | [{c['name']}]({c['source']}) | {d} | {c['pass']} |"
        )
    lines += [
        "",
        "## Investigation queue",
        "",
        f"{len(eligible)} objects cleared the first screen. The VSX review found a single variable within 1″ for {len(audits)}. {len(open_cases)} remain for counterpart investigation. The packet contains up to three cases, with open leads first and alias audits next, ordered within each group by anomaly score.",
        "",
        "| Object | Score | Usable / positive UTC dates | Broker match | SIMBAD within 3″ | First-pass triage |",
        "|---|---:|---:|---|---|---|",
    ]
    for r in result["objects"]:
        sim = ", ".join(n["name"] for n in r["simbad_neighbors"]) or "—"
        lines.append(
            f"| [{r['id']}]({PORTAL}{r['id']}) | {r['score']:.4f} | {r['usable_nights']} / {r['positive_nights']} | {r['cdsxmatch']} | {sim} | {r['status']} |"
        )
    transfer = (
        [r for r in result["objects"] if r["status"] == "moving-star alternative"]
        if result["gate_passed"]
        else []
    )
    if transfer:
        lines += ["", "## Moving-star checks outside the control set", ""]
    for r in transfer:
        lines += [f"### {r['id']}", ""]
        for m in r["gaia_matches"]:
            lines.append(
                f"Gaia DR3 {m['source_id']} moves {m['pm_mas_yr']} mas/year. The alert position is {m['distance_to_track_arcsec']}″ from the track segment between J2016.0 and the alert epoch. [Fink measurements]({PORTAL}{r['id']})"
            )
            if m["negative_rows"]:
                lines += [
                    "",
                    f"{m['negative_rows']} of {r['usable_rows']} usable residuals are negative, spanning {m['negative_utc_dates']} UTC dates. Their median closest epoch along the Gaia track is {m['median_projected_epoch']}. The median distance from the star's position predicted at the time of each exposure is {m['median_distance_to_contemporaneous_star_arcsec']}″. An old reference-image position could produce this pattern. The original coadd dates and image alignment provide independent tests of this interpretation. [Science]({image_url(r,'Science')}) · [Template]({image_url(r,'Template')}) · [Difference]({image_url(r,'Difference')})",
                    "",
                ]
        # The frozen review artifacts belong to this alert and are excluded from scoring.
        if r["id"] == "ZTF18abgobyb" and r["candid"] == "2841171223915010017":
            lines += [
                "An independent check identified the star as VB 10. The anchor predicts a crossing around 2018-05-13, inside the 18-frame reference interval of 2018-03-21 through 2018-06-24. The April, July and October 2024 difference stamps each show a pair of opposite-sign residuals with similar orientation. Together the timing and morphology support foreground-star subtraction. [Frozen reference metadata and calculation](moving-star/report.md) · [Three-epoch gallery](../review/ZTF18abgobyb/gallery.md)",
                "",
            ]
    lines += ["", "## Variable-star catalog audit", ""]
    for r in batch:
        vsx_facts = " · ".join(
            f"[{m['name']}]({m['source']}) (VSX {m['oid']}, {m['type']}, {m['separation_arcsec']:.2f}″)"
            for m in r["vsx_matches"]
        )
        lines.extend(
            [
                f"### {r['id']}",
                "",
                f"**{r['adjudication']}**. {vsx_facts or 'No entries in the archived five-arcsecond VSX cone.'}",
                "",
                f"RA {r['ra']}°, Dec {r['dec']}°. Alert {r['candid']} at JD {r['jd']}. Its {r['usable_rows']} usable detections include positive subtraction residuals on {r['positive_nights']} dates. Inspect for an offset nucleus, dipole or recurring variable. Compare with forced photometry and a current coordinate search in TNS before interpreting the signal.",
                "",
                "Epoch cutouts · "
                + " · ".join(
                    f"[{kind}]({image_url(r,kind)})"
                    for kind in ("Science", "Template", "Difference")
                )
                + ".",
                "",
            ]
        )
    lines.extend(
        [
            "The controls exercise one foreground-contamination hypothesis. This small purposive sample leaves sensitivity to other contaminants, population precision, and the yield of new objects unmeasured. The ranked objects are historical leads for checking.",
            "",
        ]
    )
    (folder / "report.md").write_text("\n".join(lines))
    return packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("replay")
    run.add_argument("--snapshot", type=Path, default=ROOT / "data/pilot")
    run.add_argument("--out", type=Path, default=ROOT / "results")
    new = commands.add_parser("fetch")
    new.add_argument("--start", type=date.fromisoformat, required=True)
    new.add_argument(
        "--stop", type=date.fromisoformat, required=True, help="exclusive UTC date"
    )
    new.add_argument("--limit", type=int, default=20)
    new.add_argument("--snapshot", type=Path, required=True)
    follow = commands.add_parser(
        "vsx", help="complete or verify the second catalog pass"
    )
    follow.add_argument("--snapshot", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "fetch":
        refresh(args.snapshot, args.start, args.stop, args.limit)
        reconcile_snapshot(args.snapshot, fetch=True)
        print(f"Snapshot saved to {args.snapshot}")
    elif args.command == "vsx":
        reconcile_snapshot(args.snapshot, fetch=True)
        print(f"VSX checks archived in {args.snapshot / 'vsx'}")
    else:
        p = outputs(reconcile_snapshot(args.snapshot), args.out)
        print(
            f"Gate passed: {p['gate_passed']}. Cases emitted: {len(p['cases'])}. Results: {args.out}"
        )
        if not p["gate_passed"]:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
