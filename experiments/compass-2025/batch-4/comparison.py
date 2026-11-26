"""Measure how the unclassified transient compares with a frozen panel of bright supernovae."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import io
from datetime import date, datetime, timezone
from xml.etree import ElementTree
import zipfile
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlencode

from onset import HERE, capture, compass_archive, pilot, inputs, group_forced
import vsx

FOLDER = HERE / "comparison-sources"
TYPES = {
    "SN Ia": "Ia",
    "SN II": "II",
    "SN IIP": "II",
    "SN Ib": "Ib/c",
    "SN Ic": "Ib/c",
    "SN Ib/c": "Ib/c",
    "SN IIb": "IIb",
    "SN IIn": "IIn",
}
URL = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?" + urlencode(
    dict(
        f="s",
        subsample="sn",
        quality="y",
        startpeakdate="2024-01-01",
        endpeakdate="2024-12-31",
        endpeakmag="17.5",
        showpeakdate="y",
        showpeakmag="y",
        showz="y",
        format="csv",
    )
)


DISTANT_URL = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?" + urlencode(
    dict(
        f="s",
        subsample="trans",
        quality="y",
        startpeakdate="2018-01-01",
        endpeakdate="2024-12-31",
        startz="0.25",
        endz="0.80",
        endpeakmag="20.5",
        showpeakdate="y",
        showpeakmag="y",
        showz="y",
        format="csv",
    )
)
DISTANT_TYPES = {**TYPES, "SLSN-I": "SLSN-I", "SLSN-II": "SLSN-II", "TDE": "TDE"}


def catalog(fetch=False, distant=False):
    name, url = ("bts-distant", DISTANT_URL) if distant else ("bts", URL)
    return capture(
        FOLDER,
        name,
        url,
        {"ZTFID", "IAUID", "type", "RA", "Dec", "peakt", "peakmag", "redshift"},
        fetch,
    )


def panel(rows, distant=False):
    types = DISTANT_TYPES if distant else TYPES
    selected = defaultdict(list)
    seen = set()
    for r in rows:
        name = r["ZTFID"]
        if (
            r["type"] not in types
            or not r["IAUID"].startswith(("SN", "AT"))
            or not r["peakt"].strip()
        ):
            continue
        if name in seen:
            raise ValueError("Repeated control ID")
        t = float(r["peakt"]) + 2458000
        if (
            not ((2018 if distant else 2024) <= pilot.day_of(t).year <= 2024)
            or float(r["peakmag"]) > (20.5 if distant else 17.5)
            or (distant and not 0.25 <= float(r["redshift"]) <= 0.80)
        ):
            raise ValueError("BTS selection inconsistent with request")
        seen.add(name)
        selected[types[r["type"]]].append(r)
    return [
        r
        for group in sorted(selected)
        for r in sorted(
            selected[group], key=lambda r: hashlib.sha256(r["ZTFID"].encode()).digest()
        )[:24]
    ]


def registry(controls, fetch):
    names = {r["IAUID"][2:] for r in controls}
    url = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
    path = FOLDER / "tns.csv.gz"
    receipt = FOLDER / "tns.json"
    if fetch and not receipt.exists():
        raw = pilot.curl(
            url,
            [
                "-A",
                'tns_marker{"tns_id":0,"type":"guest","name":"Codex astronomy archive query"}',
            ],
        )
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            with z.open("tns_public_objects.csv") as stream:
                rows = csv.reader(
                    io.TextIOWrapper(stream, encoding="utf-8-sig", newline="")
                )
                stage = next(rows)
                fields = next(rows)
                k = fields.index("name")
                dest = io.StringIO(newline="")
                writer = csv.writer(dest, lineterminator="\n")
                writer.writerow(fields)
                for row in rows:
                    if len(row) > k and row[k] in names:
                        writer.writerow(row)
        cut = dest.getvalue().encode()
        path.write_bytes(gzip.compress(cut, mtime=0))
        pilot.json_write(
            receipt,
            dict(
                url=url,
                stage=stage,
                requested=sorted(names),
                sha256=hashlib.sha256(cut).hexdigest(),
                full_zip_sha256=hashlib.sha256(raw).hexdigest(),
                retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ),
        )
    cut = gzip.decompress(path.read_bytes())
    meta = json.loads(receipt.read_text())
    if (
        meta["url"] != url
        or sorted(names) != meta["requested"]
        or hashlib.sha256(cut).hexdigest() != meta["sha256"]
    ):
        raise ValueError("TNS control snapshot changed")
    return {r["name"]: r for r in csv.DictReader(io.StringIO(cut.decode()))}


def sexagesimal(text):
    a, b, c = map(float, text.split(":"))
    return (-1 if text.startswith("-") else 1) * (abs(a) + b / 60 + c / 3600)


def audit(r, records):
    t = records.get(r["IAUID"][2:])
    if not t:
        return {"ok": False, "reason": "TNS row missing"}
    ra, dec = sexagesimal(r["RA"]) * 15, sexagesimal(r["Dec"])
    separation = math.hypot(
        *pilot.offsets(ra, dec, float(t["ra"]), float(t["declination"]))
    )
    gap = (
        pilot.day_of(float(r["peakt"]) + 2458000)
        - date.fromisoformat(t["discoverydate"][:10])
    ).days
    same = DISTANT_TYPES.get(t["type"]) == DISTANT_TYPES[r["type"]]
    return {
        "ok": separation <= 3 and abs(gap) <= 180 and same,
        "separation_arcsec": separation,
        "peak_minus_discovery_days": gap,
        "tns_type": t["type"],
    }


def control_rows(name, fetch):
    url = f"https://api.alerce.online/ztf/v1/objects/{name}/detections"
    return compass_archive.capture_list(FOLDER / "alerts", name, url, fetch)


def control_curve(raw):
    exposures = defaultdict(list)
    for r in raw:
        if r["fid"] in (1, 2):
            exposures[r["pid"]].append(r)
    dates = defaultdict(list)
    conflicts = 0
    for rows in exposures.values():
        r = max(rows, key=lambda r: (bool(r.get("has_stamp")), str(r["candid"])))
        if any(
            abs(x["mjd"] - r["mjd"]) > 1 / 86400
            or x["fid"] != r["fid"]
            or x["isdiffpos"] != r["isdiffpos"]
            or abs(x["magpsf"] - r["magpsf"]) > 0.01
            for x in rows
        ):
            conflicts += 1
            continue
        if (
            r["isdiffpos"] == 1
            and 0 < r["sigmapsf"] <= 0.2
            and math.isfinite(r["magpsf"])
        ):
            dates[int(r["mjd"]), r["fid"]].append(r["magpsf"])
    return {k: median(v) for k, v in dates.items()}, conflicts


def windows(curve, near):
    reds = [(m, d) for (d, f), m in curve.items() if f == 2 and abs(d - near) <= 15]
    if not reds:
        return None, {}
    mag0, peak = min(reds)
    red = {
        d - peak: m - mag0
        for (d, f), m in curve.items()
        if f == 2 and 25 <= d - peak <= 40
    }
    col = {
        d - peak: m - curve[d, 2]
        for (d, f), m in curve.items()
        if f == 1 and (d, 2) in curve
    }
    return peak, {
        "red": red,
        "early_color": {d: v for d, v in col.items() if -10 <= d <= 2},
        "late_color": {d: v for d, v in col.items() if 20 <= d <= 35},
    }


def match(dates, target):
    # Dynamic programming maximizes distinct matches before minimizing total separation.
    a, b = sorted(target), sorted(dates)
    best = {(len(a), j): (0, 0, ()) for j in range(len(b) + 1)}
    best.update({(i, len(b)): (0, 0, ()) for i in range(len(a) + 1)})
    for i in reversed(range(len(a))):
        for j in reversed(range(len(b))):
            opts = [best[i + 1, j], best[i, j + 1]]
            dist = abs(a[i] - b[j])
            if dist <= 2:
                n, neg, chosen = best[i + 1, j + 1]
                opts.append((n + 1, neg - dist, (b[j],) + chosen))
            best[i, j] = max(opts)
    return best[0, 0][2]


def features(w, anchor=None):
    if not w:
        return None
    med = {}
    used = {}
    for key, values in w.items():
        dates = sorted(values) if anchor is None else sorted(match(values, anchor[key]))
        if len(dates) < 2 or (anchor is not None and len(dates) * 2 < len(anchor[key])):
            return None
        med[key] = median(values[d] for d in dates)
        used[key] = dates
    return dict(
        decline=med["red"],
        reddening=med["late_color"] - med["early_color"],
        early_gr=med["early_color"],
        late_gr=med["late_color"],
        offsets=used,
    )


def dist(a, b):
    return math.hypot(a["decline"] - b["decline"], a["reddening"] - b["reddening"])


def compare(rows, target, key):
    good = [r for r in rows if r[key] is not None]
    groups = sorted({r["family"] for r in rows})

    def neighbors(q, others):
        return sorted(others, key=lambda r: (dist(q, r[key]), r["id"]))

    nearest = [
        {
            "id": r["id"],
            "tns": r["tns"],
            "family": r["family"],
            "distance_mag": dist(target, r[key]),
        }
        for r in neighbors(target, good)[:3]
    ]
    confusion = {g: dict.fromkeys(groups, 0) for g in groups}
    for x in good:
        if len(good) < 2:
            break
        guess = neighbors(x[key], [r for r in good if r["id"] != x["id"]])[0]["family"]
        confusion[x["family"]][guess] += 1
    stats = {}
    for g in groups:
        v = [r[key] for r in good if r["family"] == g]
        stats[g] = {"usable": len(v), "selected": sum(r["family"] == g for r in rows)}
        for col in ("decline", "reddening"):
            if v:
                vals = sorted(r[col] for r in v)
                stats[g][col] = {
                    "min": min(vals),
                    "median": median(vals),
                    "max": max(vals),
                }
    return {"families": stats, "nearest": nearest, "leave_one_out": confusion}


def tde_catalog(fetch):
    url = (
        vsx.SERVICE
        + "?"
        + urlencode(
            {
                "-source": "J/ApJ/942/9/table1,J/ApJ/942/9/table5",
                "-out.all": 1,
                "-out.max": 2000,
            }
        )
    )
    path = FOLDER / "tde.xml.gz"
    rec = FOLDER / "tde.json"
    if fetch and not rec.exists():
        raw = pilot.curl(url, [])
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            rec,
            dict(
                url=url,
                sha256=hashlib.sha256(raw).hexdigest(),
                retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ),
        )
    raw = gzip.decompress(path.read_bytes())
    meta = json.loads(rec.read_text())
    if meta["url"] != url or meta["sha256"] != hashlib.sha256(raw).hexdigest():
        raise ValueError("TDE catalog changed")
    tables = {}
    for tb in ElementTree.fromstring(raw).findall(".//{*}TABLE"):
        fields = [f.get("name") for f in tb.findall("{*}FIELD")]
        tables[tb.get("name").split("/")[-1]] = [
            dict(zip(fields, [td.text for td in tr], strict=True))
            for tr in tb.findall(".//{*}TR")
        ]
    if len(tables["table1"]) != 30 or len(tables["table5"]) >= 2000:
        raise ValueError("Unexpected TDE catalog coverage")
    out = []
    for src in tables["table1"]:
        name = src["ZTF"]
        spec = [
            r
            for r in tables["table5"]
            if r["ZTF"] == name and r["Phase"] and abs(float(r["Phase"])) <= 20
        ]
        row = dict(
            id=name,
            tns=src["Name"],
            type=src["Spec"],
            family="TDE",
            redshift=src["zspec"],
            reference="Hammerstein et al. 2023",
        )
        if spec and src["Spec"] in ("TDE-H", "TDE-He", "TDE-H+He"):
            obs = min(spec, key=lambda r: abs(float(r["Phase"])))
            row["peak_estimate"] = (
                pilot.utc_jd(date.fromisoformat(obs["Date"]))
                - 2400000.5
                - float(obs["Phase"])
            )
        out.append(row)
    return out


def run(fetch=False):
    raw, _, _ = inputs(False)
    curve = {
        (r["day"], r["fid"]): r["mag"]
        for r in group_forced(raw)
        if r["mag"] and r["magerr"] <= 0.2
    }
    peak, target_windows = windows(curve, 60790)
    target = features(target_windows)
    if target is None:
        raise ValueError("Target lacks registered windows")
    panels = {
        "2024_bright": panel(catalog(fetch)),
        "distant": panel(catalog(fetch, True), True),
    }
    tns = registry(panels["2024_bright"] + panels["distant"], fetch)

    def one_bts(r):
        a = audit(r, tns)
        if r["IAUID"] == "AT2023vto":
            a.update(
                ok=False,
                reason="Published He TDE spectrum conflicts with catalog SLSN-II label",
            )
        return {
            "id": r["ZTFID"],
            "tns": r["IAUID"],
            "type": r["type"],
            "family": DISTANT_TYPES[r["type"]],
            "redshift": r["redshift"],
            "audit": a,
            "peak_estimate": float(r["peakt"]) + 57999.5,
        }

    panels = {key: [one_bts(r) for r in values] for key, values in panels.items()}
    panels["nuclear_TDE"] = tde_catalog(fetch)

    def measure(row):
        if not row.get("audit", {"ok": True})["ok"] or "peak_estimate" not in row:
            return {**row, "all": None, "matched": None}
        c, bad = control_curve(control_rows(row["id"], fetch))
        ref, ws = windows(c, row["peak_estimate"])
        return {
            **row,
            "red_peak_mjd_day": ref,
            "conflicting_exposures": bad,
            "all": features(ws),
            "matched": features(ws, target_windows),
        }

    answer = {}
    with ThreadPoolExecutor(max_workers=3 if fetch else 1) as pool:
        for name, listed in panels.items():
            rows = list(pool.map(measure, listed))
            answer[name] = {
                "controls": rows,
                "class_comparison": {
                    key: compare(rows, target, key) for key in ("all", "matched")
                },
            }
    return {"target_red_peak_mjd_day": peak, "target": target, "panels": answer}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "results/comparison.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
