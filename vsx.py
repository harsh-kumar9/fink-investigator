"""Resolve additional variable-star counterparts for the gated first-pass cohort."""

import gzip
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

SERVICE = "https://vizier.cfa.harvard.edu/viz-bin/votable"


def query_url(case):
    return (
        SERVICE
        + "?"
        + urlencode(
            {
                "-source": "B/vsx/vsx",
                "-c": f"{case['ra']} {case['dec']:+}",
                "-c.rs": "5",
                "-out": "OID,Name,Type,RAJ2000,DEJ2000",
                "-out.add": "_r",
                "-out.max": "50",
            }
        )
    )


def parse(raw, case):
    root = ET.fromstring(raw)
    table = root.find(".//{*}TABLE[@name='B/vsx/vsx']")
    if table is None or any(
        i.get("value", "").upper() in {"ERROR", "OVERFLOW"}
        or i.get("name", "").lower() == "error"
        for i in root.findall(".//{*}INFO")
    ):
        raise ValueError("VSX lookup failed or overflowed")
    cols = [f.attrib["name"] for f in table.findall("{*}FIELD")]
    if not {"_r", "OID", "Name", "Type", "RAJ2000", "DEJ2000"}.issubset(cols):
        raise ValueError("missing expected VSX fields")
    matches = []
    for tr in table.findall(".//{*}TR"):
        row = dict(
            zip(
                cols,
                [(td.text or "").strip() for td in tr.findall("{*}TD")],
                strict=True,
            )
        )
        ra, dec = float(row["RAJ2000"]), float(row["DEJ2000"])
        dx = ((ra - case["ra"] + 180) % 360 - 180) * 3600 * math.cos(math.radians(dec))
        dy = (dec - case["dec"]) * 3600
        sep = float(row["_r"])
        if (
            not row["OID"].isdigit()
            or not 0 <= sep <= 5.01
            or abs(math.hypot(dx, dy) - sep) > 0.10
        ):
            raise ValueError("VSX response has inconsistent coordinates or identifier")
        matches.append(
            {
                "oid": row["OID"],
                "name": row["Name"],
                "type": row["Type"],
                "ra": ra,
                "dec": dec,
                "separation_arcsec": sep,
                "source": "https://vsx.aavso.org/index.php?view=detail.top&oid="
                + row["OID"],
            }
        )
    if len(matches) >= 50:
        raise ValueError("VSX cone response reached row limit")
    return sorted(matches, key=lambda m: m["separation_arcsec"])


def reconcile(folder, result, fetch=None):
    cases = [
        r
        for r in result["objects"]
        if r["status"] == "review" and result["gate_passed"]
    ]
    archive = folder / "vsx"
    if fetch:
        archive.mkdir(exist_ok=True)
    manifest_path = archive / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for r in cases:
        key = r["id"]
        url = query_url(r)
        path = archive / (key + ".xml.gz")
        if key not in manifest and fetch:
            raw = fetch(url)
            parse(raw, r)
            path.write_bytes(gzip.compress(raw, mtime=0))
            manifest[key] = {
                "url": url,
                "candid": r["candid"],
                "method": "GET",
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            }
            tmp = archive / "manifest.tmp"
            tmp.write_text(json.dumps(manifest, indent=2) + "\n")
            os.replace(tmp, manifest_path)
        if key not in manifest:
            raise ValueError(
                "missing VSX snapshot; run the vsx command on this snapshot"
            )
        raw = gzip.decompress(path.read_bytes())
        m = manifest[key]
        if (
            m["url"] != url
            or m["candid"] != r["candid"]
            or hashlib.sha256(raw).hexdigest() != m["sha256_uncompressed"]
        ):
            raise ValueError(f"VSX query or checksum mismatch for {key}")
        matches = parse(raw, r)
        r["vsx_matches"] = matches
        r["vsx_query"] = url
        if len(matches) == 1 and matches[0]["separation_arcsec"] <= 1.0:
            r["adjudication"] = "closed as a new source, retained for alias review"
        elif matches:
            r["adjudication"] = "inspect possible VSX counterpart"
        else:
            r["adjudication"] = "open for image and catalog checks"
    return result
