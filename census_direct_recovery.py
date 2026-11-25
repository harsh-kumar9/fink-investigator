"""Resolve timed-out IRSA cones into DR24 object IDs and request the same public data by ID."""

import argparse
import csv
import io
from datetime import datetime, timezone
import gzip
import hashlib
import json
from urllib.parse import urlencode

import census
import pilot

SERVICE = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"


def recovered(case, directory, fetch=False):
    """Replay the full bounded cone first, then the light curves of every returned series."""
    if fetch:
        directory.mkdir(exist_ok=True, parents=True)
    query = (
        "SELECT TOP 100 oid,ra,dec,nobs,fid FROM ztf_objects_dr24 "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{case['ra']},{case['dec']},0.0004166667))=1 "
        "AND fid=2 AND nobs>=30"
    )
    tap_url = (
        pilot.IRSA
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "CSV", "QUERY": query}
        )
    )
    receipt = directory / (case["id"] + ".json")
    meta = json.loads(receipt.read_text()) if receipt.exists() else {}

    def capture(suffix, url):
        path = directory / (case["id"] + "." + suffix + ".csv.gz")
        if suffix not in meta and fetch:
            raw = pilot.curl(url, ["--max-time", "480", "--retry", "0"])
            if raw.lstrip().startswith(b"<"):
                raise ValueError(f"IRSA table error: {raw[:300]!r}")
            path.write_bytes(gzip.compress(raw, mtime=0))
            meta[suffix] = {
                "url": url,
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
            }
            pilot.json_write(receipt, meta)
        raw = gzip.decompress(path.read_bytes())
        if (
            meta[suffix]["url"] != url
            or meta[suffix]["sha256_uncompressed"] != hashlib.sha256(raw).hexdigest()
        ):
            raise ValueError("IRSA ID recovery request or checksum mismatch")
        return raw

    entries = list(csv.DictReader(io.StringIO(capture("cone", tap_url).decode())))
    if not 0 < len(entries) < 100 or any(
        not r["oid"].isdigit() or int(r["fid"]) != 2 or int(r["nobs"]) < 30
        for r in entries
    ):
        raise ValueError("unexpected IRSA direct lightcurve IDs")
    ids = {r["oid"] for r in entries}
    options = [("ID", oid) for oid in sorted(ids)]
    options.extend(
        [
            ("COLLECTION", "ztf_dr24"),
            ("BAD_CATFLAGS_MASK", "65535"),
            ("TIME", f"-Inf {case['jd']-2400000.5:.8f}"),
            ("FORMAT", "csv"),
        ]
    )
    url = SERVICE + "?" + urlencode(options)
    raw = capture("photometry", url)
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if not {"oid", "expid", "mjd", "mag", "magerr", "catflags", "filtercode"}.issubset(
        reader.fieldnames or []
    ):
        raise ValueError("Invalid IRSA recovered light-curve columns")
    if {r["oid"] for r in reader} != ids:
        raise ValueError("IRSA photometry did not cover the exact recovery ID set")
    return raw, url


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("id")
    args = parser.parse_args()
    _, cases = census.assess(False)
    case = next(
        r
        for r in cases
        if r["id"] == args.id and r.get("adjudication") and not r["vsx_matches"]
    )
    raw, url = recovered(case, census.CENSUS / "sources/direct/recovery", True)
    print(case["id"], len(raw), url)
