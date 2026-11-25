"""Check a moving-star alert's inferred reference date against the original ZTF coadd."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode

import pilot


def investigate(snapshot, name, folder, fetch):
    result = pilot.evaluate(snapshot)
    case = next(
        (
            r
            for r in result["objects"]
            if r["id"] == name and r["status"] == "moving-star alternative"
        ),
        None,
    )
    if not result["gate_passed"] or case is None or len(case["gaia_matches"]) != 1:
        raise ValueError(
            "select a gated, non-control object with one moving-star match"
        )
    match = case["gaia_matches"][0]
    sid = match["source_id"]
    _, raw = pilot.snapshot(snapshot)
    star = next(
        s
        for key, rows in raw.items()
        if key.startswith("gaia-")
        for s in rows
        if s["source_id"] == sid
    )
    return inspect_match(case, star, folder, fetch)


def inspect_match(case, star, folder, fetch):
    name = case["id"]
    match = case["gaia_matches"][0]
    sid = match["source_id"]
    if not sid.isdigit():
        raise ValueError("invalid Gaia source ID")
    query = (
        "SELECT TOP 10 b.main_id,b.otype,b.pmra,b.pmdec FROM basic AS b "
        f"JOIN ident AS i ON i.oidref=b.oid WHERE i.id='Gaia DR3 {sid}'"
    )
    url = (
        pilot.SIMBAD
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
        )
    )
    ref_url = (
        pilot.FINK
        + "/objects?"
        + urlencode(
            {
                "objectId": name,
                "columns": "i:candid,i:jd,i:field,i:rcid,i:fid,i:jdstartref,i:jdendref,i:nframesref",
            }
        )
    )
    manifest_file = folder / "manifest.json"
    folder.mkdir(exist_ok=True, parents=True)
    receipt = json.loads(manifest_file.read_text()) if manifest_file.exists() else {}

    def source(key, url):
        path = folder / (key + ".gz")
        if key not in receipt and fetch:
            data = pilot.curl(url, [])
            path.write_bytes(gzip.compress(data, mtime=0))
            receipt[key] = {
                "url": url,
                "sha256_uncompressed": hashlib.sha256(data).hexdigest(),
                "retrieved_utc": pilot.datetime.now(pilot.timezone.utc).isoformat(
                    timespec="seconds"
                ),
            }
        data = gzip.decompress(path.read_bytes())
        if (
            receipt[key]["url"] != url
            or hashlib.sha256(data).hexdigest() != receipt[key]["sha256_uncompressed"]
        ):
            raise ValueError("archived source mismatch")
        return data

    ref = next(
        r
        for r in json.loads(source("reference.json", ref_url))
        if str(r["i:candid"]) == case["candid"]
    )
    objects = pilot.parse_csv(source("identifier.csv", url))
    names = ", ".join(o["main_id"] for o in objects) or "unresolved"
    lo, hi = ref["i:jdstartref"], ref["i:jdendref"]
    x, y = pilot.offsets(case["ra"], case["dec"], float(star["ra"]), float(star["dec"]))
    vx, vy = float(star["pmra"]) / 1000, float(star["pmdec"]) / 1000
    inferred = float(star["ref_epoch"]) + (x * vx + y * vy) / (vx * vx + vy * vy)
    inferred_jd = 2451545 + (inferred - 2000) * 365.25
    verdict = "inside" if lo <= inferred_jd <= hi else "outside"
    interpretation = (
        "Together the sign, positional offset and coadd dates support subtraction at the star's former position."
        if verdict == "inside" and match["negative_rows"] == case["usable_rows"] and case["usable_rows"] > 0
        else "The interval and sign comparison leaves the proposed reference subtraction uncertain. Multi-epoch image registration can test it."
    )
    text = f"""# Reference-image test for {name}

The calibrated foreground screen flagged Gaia DR3 {sid} ({names}). The anomaly packet has candid {case['candid']} and epoch JD {case['jd']}. [SIMBAD identity query]({url}) · [Fink object]({pilot.PORTAL+name})

Projecting the alert position onto the straight Gaia track gives {inferred:.3f}, about {pilot.day_of(inferred_jd)}. This is a geometric estimate of when the star crossed the residual position. Its perpendicular offset is {match['distance_to_track_arcsec']}″.

The alert's {ref['i:nframesref']}-frame ZTF reference spans {pilot.day_of(lo)} through {pilot.day_of(hi)} (JD {lo}–{hi}, field {ref['i:field']}, readout channel {ref['i:rcid']}, filter {ref['i:fid']}). The inferred crossing falls **{verdict} the actual coadd interval**. [Fink reference metadata]({ref_url}) · [ZTF field definitions](https://zwickytransientfacility.github.io/ztf-avro-alert/schema.html)

{match['negative_rows']} of {case['usable_rows']} quality-selected detections through this alert are negative residuals. {interpretation} The effective epoch depends on the weights of individual reference exposures. A superposed variable remains possible. Multi-epoch image registration and the reference input exposures can discriminate between these explanations.

The estimate uses Gaia's catalog `pmra` (already multiplied by cos δ), `pmdec`, and reference epoch. It omits parallax and perspective acceleration. This targeted follow-up used an object outside the seven calibration controls. It does not independently measure the false-positive rate of the screen.
"""
    pilot.json_write(manifest_file, receipt)
    (folder / "report.md").write_text(text)
    return verdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("object_id")
    parser.add_argument("--snapshot", type=Path, default=pilot.ROOT / "data/pilot")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--fetch", action="store_true", help="archive missing public source responses"
    )
    args = parser.parse_args()
    print(investigate(args.snapshot, args.object_id, args.out, args.fetch))
