"""Compare MGAB-V488 with nearby stars on its December 2013 PS1 exposures."""

import argparse
import csv
import io
import json
import math
from statistics import median
from urllib.parse import urlencode

from ps1 import DEC, RA, archive, sep, table

DAYS = [56632, 56638]


def position(row):
    return float(row["raMean"]), float(row["decMean"])


def forced_measurements(obj, fetch):
    cols = [
        "objid",
        "uniquepspsfwid",
        "forcedwarpid",
        "detectid",
        "filterid",
        "obstime",
        "fpsfflux",
        "fpsffluxerr",
        "fpsfqfperfect",
        "fapflux",
        "fpsfchisq",
        "finfoflag",
        "finfoflag2",
        "finfoflag3",
    ]
    query = (
        "SELECT TOP 1000 "
        + ",".join(cols)
        + " FROM ps1_dr2.forced_warp_measurement WHERE objid="
        + obj
        + " AND ((obstime BETWEEN 56632 AND 56633) OR (obstime BETWEEN 56638 AND 56639)) AND filterid IN (2,3) ORDER BY obstime,uniquepspsfwid"
    )
    url = "https://mast.stsci.edu/vo-tap/api/v0.1/mast_catalogs/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    reader = csv.DictReader(io.StringIO(archive.source("ps1-forced.csv", url, fetch).decode()))
    rows = list(reader)
    if (
        reader.fieldnames != cols
        or len({r["uniquepspsfwid"] for r in rows}) != len(rows)
        or any(r["objid"] != obj for r in rows)
        or len(rows) >= 1000
    ):
        raise ValueError("Unexpected PS1 forced measurement table")
    result = []
    for r in rows:
        flux = float(r["fpsfflux"])
        err = float(r["fpsffluxerr"])
        ap = float(r["fapflux"])
        q = float(r["fpsfqfperfect"])
        if flux <= 0 or err <= 0:
            continue
        result.append(
            {
                "mjd": float(r["obstime"]),
                "band": "grizy"[int(r["filterid"]) - 1],
                "warp": r["forcedwarpid"],
                "mag": round(-2.5 * math.log10(flux / 3631), 4),
                "magerr": round(2.5 * err / flux / math.log(10), 4),
                "aperture_mag": round(-2.5 * math.log10(ap / 3631), 4) if ap > 0 else None,
                "psf_fraction_unmasked": q,
                "passes_strict": flux / err >= 10 and q >= 0.90,
                "flags": [int(r[k]) for k in ("finfoflag", "finfoflag2", "finfoflag3")],
            }
        )
    return result


def run(fetch=False):
    means = table(
        "mean",
        "field",
        {"ra": RA, "dec": DEC, "radius": 90 / 3600},
        ["objID", "raMean", "decMean", "nDetections", "rMeanPSFMag", "rMeanKronMag"],
        fetch,
    )
    nearby = [r for r in means if int(r["nDetections"]) >= 3]
    (target,) = [r for r in nearby if sep(*position(r)) < 1]
    eligible = []
    for r in nearby:
        rr = float(r["rMeanPSFMag"])
        ra, dec = position(r)
        if (
            r == target
            or int(r["nDetections"]) < 25
            or not 17 <= rr <= 20
            or abs(rr - float(r["rMeanKronMag"])) >= 0.08
        ):
            continue
        if any(k != r and math.hypot(*archive.pilot.offsets(ra, dec, *position(k))) < 3 for k in nearby):
            continue
        eligible.append((sep(ra, dec), r["objID"]))
    controls = sorted(eligible)[:10]
    ids = [target["objID"]] + [k for d, k in controls]
    cols = [
        "objid",
        "imageid",
        "detectid",
        "filterid",
        "obstime",
        "ra",
        "dec",
        "psfflux",
        "psffluxerr",
        "psfqfperfect",
        "psfmajorfwhm",
        "psfminorfwhm",
        "psfchisq",
        "apflux",
        "infoflag",
        "infoflag2",
        "infoflag3",
    ]
    query = (
        "SELECT TOP 1000 "
        + ",".join(cols)
        + " FROM ps1_dr2.detection WHERE objid IN ("
        + ",".join(ids)
        + ") AND ((obstime BETWEEN 56632 AND 56633) OR (obstime BETWEEN 56638 AND 56639)) AND filterid IN (2,3) ORDER BY objid,obstime,detectid"
    )
    url = "https://mast.stsci.edu/vo-tap/api/v0.1/mast_catalogs/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    raw = archive.source("ps1-field-epochs.csv", url, fetch)
    reader = csv.DictReader(io.StringIO(raw.decode()))
    rows = list(reader)
    if (
        reader.fieldnames != cols
        or len(rows) >= 1000
        or len({r["detectid"] for r in rows}) != len(rows)
        or any(r["objid"] not in ids for r in rows)
    ):
        raise ValueError("Unexpected field epoch response")
    positions = {r["objID"]: position(r) for r in nearby}

    def mag(r):
        return -2.5 * math.log10(float(r["psfflux"]) / 3631)

    def offset(r):
        return math.hypot(*archive.pilot.offsets(float(r["ra"]), float(r["dec"]), *positions[r["objid"]]))

    def quality(r):
        return (
            float(r["psffluxerr"]) > 0
            and float(r["psfflux"]) / float(r["psffluxerr"]) >= 10
            and float(r["psfqfperfect"]) >= 0.90
            and offset(r) <= 1.2
        )

    target_rows = [r for r in rows if r["objid"] == target["objID"]]
    images = {r["imageid"]: float(r["obstime"]) for r in target_rows if quality(r)}
    accepted = [
        r
        for r in rows
        if quality(r)
        and r["imageid"] in images
        and abs(float(r["obstime"]) - images[r["imageid"]]) * 86400 <= 5
    ]
    source_summary = []
    for r in target_rows:
        source_summary.append(
            {
                "imageid": r["imageid"],
                "detectid": r["detectid"],
                "mjd": float(r["obstime"]),
                "band": "grizy"[int(r["filterid"]) - 1],
                "offset_from_ps1_mean_arcsec": round(offset(r), 3),
                "offset_from_ztf_anchor_arcsec": round(sep(float(r["ra"]), float(r["dec"])), 3),
                "mag": round(mag(r), 4),
                "magerr": round(2.5 * float(r["psffluxerr"]) / float(r["psfflux"]) / math.log(10), 4),
                "aperture_mag": (
                    round(-2.5 * math.log10(float(r["apflux"]) / 3631), 4) if float(r["apflux"]) > 0 else None
                ),
                "psf_fraction_unmasked": round(float(r["psfqfperfect"]), 3),
                "psf_fwhm_arcsec": [round(float(r[k]), 3) for k in ("psfmajorfwhm", "psfminorfwhm")],
                "psf_reduced_chisq": round(float(r["psfchisq"]), 3),
                "flags": [int(r[k]) for k in ("infoflag", "infoflag2", "infoflag3")],
                "defect_flag": bool(int(r["infoflag"]) & 0x800),
                "masked_for_mean_psf": bool(int(r["infoflag3"]) & 0x01000000),
                "burntool_flag": bool(int(r["infoflag2"]) & 0x20),
                "passes_strict": quality(r),
            }
        )
    results = {}
    for fid, band in [("2", "r"), ("3", "i")]:
        changes = {}
        for obj in ids:
            mags = [
                [
                    mag(r)
                    for r in accepted
                    if r["objid"] == obj and r["filterid"] == fid and math.floor(float(r["obstime"])) == d
                ]
                for d in DAYS
            ]
            if all(mags):
                changes[obj] = round(median(mags[1]) - median(mags[0]), 4)
        dm = changes.pop(target["objID"])
        res = list(changes.values())
        center = median(res) if res else None
        results[band] = {
            "target_change_mag": dm,
            "control_changes_mag": changes,
            "median_control_mag": round(center, 5) if res else None,
            "control_MAD_mag": round(median(abs(delta - center) for delta in res), 4) if res else None,
            "target_minus_comparison_median_mag": round(dm - center, 4) if res else None,
        }
    return {
        "mean_rows_90arcsec": len(means),
        "eligible_controls": len(eligible),
        "selected_controls": [{"id": obj, "distance_arcsec": round(d, 2)} for d, obj in controls],
        "target_exposures": source_summary,
        "same_image_changes": results,
        "forced_warp_photometry": forced_measurements(target["objID"], fetch),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
