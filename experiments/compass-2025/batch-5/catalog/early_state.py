"""Check the 2013 MGAB-V306 state with PS1 fixed-position photometry and field stars."""

import argparse
from collections import defaultdict
import csv
import io
import json
import math
from statistics import median
from urllib.parse import urlencode

from ps1 import STARS, archive, fetch_table

NAME = "MGAB-V306"
DAY = 56599


def run(fetch):
    ra, dec = STARS[NAME][:2]
    cat = fetch_table(
        "mean",
        "field-" + NAME,
        {"ra": ra, "dec": dec, "radius": 90 / 3600},
        ["objID", "raMean", "decMean", "nDetections", "rMeanPSFMag", "rMeanKronMag"],
        fetch,
    )
    known = [r for r in cat if int(r["nDetections"]) >= 3]

    def pos(r):
        return float(r["raMean"]), float(r["decMean"])

    def distance(a, b):
        return math.hypot(*archive.pilot.offsets(*a, *b))

    (target,) = [r["objID"] for r in known if distance(pos(r), (ra, dec)) <= 1]
    eligible = [
        (distance(pos(r), (ra, dec)), r["objID"])
        for r in known
        if r["objID"] != target
        and int(r["nDetections"]) >= 25
        and 17 <= float(r["rMeanPSFMag"]) <= 20
        and abs(float(r["rMeanPSFMag"]) - float(r["rMeanKronMag"])) < 0.08
        and not any(s != r and distance(pos(s), pos(r)) < 3 for s in known)
    ]
    controls = sorted(eligible)[:10]
    ids = [target] + [obj for d, obj in controls]
    cols = [
        "objid",
        "uniquepspsfwid",
        "forcedwarpid",
        "filterid",
        "obstime",
        "fpsfflux",
        "fpsffluxerr",
        "fpsfqfperfect",
        "fpsfchisq",
        "fapflux",
        "finfoflag",
        "finfoflag2",
        "finfoflag3",
    ]
    query = (
        "SELECT TOP 4000 "
        + ",".join(cols)
        + " FROM ps1_dr2.forced_warp_measurement WHERE filterid IN (4,5) AND objid IN ("
        + ",".join(ids)
        + ") ORDER BY objid,obstime,uniquepspsfwid"
    )
    url = "https://mast.stsci.edu/vo-tap/api/v0.1/mast_catalogs/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    raw = archive.source("field-forced.csv", url, fetch)
    reader = csv.DictReader(io.StringIO(raw.decode()))
    rows = list(reader)
    if (
        reader.fieldnames != cols
        or len(rows) >= 4000
        or len({r["uniquepspsfwid"] for r in rows}) != len(rows)
        or any(r["objid"] not in ids for r in rows)
    ):
        raise ValueError("Unexpected repeated or truncated warp query")
    result = {
        "selected_controls": controls,
        "mean_catalog_90arcsec": len(cat),
        "target": target,
        "comparisons": {},
        "source_exposures": {},
        "raw_resampled_rows": len(rows),
    }

    def clean(r):
        return (
            float(r["fpsffluxerr"]) > 0 and float(r["fpsfflux"]) != -999 and float(r["fpsfqfperfect"]) >= 0.85
        )

    for f, b in [("4", "z"), ("5", "y")]:
        target_rows = [r for r in rows if r["objid"] == target and r["filterid"] == f and clean(r)]
        groups = defaultdict(list)
        for r in target_rows:
            groups[float(r["obstime"])].append(r)
        visits = []
        for t, warp in sorted(groups.items()):
            flux = median(float(r["fpsfflux"]) for r in warp)
            err = max(float(r["fpsffluxerr"]) for r in warp)
            snr = flux / err
            visits.append(
                {
                    "mjd": t,
                    "warps": [r["forcedwarpid"] for r in warp],
                    "flux_uJy": round(flux * 1e6, 3),
                    "error_uJy": round(err * 1e6, 3),
                    "mag_if_SNR5": round(-2.5 * math.log10(flux / 3631), 4) if snr >= 5 else None,
                    "depth_5sigma": round(-2.5 * math.log10(5 * err / 3631), 4),
                    "unmasked_fraction": min(float(r["fpsfqfperfect"]) for r in warp),
                }
            )
        result["source_exposures"][b] = visits
        target_warps = {r["forcedwarpid"] for r in target_rows}
        changes = {}
        for obj in ids:
            # Require the very same resampled image as an accepted measurement of the target.
            bytime = defaultdict(list)
            for r in rows:
                if not (
                    r["objid"] == obj
                    and r["filterid"] == f
                    and r["forcedwarpid"] in target_warps
                    and clean(r)
                ):
                    continue
                flux, err = float(r["fpsfflux"]), float(r["fpsffluxerr"])
                if flux / err >= 5:
                    bytime[float(r["obstime"])].append(-2.5 * math.log10(flux / 3631))
            byday = defaultdict(list)
            for t, mags in bytime.items():
                byday[math.floor(t)].append(median(mags))
            compare = [median(mags) for d, mags in byday.items() if d != DAY]
            if DAY in byday and len(compare) >= 3:
                changes[obj] = {
                    "event_minus_other_nights_median_mag": round(median(byday[DAY]) - median(compare), 4),
                    "other_nights": len(compare),
                }
        result["comparisons"][b] = changes
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
