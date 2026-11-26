"""Check the single 2012 y-band forced measurement against field stars on the same warp."""

import argparse
import json
import math
from statistics import median
from urllib.parse import urlencode

from check import HERE, YEAR, NAME, compass_archive, compass_direct, pilot, query

WARP = "55762501"


def run(fetch=False):
    anchors, _, _ = compass_archive.load(YEAR / "snapshot")
    anchor = next(r for r in anchors if pilot.oid(r) == NAME)
    ra, dec = anchor["i:ra"], anchor["i:dec"]
    catalog = query(
        "stack",
        "field-stack",
        {"ra": ra, "dec": dec, "radius": 60 / 3600},
        ["objID", "raStack", "decStack", "primaryDetection", "yPSFMag"],
        fetch,
    )
    controls = []
    for star in catalog:
        offset = math.hypot(
            *pilot.offsets(float(star["raStack"]), float(star["decStack"]), ra, dec)
        )
        if (
            4 < offset <= 60
            and 18 < float(star["yPSFMag"]) < 23
            and star["primaryDetection"] == "1"
        ):
            controls.append((offset, star["objID"]))
    controls = sorted(controls)[:10]
    (source,) = json.loads((HERE / "results.json").read_text())[
        "matches_within_1.5arcsec"
    ]
    selected = [source] + [k for d, k in controls]
    columns = [
        "objid",
        "uniquepspsfwid",
        "forcedwarpid",
        "obstime",
        "fpsfflux",
        "fpsffluxerr",
        "fpsfqfperfect",
        "fpsfchisq",
        "fsky",
        "fapflux",
        "fkronflux",
        "finfoflag3",
    ]
    adql = (
        "SELECT TOP 1000 "
        + ",".join(columns)
        + " FROM ps1_dr2.forced_warp_measurement "
        + "WHERE filterid=5 AND objid IN ("
        + ",".join(selected)
        + ") ORDER BY objid,obstime,uniquepspsfwid"
    )
    url = "https://mast.stsci.edu/vo-tap/api/v0.1/mast_catalogs/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}
    )
    rows = compass_direct.capture(
        HERE / "sources", "field-forced", url, set(columns), fetch
    )
    if (
        len(rows) >= 1000
        or len({r["uniquepspsfwid"] for r in rows}) != len(rows)
        or any(r["objid"] not in selected for r in rows)
    ):
        raise ValueError("Unexpected forced field sample")
    out = []
    for offset, obj in [(None, source)] + controls:
        quality = [
            r
            for r in rows
            if r["objid"] == obj
            and float(r["fpsfqfperfect"]) >= 0.85
            and float(r["fpsffluxerr"]) > 0
            and float(r["fpsfflux"]) != -999
        ]
        event = [r for r in quality if r["forcedwarpid"] == WARP]
        other = [r for r in quality if r["forcedwarpid"] != WARP]
        record = {
            "id": obj,
            "target": obj == source,
            "offset_arcsec": round(offset, 2) if offset else 0,
            "other_exposures": len(other),
            "event_exposures": len(event),
        }
        if len(event) == 1 and len(other) >= 3:
            (r,) = event
            residual = float(r["fpsfflux"]) - median(
                float(v["fpsfflux"]) for v in other
            )
            record.update(
                flux_excess_uJy=round(residual * 1e6, 2),
                excess_over_event_error=round(residual / float(r["fpsffluxerr"]), 2),
                event_sky_uJy_arcsec2=round(float(r["fsky"]) * 1e6, 2),
                typical_sky_uJy_arcsec2=round(
                    median(float(v["fsky"]) for v in other) * 1e6, 2
                ),
                event_reduced_chisq=round(float(r["fpsfchisq"]), 2),
                typical_reduced_chisq=round(
                    median(float(v["fpsfchisq"]) for v in other), 2
                ),
                survey_PSF_IRLS_mask=bool(int(r["finfoflag3"]) & 0x01000000),
            )
        out.append(record)
    return {"forced_warp": WARP, "sources": out}


def figure(result):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.2, 2.8), layout="constrained")
    for target, color, symbol, name in [
        (False, "#475569", "o", "Field stars"),
        (True, "#a53d15", "D", "ASASSN-25ba"),
    ]:
        rows = [
            r
            for r in result["sources"]
            if r["target"] == target and "flux_excess_uJy" in r
        ]
        ax.scatter(
            [r["offset_arcsec"] for r in rows],
            [r["excess_over_event_error"] for r in rows],
            s=40,
            marker=symbol,
            color=color,
            label=name,
        )
    ax.axhline(0, color="0.4", linewidth=0.9)
    ax.set(
        xlabel="Distance from ASASSN-25ba (arcsec)",
        ylabel="Flux excess / error",
        ylim=(-0.6, 11.5),
        title="Pan-STARRS y, 2012 October 19",
    )
    ax.legend(loc="lower right", frameon=False)
    fig.savefig(HERE / "field.png", dpi=145)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()
    result = run(args.fetch)
    file = HERE / "field.json"
    if args.fetch:
        pilot.json_write(file, result)
    elif json.loads(file.read_text()) != result:
        raise ValueError("Forced field comparison changed")
    if args.plot:
        figure(result)
    print(json.dumps(result, indent=2))
