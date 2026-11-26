"""Replay ZTF detections, signed forced photometry and upper limits for the short-history case."""

import argparse
from collections import Counter, defaultdict
from statistics import median
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
YEAR = HERE.parent
sys.path.insert(0, str(HERE.parents[2]))
import compass_archive
from compass_direct import capture
import pilot

NAME = "ZTF25aanlbls"
FOLDER = HERE / "onset-sources"
FIELDS = "oid,pid,mjd,fid,mag,e_mag,isdiffpos,field,rcid,rfid,procstatus,dubious,magzpsci,diffmaglim,sciinpseeing"


def forced_inputs(fetch, name, folder):
    query = f"SELECT TOP 10000 {FIELDS} FROM ztf.forced_photometry WHERE oid='{name}'"
    url = "https://tap.alerce.online/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    )
    forced = capture(folder, "forced", url, set(FIELDS.split(",")), fetch)
    if (
        len(forced) >= 10000
        or len({r["pid"] for r in forced}) != len(forced)
        or any(r["oid"] != name for r in forced)
    ):
        raise ValueError("Wrong or truncated forced series")
    return forced


def inputs(fetch, name=NAME, folder=FOLDER):
    forced = forced_inputs(fetch, name, folder)
    base = f"https://api.alerce.online/ztf/v1/objects/{name}"
    det = compass_archive.capture_list(
        folder, "detections", base + "/detections", fetch
    )
    lim = compass_archive.capture_list(
        folder, "limits", base + "/non_detections", fetch
    )
    return forced, det, lim


def ujy(mag):
    # Conventional AB zero point for microjanskys. The scale cancels in flux ratios.
    return 10 ** (0.4 * (23.9 - mag))


def group_forced(rows):
    groups = defaultdict(list)
    for r in rows:
        if (
            r["procstatus"] != "0"
            or r["dubious"] != "False"
            or not r["mag"]
            or not r["e_mag"]
        ):
            continue
        sign, mag, err = int(r["isdiffpos"]), float(r["mag"]), float(r["e_mag"])
        if sign not in (-1, 1) or not (math.isfinite(mag) and 0 < err < 99):
            continue
        t = float(r["mjd"])
        flux = sign * ujy(mag)
        sigma = abs(flux) * err * math.log(10) / 2.5
        groups[(int(t), int(r["fid"]))].append((t, flux, sigma, r["pid"], r["rfid"]))
    nights = []
    for (day, fid), epochs in sorted(groups.items()):
        sw = sum(1 / r[2] ** 2 for r in epochs)
        signal = sum(r[1] / r[2] ** 2 for r in epochs) / sw
        sigma = sw**-0.5
        mjd = sum(r[0] / r[2] ** 2 for r in epochs) / sw
        nights.append(
            dict(
                day=day,
                date=str(pilot.day_of(mjd + 2400000.5)),
                mjd=mjd,
                fid=fid,
                flux=signal,
                error=sigma,
                snr=signal / sigma,
                mag=23.9 - 2.5 * math.log10(signal) if signal > 0 else None,
                magerr=2.5 / math.log(10) * sigma / signal if signal > 0 else None,
                pids=[r[3] for r in epochs],
                references=[r[4] for r in epochs],
            )
        )
    return nights


def half_peak(nights, fid):
    v = [r for r in nights if r["fid"] == fid]
    peak = max((r for r in v if r["snr"] >= 5), key=lambda r: r["flux"])
    half = 0.5 * peak["flux"]
    bands = {}
    for z in (0, 2):
        # Compare to a fixed fraction of the measured peak. The peak has its own uncertainty.
        below = [r for r in v if r["flux"] + z * r["error"] < half]
        above = [r for r in v if r["flux"] - z * r["error"] >= half]
        before = max((r["mjd"] for r in below if r["mjd"] < peak["mjd"]), default=None)
        after = min((r["mjd"] for r in below if r["mjd"] > peak["mjd"]), default=None)
        high = [
            r["mjd"]
            for r in above
            if (before is None or r["mjd"] > before)
            and (after is None or r["mjd"] < after)
        ]
        bands[str(z)] = {
            "before_below": before,
            "first_above": min(high, default=None),
            "last_above": max(high, default=None),
            "after_below": after,
            "span_above_days": max(high) - min(high) if high else None,
        }
    return {"observed_peak": peak, "half_ujy": half, "crossings_mjd": bands}


def difference_detections(phot):
    # Brokers can attach a historical alert more than once. Prefer the original stamp;
    # historical rows can have rounded magnitudes. Never count a processing ID twice.
    groups = defaultdict(list)
    for r in phot:
        groups[str(r["pid"])].append(r)
    selected = []
    for pid, items in groups.items():
        reference = max(
            items, key=lambda r: (bool(r.get("has_stamp")), str(r["candid"]))
        )
        if any(
            r["fid"] != reference["fid"]
            or abs(r["mjd"] - reference["mjd"]) > 1 / 86400
            or r["isdiffpos"] != reference["isdiffpos"]
            or abs(r["magpsf"] - reference["magpsf"]) > 0.01
            for r in items
        ):
            raise ValueError(f"Inconsistent repeated alert {pid}")
        selected.append(reference)
    return sorted(selected, key=lambda r: r["mjd"])


def run(fetch=False):
    raw, det, limits = inputs(fetch)
    nights = group_forced(raw)
    exposures = difference_detections(det)
    _, hist, _ = compass_archive.load(YEAR / "snapshot")
    broker_by_id = {r["candid"]: r for r in det}
    for f in hist[NAME]:
        r = broker_by_id[str(f["i:candid"])]
        if (
            r["fid"] != f["i:fid"]
            or abs(r["mjd"] + 2400000.5 - f["i:jd"]) > 1 / 86400
            or abs(r["magpsf"] - f["i:magpsf"]) > 0.001
        ):
            raise ValueError("Fink exposure mismatch")
    by_day = defaultdict(dict)
    for r in nights:
        if r["snr"] >= 5:
            by_day[r["day"]][r["fid"]] = r
    colors = [
        {
            "day": d,
            "date": v[1]["date"],
            "g_minus_r": v[1]["mag"] - v[2]["mag"],
            "error": math.hypot(v[1]["magerr"], v[2]["magerr"]),
            "hours_apart": 24 * abs(v[1]["mjd"] - v[2]["mjd"]),
        }
        for d, v in sorted(by_day.items())
        if 1 in v and 2 in v
    ]
    summary = {}
    for band in (1, 2):
        for key, a, b in [
            ("before", "2025-03-01", "2025-04-18"),
            ("discovery", "2025-04-19", "2025-04-23"),
            ("late_may", "2025-05-20", "2025-05-31"),
        ]:
            rr = [r for r in nights if r["fid"] == band and a <= r["date"] <= b]
            if rr:
                sw = sum(1 / r["error"] ** 2 for r in rr)
                summary[f"{band}_{key}"] = {
                    "dates": len(rr),
                    "mean_flux_ujy": sum(r["flux"] / r["error"] ** 2 for r in rr) / sw,
                    "formal_mean_error": sw**-0.5,
                }
    return {
        "alert_rows": len(det),
        "alert_exposures": len(exposures),
        "fink_exposures": len(hist[NAME]),
        "good_alert_exposures": sum(
            r["isdiffpos"] == 1 and 0 < r["sigmapsf"] <= 0.20 for r in exposures
        ),
        "limit_rows": len(limits),
        "forced_rows": len(raw),
        "used_forced": sum(len(r["pids"]) for r in nights),
        "first_and_last_limit_mjd": (
            [min(r["mjd"] for r in limits), max(r["mjd"] for r in limits)]
            if limits
            else None
        ),
        "nightly_forced": nights,
        "near_simultaneous_colors": colors,
        "half_peak": {str(f): half_peak(nights, f) for f in (1, 2)},
        "window_flux": summary,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "results/onset.json")
    args = p.parse_args()
    pilot.json_write(args.out, run(args.fetch))
