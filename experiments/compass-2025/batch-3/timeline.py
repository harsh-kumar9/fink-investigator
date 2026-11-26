"""Crossmatch the Fink alerts, broker forced light, and DR24 source detections."""

import argparse
import json
from pathlib import Path
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import compass_archive
import compass_direct as direct
import pilot

NAME = "ZTF25aaohzdc"
TAP = "https://tap.alerce.online/tap/sync"
BASE = f"https://api.alerce.online/ztf/v1/objects/{NAME}"
FIELDS = "oid,pid,mjd,fid,mag,e_mag,mag_corr,e_mag_corr_ext,isdiffpos,corrected,dubious,field,rcid,procstatus,distnr,magnr"


def run(fetch=False):
    folder = HERE / "onset-sources"
    query = f"SELECT TOP 10000 {FIELDS} FROM ztf.forced_photometry WHERE oid='{NAME}'"
    url = (
        TAP
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    )
    forced = direct.capture(folder, "alerce-forced", url, set(FIELDS.split(",")), fetch)
    if (
        len(forced) >= 10000
        or any(r["oid"] != NAME for r in forced)
        or len({r["pid"] for r in forced}) != len(forced)
    ):
        raise ValueError("Incorrect or truncated forced series")
    alert = compass_archive.capture_list(
        folder, "alerce-detections", BASE + "/detections", fetch
    )
    limits = compass_archive.capture_list(
        folder, "alerce-limits", BASE + "/non_detections", fetch
    )
    ledger = json.loads((HERE / "results/ledger.json").read_text())
    (case,) = [r for r in ledger["objects"] if r["id"] == NAME]
    _, rows, _ = direct.direct_data(case, HERE / "direct" / NAME)
    direct_exps = {r["expid"]: r for r in rows}
    if len(direct_exps) != len(rows):
        raise ValueError("Duplicate exposure in direct observations")
    _, histories, _ = compass_archive.load(HERE.parent / "snapshot")
    fink = {str(r["i:candid"]): r for r in histories[NAME]}
    data = {}
    for r in alert:
        pid = str(r["pid"])
        if pid in data:
            raise ValueError("Duplicate broker alert")
        out = {
            "pid": pid,
            "mjd": r["mjd"],
            "band": {1: "g", 2: "r"}[r["fid"]],
            "candid": r["candid"],
            "difference_mag": r["magpsf"],
            "alert_corrected_mag": (
                r["magpsf_corr"] if r["corrected"] and not r["dubious"] else None
            ),
        }
        if r["candid"] in fink:
            row = fink[r["candid"]]
            if (
                row["i:fid"] != r["fid"]
                or abs(row["i:jd"] - 2400000.5 - r["mjd"]) > 1 / 86400
                or row["i:magpsf"] != r["magpsf"]
            ):
                raise ValueError("Fink exposure mismatch")
            out["in_fink"] = True
        data[pid] = out
    if set(fink) - {r["candid"] for r in alert}:
        raise ValueError("Missing Fink alert in broker series")
    for r in forced:
        pid = r["pid"]
        t = float(r["mjd"])
        band = {1: "g", 2: "r"}[int(r["fid"])]
        out = data.setdefault(pid, {"pid": pid, "mjd": t, "band": band})
        if abs(t - out["mjd"]) > 1 / 86400 or out["band"] != band:
            raise ValueError("Forced epoch mismatch")
        if (
            r["corrected"] == "True"
            and r["dubious"] == "False"
            and r["procstatus"] == "0"
            and float(r["distnr"]) < 1.5
        ):
            out.update(
                forced_mag=float(r["mag_corr"]), forced_err=float(r["e_mag_corr_ext"])
            )
    # Alert-linked forced catalogs omit exposures after the broker stops alerting.
    # Include DR24-only measurements in the registered event window as well.
    linked = {str(int(pid) // 10000) for pid in data}
    for expid, row in direct_exps.items():
        mjd = float(row["mjd"])
        date = str(pilot.day_of(mjd + 2400000.5))
        if expid not in linked and "2025-04-01" <= date <= "2025-06-15":
            data["direct-" + expid] = {
                "expid": expid,
                "mjd": mjd,
                "band": row["filtercode"][1],
            }
    for key, out in data.items():
        if "expid" not in out:
            out["expid"] = str(int(key) // 10000)
        row = direct_exps.get(out["expid"])
        if row:
            if (
                row["filtercode"] != "z" + out["band"]
                or abs(float(row["mjd"]) - out["mjd"]) > 1 / 86400
            ):
                raise ValueError("Direct exposure mismatch")
            out.update(
                direct_mag=float(row["mag"]),
                direct_err=float(row["magerr"]),
                direct_catflags=int(row["catflags"]),
            )
        out["date"] = str(pilot.day_of(out["mjd"] + 2400000.5))
    return {
        "fink_alerts": len(fink),
        "broker_alerts": len(alert),
        "broker_forced": len(forced),
        "returned_limits": limits,
        "measurements": sorted(data.values(), key=lambda r: r["mjd"]),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch))
