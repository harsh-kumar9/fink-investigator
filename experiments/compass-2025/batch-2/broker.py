"""Compare forced, alert and direct measurements on the exact January exposures."""

import argparse
from datetime import date
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import compass_archive
import compass_direct as direct
import pilot

NAME = "ZTF25aacfjde"
BASE = "https://api.alerce.online/ztf/v1/objects/" + NAME
TAP = "https://tap.alerce.online/tap/sync"
FIELDS = "oid,pid,mjd,fid,mag,e_mag,mag_corr,e_mag_corr_ext,isdiffpos,corrected,dubious,field,rcid,sciinpseeing,diffmaglim,procstatus,distnr,magnr"


def run(fetch=False):
    folder = HERE / "onset-field"
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
        raise ValueError("Wrong or truncated forced series")
    _, histories, _ = compass_archive.load(HERE.parent / "snapshot")
    fink = {str(r["i:candid"]): r for r in histories[NAME]}
    alert = raw_json(folder, "alerce-detections", BASE + "/detections", fetch)
    limits = raw_json(folder, "alerce-limits", BASE + "/non_detections", fetch)
    obs = json.loads((HERE / "results" / "ledger.json").read_text())
    (case,) = [r for r in obs["objects"] if r["id"] == NAME]
    _, raw, _ = direct.direct_data(case, HERE / "direct" / NAME)
    by_expid = {r["expid"]: r for r in raw}
    out = {}
    for d in (date(2025, 1, 17), date(2025, 1, 19)):
        rows = []
        for f in sorted(forced, key=lambda r: float(r["mjd"])):
            t = float(f["mjd"])
            if pilot.day_of(t + 2400000.5) != d:
                continue
            if (
                f["procstatus"] != "0"
                or f["dubious"] != "False"
                or f["corrected"] != "True"
                or float(f["distnr"]) >= 1.5
            ):
                raise ValueError("Unreliable forced measurement")
            # A ZTF PID appends the readout channel and processing code to the exposure ID.
            expid = str(int(f["pid"]) // 10000)
            r = by_expid.get(expid)
            if (
                r is None
                or int(r["catflags"])
                or abs(float(r["mjd"]) - t) > 1 / 86400
                or r["filtercode"] != {1: "zg", 2: "zr"}[int(f["fid"])]
            ):
                raise ValueError("Missing exact direct exposure")
            row = {
                "band": r["filtercode"][1:],
                "expid": expid,
                "mjd": t,
                "direct_mag": float(r["mag"]),
                "direct_err": float(r["magerr"]),
                "forced_mag": float(f["mag_corr"]),
                "forced_err": float(f["e_mag_corr_ext"]),
                "seeing_arcsec": float(f["sciinpseeing"]),
                "reference_mag": float(f["magnr"]),
            }
            found = [a for a in alert if str(a["pid"]) == f["pid"]]
            limit = [
                a
                for a in limits
                if abs(float(a["mjd"]) - t) < 1 / 86400 and a["fid"] == int(f["fid"])
            ]
            if found:
                (a,) = found
                original = fink[a["candid"]]
                if (
                    original["i:magpsf"] != a["magpsf"]
                    or not a["corrected"]
                    or a["dubious"]
                    or a["isdiffpos"] != 1
                ):
                    raise ValueError("Inconsistent alert correction")
                delta = a["magpsf"]
                recon = -2.5 * math.log10(
                    10 ** (-0.4 * row["reference_mag"]) + 10 ** (-0.4 * delta)
                )
                row.update(
                    difference_mag=delta,
                    alert_corrected=a["magpsf_corr"],
                    alert_corrected_err=a["sigmapsf_corr_ext"],
                    reconstructed_total_mag=recon,
                )
            if limit:
                (item,) = limit
                row["limiting_difference_mag"] = item["diffmaglim"]
            rows.append(row)
        if len(rows) != 2 or len({r["band"] for r in rows}) != 2:
            raise ValueError("Expected two forced exposures on each date")
        out[str(d)] = {r["band"]: r for r in rows}
    change = {}
    for b in "gr":
        a = out["2025-01-17"][b]
        c = out["2025-01-19"][b]
        change[b] = {
            "forced_brightening_mag": a["forced_mag"] - c["forced_mag"],
            "errors_in_quadrature_mag": math.hypot(a["forced_err"], c["forced_err"]),
        }
    return {
        "forced_epochs_returned": len(forced),
        "detection_count": len(alert),
        "difference_limits": len(limits),
        "dates": out,
        "january_change": change,
    }


def raw_json(folder, name, url, fetch):
    path = folder / (name + ".json.gz")
    receipt = folder / (name + "-receipt.json")
    if fetch and not receipt.exists():
        raw = pilot.curl(url, [])
        if not isinstance(json.loads(raw), list):
            raise ValueError("Unexpected broker response")
        from datetime import datetime, timezone

        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(
            receipt,
            {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            },
        )
    raw = gzip.decompress(path.read_bytes())
    m = json.loads(receipt.read_text())
    if m["url"] != url or hashlib.sha256(raw).hexdigest() != m["sha256_uncompressed"]:
        raise ValueError("Broker receipt mismatch")
    return json.loads(raw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2, allow_nan=False))
