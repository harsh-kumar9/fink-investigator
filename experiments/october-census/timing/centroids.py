"""Locate the missing light in negative ZTF difference-image detections."""

from collections import defaultdict
import json
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from epochs import CENSUS, HERE, NAME, alert_source, heliocentric


def run():
    _, rows = alert_source()
    result = json.loads((HERE / "dips.json").read_text())
    fit = result["fits"][result["training_winner"]]
    usable = [
        r
        for r in rows
        if r["i:fid"] in {1, 2}
        and (r["i:drb"] or 0) >= 0.8
        and 0 < (r["i:sigmapsf"] or 99) <= 0.2
        and str(r["i:isdiffpos"]).lower() in {"f", "false", "0"}
    ]
    hjd = heliocentric(np.array([r["i:jd"] - 2400000.5 for r in usable]))
    phase = (hjd - fit["center_hjd"]) / fit["period_days"]
    inside = abs(phase - np.round(phase)) <= fit["duration_days"] / (
        2 * fit["period_days"]
    )
    groups = defaultdict(list)
    for r, keep in zip(usable, inside):
        if keep:
            groups[int(r["i:jd"] - 2400000.5)].append(
                [r["i:jd"], r["i:ra"], r["i:dec"]]
            )
    nightly = np.array([np.median(v, axis=0) for _, v in sorted(groups.items())])
    years = 2000 + (nightly[:, 0] - 2451545.0) / 365.25 - 2016
    positions = SkyCoord(nightly[:, 1], nightly[:, 2], unit=u.deg)
    catalog = json.loads((CENSUS / "results/extended.json").read_text())
    stars = next(r["tables"]["I/355/gaiadr3"] for r in catalog if r["id"] == NAME)
    measured = []
    for s in stars:
        ra, dec = float(s["RA_ICRS"]), float(s["DE_ICRS"])
        predicted = SkyCoord(
            ra + years * float(s["pmRA"]) / 3.6e6 / np.cos(np.deg2rad(dec)),
            dec + years * float(s["pmDE"]) / 3.6e6,
            unit=u.deg,
        )
        offset = positions.separation(predicted).arcsec
        measured.append(
            {
                "source": s["Source"],
                "median_arcsec": float(np.median(offset)),
                "p90_arcsec": float(np.percentile(offset, 90)),
            }
        )
    return {
        "negative_in_window_detections": int(inside.sum()),
        "dates": len(nightly),
        "first_jd": nightly[0, 0],
        "last_jd": nightly[-1, 0],
        "gaia": measured,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
