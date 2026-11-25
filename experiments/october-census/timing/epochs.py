"""Compare the locked ZTF waveform with earlier PS1 and later Fink epochs."""

from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from statistics import median

from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy.utils import iers
import numpy as np
from scipy.stats import spearmanr

from analyze import CENSUS, HERE, nights

NAME = "ZTF18abjtfpx"
POS = SkyCoord(280.305911, 9.7241705, unit=u.deg)
iers.conf.auto_download = False


class Ephemeris:
    def __init__(self):
        cohort = json.loads((CENSUS / "results/direct.json").read_text())["exploration"]
        self.data = nights(next(row for row in cohort if row["id"] == NAME))
        rows = json.loads((HERE / "results.json").read_text())["period_search"]
        result = next(row for row in rows if row["id"] == NAME)
        self.frequency = 1 / result["period_days"]
        self.n = result["training_nights"]
        self.epoch = self.data[0, 1]
        self.train_t = self.data[: self.n, 1] - self.epoch
        self.model = LombScargle(self.train_t, self.data[: self.n, 2], nterms=2)

    def mag(self, hjd):
        return self.model.model(np.asarray(hjd) - self.epoch, self.frequency)

    def phase(self, hjd):
        return ((np.asarray(hjd) - self.epoch) * self.frequency) % 1


def heliocentric(mjd):
    """Geocentric UTC plus the heliocentric light-time (topocenter changes <22 ms)."""
    time = Time(mjd, format="mjd", scale="utc")
    vector = get_body_barycentric("earth", time) - get_body_barycentric("sun", time)
    return time.utc.jd + (vector.dot(POS.cartesian) / c).to_value(u.day)


def table(name):
    receipt = json.loads((HERE / "ps1-sources.json").read_text())[name]
    path = HERE / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != receipt["sha256"]:
        raise ValueError(f"Changed PS1 response {name}")
    with path.open(newline="") as file:
        return list(csv.DictReader(file))


def ps1(eph):
    sources = []
    for row in table("ps1-mean.csv"):
        coord = SkyCoord(float(row["raMean"]), float(row["decMean"]), unit=u.deg)
        sources.append(
            {
                "id": row["objID"],
                "offset_arcsec": POS.separation(coord).arcsec,
                "detections": int(row["nDetections"]),
            }
        )
    sources.sort(key=lambda r: r["offset_arcsec"])
    target = next(
        r for r in sources if r["offset_arcsec"] <= 1 and r["detections"] >= 2
    )
    grouped = defaultdict(list)
    raw = table("ps1-detections.csv")
    if any(r["objID"] != target["id"] for r in raw) or len(
        {r["detectID"] for r in raw}
    ) != len(raw):
        raise ValueError("PS1 object or detection ID mismatch")
    for row in raw:
        flag, flag3 = int(row["infoFlag"]), int(row["infoFlag3"])
        if (
            not (flag & 1)
            or flag & (8 | 16 | 128 | 1024 | 2048 | 4096)
            or flag3 & (1 | 4 | 8192)
        ):
            continue
        flux, err = float(row["psfFlux"]), float(row["psfFluxErr"])
        if flux <= 0 or not 0 < err / flux <= 0.1 or float(row["psfQfPerfect"]) < 0.85:
            continue
        t = float(row["obsTime"])
        grouped[(int(row["filterID"]), math.floor(t))].append(
            (t, -2.5 * math.log10(flux / 3631))
        )
    bands = {}
    for band, name in enumerate("grizy", start=1):
        pairs = [
            [median(t for t, _ in vals), median(y for _, y in vals)]
            for (f, day), vals in sorted(grouped.items())
            if f == band
        ]
        bands[name] = {"nights": len(pairs), "mjd_mag": pairs}
        if len(pairs) >= 5:
            t, y = np.array(pairs).T
            bands[name]["spearman"] = float(
                spearmanr(eph.mag(heliocentric(t)), y).statistic
            )
    return {
        "sources": sources,
        "retrieved": len(raw),
        "retained": sum(map(len, grouped.values())),
        "bands": bands,
    }


def alert_source():
    folder = CENSUS / "sources/fink"
    manifest = json.loads((folder / "manifest.json").read_text())["sources"]
    filename, meta = next(
        (f, m)
        for f, m in manifest.items()
        if NAME in m["params"]["objectId"].split(",")
    )
    raw = gzip.decompress((folder / filename).read_bytes())
    if hashlib.sha256(raw).hexdigest() != meta["sha256_uncompressed"]:
        raise ValueError("Changed Fink archive")
    with (CENSUS / "results/outcomes.csv").open() as file:
        source = next(r for r in csv.DictReader(file) if NAME in r.values())
    anchor = float(source["jd"])
    rows = [
        r
        for r in json.loads(raw)
        if r["i:objectId"] == NAME and r["i:jd"] and r["i:jd"] <= anchor
    ]
    return filename, rows


def later_alerts(eph):
    filename, rows = alert_source()
    grouped = defaultdict(list)
    for r in rows:
        if (
            r["i:fid"] != 2
            or (r["i:drb"] or 0) < 0.8
            or not 0 < (r["i:sigmapsf"] or 99) <= 0.2
        ):
            continue
        t = r["i:jd"] - 2400000.5
        sign = str(r["i:isdiffpos"]).lower()
        if math.floor(t) <= math.floor(eph.data[-1, 0]) or sign not in {
            "f",
            "false",
            "0",
            "t",
            "true",
            "1",
        }:
            continue
        flux = 10 ** (-0.4 * (r["i:magpsf"] - 20))
        grouped[math.floor(t)].append(
            (t, -flux if sign in {"f", "false", "0"} else flux)
        )
    pairs = [
        [median(t for t, _ in vals), median(f for _, f in vals)]
        for _, vals in sorted(grouped.items())
    ]
    out = {"source": filename, "nights": len(pairs), "mjd_signed_flux": pairs}
    if len(pairs) >= 5:
        t, flux = np.array(pairs).T
        hjd = heliocentric(t)
        out.update(
            first_mjd=t[0],
            last_mjd=t[-1],
            occupied_phase_deciles=len(set((10 * eph.phase(hjd)).astype(int))),
            spearman=float(spearmanr(flux, -eph.mag(hjd)).statistic),
        )
    return out


if __name__ == "__main__":
    fit = Ephemeris()
    print(
        json.dumps(
            {"ps1": ps1(fit), "later_fink": later_alerts(fit)},
            indent=2,
            allow_nan=False,
        )
    )
