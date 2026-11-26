"""Transfer the archived Sesar et al. light-curve shapes to ZTF exposure times."""

import argparse
from collections import defaultdict
from datetime import date
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from xml.etree import ElementTree

import numpy as np

from measure import BATCH, HERE, baseline, compass_direct, load, pilot
from pulsation import NAME, PERIOD

PHASES = np.arange(40) / 40


def catalog(name):
    root = BATCH / "catalog"
    raw = gzip.decompress((root / (name + ".gz")).read_bytes())
    receipt = json.loads((root / (name + ".json")).read_text())
    if hashlib.sha256(raw).hexdigest() != receipt["sha256_uncompressed"]:
        raise ValueError("Catalog data failed checksum")
    return raw


def nightly(case):
    tables = ElementTree.fromstring(catalog("MGAB-V488-literature.xml"))
    table = next(
        t
        for t in tables.findall(".//{*}TABLE")
        if t.get("name") == "J/AJ/153/204/table5"
    )
    cols = [f.get("name") for f in table.findall("{*}FIELD")]
    (row,) = table.findall(".//{*}TR")
    model = dict(
        zip(cols, [td.text.strip() if td.text else None for td in row], strict=True)
    )
    if float(model["Per"]) != PERIOD:
        raise ValueError("Published period changed")
    shapes = {
        b: np.loadtxt(io.StringIO(v)).T
        for b, v in json.loads(catalog("ps1-templates")).items()
    }
    _, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
    dates = defaultdict(dict)
    for r in rows:
        t, m, err = map(float, [r["hjd"], r["mag"], r["magerr"]])
        if (
            int(r["catflags"])
            or not 0 < err <= 0.2
            or not all(map(math.isfinite, (t, m)))
        ):
            continue
        d = pilot.day_of(float(r["mjd"]) + 2400000.5)
        if d > date.fromisoformat(case["date"]):
            continue
        b = r["filtercode"][-1]
        phase = ((t - 2400000) / PERIOD + float(model["phi0"]) + PHASES) % 1
        xp, yp = shapes[b]
        x = float(model[b + "Amp"]) * np.interp(phase, xp, yp, period=1)
        dates[b, d][r["expid"]] = np.r_[m, x]
    return [
        (b, d, np.mean(list(exps.values()), axis=0))
        for (b, d), exps in sorted(dates.items())
    ]


def assess(data, width):
    bins = defaultdict(list)
    for b, d, r in data:
        bins[b, ((d - baseline.EPOCH).days - 58000) // width].append(r)
    local = []
    for (b, k), vals in sorted(bins.items()):
        if len(vals) < 2:
            continue
        centered = np.array(vals) - np.mean(vals, axis=0)
        y, x = centered[:, 0], centered[:, 1:]
        block = (k * width + width // 2) // 90
        local.append((block, len(y), y @ x, np.sum(x * x, axis=0)))
    ids = sorted({k for k, n, c, s in local})
    numerator, denominator = [
        np.array(
            [np.sum([row[col] for row in local if row[0] == k], axis=0) for k in ids]
        )
        for col in (2, 3)
    ]
    amplitude = numerator.sum(axis=0) / denominator.sum(axis=0)
    rng = np.random.default_rng(20250301)
    draw = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), 2000)
    sampled = (draw @ numerator) / (draw @ denominator)
    # Adding the template at its own phase shifts this linear least-squares estimate by exactly one.
    injected = (numerator + denominator).sum(axis=0) / denominator.sum(axis=0)
    return {
        "local_bin_width_days": width,
        "eligible_filter_bins": len(local),
        "filter_nights": sum(r[1] for r in local),
        "resampling_blocks": len(ids),
        "published_phase_scale": round(float(amplitude[0]), 6),
        "phase_shifts": list(PHASES),
        "amplitude_at_phase": list(np.round(amplitude, 6)),
        "pointwise_bootstrap_p025_p975": np.round(
            np.quantile(sampled, [0.025, 0.975], axis=0), 6
        ).tolist(),
        "maximum_over_phases": round(float(max(amplitude)), 6),
        "maximizing_phase_shift": float(PHASES[np.argmax(amplitude)]),
        "bootstrap_p95_maximum": round(
            float(np.quantile(np.max(sampled, axis=1), 0.95)), 6
        ),
        "bootstrap_p025_p975_published_phase": list(
            np.round(np.quantile(sampled[:, 0], [0.025, 0.975]), 6)
        ),
        "injection_scale_range": list(np.round([min(injected), max(injected)], 6)),
        "injection_0p7_to_1p3_count_of_40": int(
            np.sum((injected >= 0.7) & (injected <= 1.3))
        ),
    }


def run():
    case, pre, _ = load(NAME)
    cutoff = list(pre["r"])[len(pre["r"]) * 2 // 3 - 1]
    data = nightly(case)
    return {
        "position": NAME,
        "period_days": PERIOD,
        "training_cutoff": str(cutoff),
        "full_pre_anchor": [assess(data, w) for w in (14, 7, 28)],
        "later_only": [
            assess([r for r in data if r[1] > cutoff], w) for w in (14, 7, 28)
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE / "template.json")
    pilot.json_write(parser.parse_args().out, run())
