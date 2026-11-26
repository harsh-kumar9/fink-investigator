"""Apply the earlier dwarf-nova waveform test to an independently selected star."""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import tarfile
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "compass-2025/batch-4/variables"))
from template import PHASES, assess, compass_direct, pilot

URL = "https://raw.githubusercontent.com/astroML/astroML-data/main/datasets/RRLyr_ugriz_templates.tar.gz"
NAME = "ZTF18abbtzhm"


def shapes(model, fetch):
    root = HERE / "catalog"
    files = {band: model["T" + band] + band + ".dat" for band in "gr"}
    path, receipt = root / "shapes.json.gz", root / "shapes.json"
    if fetch and not receipt.exists():
        raw = pilot.curl(URL, [])
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
            content = {
                b: archive.extractfile(f).read().decode() for b, f in files.items()
            }
        data = json.dumps(content).encode()
        path.write_bytes(gzip.compress(data, mtime=0))
        pilot.json_write(
            receipt,
            dict(
                url=URL,
                members=files,
                sha256=hashlib.sha256(data).hexdigest(),
                archive_sha256=hashlib.sha256(raw).hexdigest(),
                retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ),
        )
    data = gzip.decompress(path.read_bytes())
    meta = json.loads(receipt.read_text())
    if (
        meta["url"] != URL
        or meta["members"] != files
        or meta["sha256"] != hashlib.sha256(data).hexdigest()
    ):
        raise ValueError("Template receipt mismatch")
    return {b: np.loadtxt(io.StringIO(v)).T for b, v in json.loads(data).items()}


def run(fetch):
    (model,) = json.loads((HERE / "catalog/priors.json").read_text())["sesar_5arcsec"]
    period = float(model["Per"])
    curve = shapes(model, fetch)
    obj = next(
        r
        for r in json.loads((HERE / "results/ledger.json").read_text())["objects"]
        if r["id"] == NAME
    )
    _, rows, _ = compass_direct.direct_data(obj, HERE / "direct" / NAME)
    nights = defaultdict(dict)
    for row in rows:
        t, mag, err = map(float, [row["hjd"], row["mag"], row["magerr"]])
        if (
            int(row["catflags"])
            or not 0 < err <= 0.2
            or not all(map(math.isfinite, (t, mag)))
        ):
            continue
        jd = float(row["mjd"]) + 2400000.5
        b = row["filtercode"][-1]
        phase = ((t - 2400000) / period + float(model["phi0"]) + PHASES) % 1
        predicted = float(model[b + "Amp"]) * np.interp(phase, *curve[b], period=1)
        key = (
            "historical" if jd <= obj["jd"] + 1e-7 else "after_anchor",
            b,
            pilot.day_of(jd),
        )
        nights[key][row["expid"]] = np.r_[mag, predicted]
    results = {}
    for window in ("historical", "after_anchor"):
        data = [
            (b, d, np.mean(list(exps.values()), axis=0))
            for (win, b, d), exps in sorted(nights.items())
            if win == window
        ]
        results[window] = dict(
            first=str(min(d for b, d, r in data)),
            last=str(max(d for b, d, r in data)),
            fits=[assess(data, w) for w in (14, 7, 28)],
        )
    idx = int(np.argmax(results["historical"]["fits"][0]["amplitude_at_phase"]))
    transfer = results["after_anchor"]["fits"][0]
    return dict(
        id=NAME,
        model=model,
        windows=results,
        early_best_shift=PHASES[idx],
        late_scale_at_early_shift=transfer["amplitude_at_phase"][idx],
        late_interval_at_early_shift=[
            r[idx] for r in transfer["pointwise_bootstrap_p025_p975"]
        ],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--out", type=Path, default=HERE / "results/waveform.json")
    a = p.parse_args()
    pilot.json_write(a.out, run(a.fetch))
