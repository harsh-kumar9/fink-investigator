"""Draw held-out ZTF and earlier Gaia evidence for V1153 Aql."""

from datetime import date
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata

from check import BATCH, HERE, NAME, compass_direct, design, fit, nightly

case = next(
    r
    for r in json.loads((BATCH / "results/ledger.json").read_text())["objects"]
    if r["id"] == NAME
)
_, rows, _ = compass_direct.direct_data(case, BATCH / "direct" / NAME)
curves = nightly(rows, case["jd"] + 1e-7)
timing = json.loads((HERE / "results.json").read_text())["timing"]
gaia = json.loads((HERE / "gaia-results.json").read_text())
cutoff = date.fromisoformat(timing["training_last"])
f = 1 / timing["models"]["training_best"]["period_days"]
t, y = np.asarray([v for d, v in curves["r"].items() if d <= cutoff]).T
origin = t[0]
_, coeff = fit(t - origin, y, f)
fig, (ax, bx) = plt.subplots(
    2,
    1,
    figsize=(6.8, 6.7),
    layout="constrained",
    gridspec_kw={"height_ratios": [1, 1]},
)
days, values = zip(*((d, v) for d, v in curves["r"].items() if d > cutoff))
test_t, test_y = np.array(values).T
ax.scatter(days, test_y, s=13, color="#202934", label="Later ZTF red nights")
tt = np.linspace(test_t[0], test_t[-1], 1000)
dates = days[0].toordinal() + tt - test_t[0] - date(1970, 1, 1).toordinal()
ax.plot(
    dates,
    design(tt - origin, f) @ coeff,
    color="#b03d3a",
    lw=2,
    label="Fixed early ZTF fit",
)
ax.invert_yaxis()
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.38, 1.0))
ax.set_ylabel("Red magnitude")
ax.set_title(
    "A waveform fitted through August 2023 predicts later ZTF nights",
    loc="left",
    fontsize=12,
)
ax.tick_params(axis="x", labelrotation=20)
bx.plot([0, 1], [0, 1], linestyle=":", color="#77818b", lw=1)
for band, label, color, marker in (
    ("bp", "Gaia blue (BP)", "#00728b", "o"),
    ("rp", "Gaia red (RP)", "#bf5a20", "x"),
):
    tx, mag, n = np.asarray(gaia["bins"][band]).T
    pred = design(tx - origin, f) @ coeff
    norm = lambda arr: (rankdata(arr) - 1) / (len(arr) - 1)
    bx.scatter(
        norm(pred), norm(mag), marker=marker, s=37, color=color, alpha=0.85, label=label
    )
bx.set_title(
    "The same waveform ranks Gaia visits years earlier", loc="left", fontsize=12
)
bx.set_xlabel("Predicted magnitude rank (bright to faint)")
bx.set_ylabel("Gaia magnitude rank")
bx.legend(frameon=False)
for a in (ax, bx):
    a.spines[["top", "right"]].set_visible(False)
    a.grid(axis="y", alpha=0.16)
fig.savefig(HERE / "timing.png", dpi=150)
