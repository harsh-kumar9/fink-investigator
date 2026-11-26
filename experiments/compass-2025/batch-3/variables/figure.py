"""Plot direct ZTF eruption dates and ASAS-SN recurrence measurements."""

from datetime import date
import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from analyze import HERE, curves

_, data = curves("ZTF17aabshbq")
records = json.loads((HERE / "recurrence.json").read_text())["years"]
erupt = json.loads((HERE / "eruptions.json").read_text())
fig, (a, b) = plt.subplots(2, 1, figsize=(6.7, 6.5), layout="constrained")
dates = [d for d in data["r"] if date(2018, 10, 20) <= d <= date(2018, 12, 12)]
mag = [data["r"][d][1] for d in dates]
for i in range(len(dates) - 1):
    if (dates[i + 1] - dates[i]).days <= 5:
        a.plot(dates[i : i + 2], mag[i : i + 2], c="#50738c", lw=1)
a.scatter(dates, mag, s=18, color="#294b62")
a.axhline(
    erupt["training_thresholds_mag"]["r"]["bright_q25"],
    color="#006d5b",
    ls=":",
    label="Bright-state threshold",
)
for day in ["2018-11-02", "2018-11-19", "2018-12-04"]:
    d = date.fromisoformat(day)
    a.axvline(d, color="#ad5034", lw=1, alpha=0.75)
a.invert_yaxis()
a.set_ylabel("ZTF red magnitude")
a.set_title("Three bright episodes in autumn 2018", loc="left", fontsize=13)
a.xaxis.set_major_locator(mdates.DayLocator(interval=10))
a.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
a.legend(loc="lower left", frameon=False, fontsize=9)
years = [y for y, v in records.items() if v["common_anchors"] >= 10]
for y in years:
    vals = records[y]["delta_flux_medians_mjy"]
    b.plot([0, 1, 2], vals, marker="o", label=y, linewidth=1.5)
b.set_xticks([0, 1, 2], ["6–10 days", "14–20 days", "24–29 days"])
b.set_xlim(-0.12, 2.45)
b.set_ylim(0, 0.68)
b.set_ylabel("Median absolute flux\ndifference (mJy)")
b.set_title("ASAS-SN nights 14–20 days apart were more alike", loc="left", fontsize=13)
b.legend(ncols=1, frameon=False, fontsize=9, loc="upper right")
for ax in (a, b):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.15)
fig.savefig(HERE / "recurrence.png", dpi=155)
