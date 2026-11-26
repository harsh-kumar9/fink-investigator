"""Plot historical target photometry and the same-image comparison stars."""

from datetime import datetime, timedelta
import json
import math
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ORIGIN = datetime(1858, 11, 17)


def plot():
    data = json.loads((HERE / "early-state.json").read_text())
    fig, (series, field) = plt.subplots(2, 1, figsize=(5.5, 4.9), height_ratios=[1.7, 1], layout="constrained")
    fig.suptitle("MGAB-V306 on 3 November 2013", fontsize=15)
    for index, (band, color) in enumerate([("z", "#9e3e28"), ("y", "#235b87")]):
        for obs in data["source_exposures"][band]:
            day = ORIGIN + timedelta(days=obs["mjd"])
            mag = obs["mag_if_SNR5"]
            if mag is None:
                series.plot(day, obs["depth_5sigma"], marker="v", mfc="none", mec=color)
                continue
            series.errorbar(day, mag, (2.5 / math.log(10)) * obs["error_uJy"] / obs["flux_uJy"],
                            fmt="o" if band == "z" else "s", ms=6 if band == "z" else 4,
                            color=color, capsize=2)
        series.plot([], [], "o" if band == "z" else "s", color=color, label=f"{band} band")
        controls = []
        for source, stats in data["comparisons"][band].items():
            change = stats["event_minus_other_nights_median_mag"]
            if source == data["target"]:
                field.scatter(change, index, marker="D", s=55, color=color,
                              label="MGAB-V306" if index == 0 else None)
            else:
                controls.append(change)
        field.scatter(controls, [index + (i - (len(controls) - 1) / 2) * 0.026 for i in range(len(controls))],
                      s=24, color="#64748b", label="Field stars" if index == 0 else None)
    series.set(ylim=(20.35, 16.2), ylabel="PS1 forced magnitude")
    series.xaxis.set_major_locator(mdates.YearLocator())
    series.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    series.text(0.98, 0.80, "3 Nov", transform=series.transAxes, ha="right", fontsize=9)
    series.legend(frameon=False, loc="upper left", ncol=2, fontsize=9)
    field.set(yticks=[0, 1], yticklabels=["z", "y"], xlim=(-3.2, 0.28), ylim=(-0.4, 1.6),
              xlabel="3 Nov minus earlier-night median (mag)")
    field.axvline(0, color="#94a3b8", lw=1, zorder=0)
    field.legend(frameon=False, ncol=2, loc="upper center", fontsize=8.5)
    for ax in [series, field]:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#e2e8f0", lw=0.7)
        ax.set_axisbelow(True)
    fig.savefig(HERE / "november-2013.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    plot()
