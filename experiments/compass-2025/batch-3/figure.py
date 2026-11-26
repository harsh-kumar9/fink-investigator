"""Plot the March discovery, April observer photometry and the original ZTF images."""

from datetime import datetime, timedelta
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_text())


def jd_date(jd):
    return datetime(2000, 1, 1, 12) + timedelta(days=jd - 2451545.0)


def main():
    source = read("catalog/sources.json")["discovery"]
    apr = read("catalog/timing.json")["nights"]
    ztf = read("results/apertures.json")["photometry"]
    colors = {"g": "#00735f", "r": "#af3b15", "CV": "#573c8b", "ASAS": "#254e8e"}
    fig, axes = plt.subplots(2, 1, figsize=(5.45, 5.8))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.78, bottom=0.10, hspace=0.40)
    fig.suptitle(
        "ASASSN-25ba brightened again in May",
        fontsize=14,
        x=0.07,
        ha="left",
        y=0.98,
    )
    # The discovery magnitude has no explicit filter in the survey row.
    for sign, mag, day, fraction in re.findall(
        r"g([=>])([\d.]+) on 2025-03-(\d\d)\.(\d+)", source["survey_note"]
    ):
        t = datetime(2025, 3, int(day)) + timedelta(days=float("0." + fraction))
        axes[0].plot(
            t,
            float(mag),
            marker="v" if sign == ">" else "s",
            color=colors["ASAS"],
            markerfacecolor="none" if sign == ">" else colors["ASAS"],
            ms=6,
        )
    for d, v in apr.items():
        axes[0].plot(
            jd_date(0.5 * (v["first_jd"] + v["last_jd"])),
            v["median_mag"],
            "D",
            color=colors["CV"],
            ms=5.5,
        )
    for ax, (lo, hi, title) in zip(
        axes,
        [
            ("03-07", "04-28", "March discovery and April decline"),
            ("04-22", "06-12", "ZTF caught a second brightening on May 16"),
        ],
    ):
        left, right = [datetime.fromisoformat("2025-" + v) for v in (lo, hi)]
        for r in ztf:
            t = jd_date(r["obsjd"])
            if (
                r["window"] != "event"
                or r["status"] != "measured"
                or not left <= t <= right
            ):
                continue
            color = colors[r["band"]]
            if r["mag"] is None:
                ax.plot(
                    t, r["five_sigma_limit"], marker="v", color=color, mfc="none", ms=5
                )
            else:
                ax.errorbar(
                    t,
                    r["mag"],
                    yerr=r["magerr"],
                    fmt="o" if r["band"] == "g" else "^",
                    color=color,
                    ms=5,
                    capsize=2,
                )
        ax.set(xlim=(left, right), ylim=(22.2, 14.1), ylabel="Magnitude")
        ax.set_title(title, loc="left", fontsize=11, fontweight="semibold", pad=6)
        ax.set_yticks([15, 17, 19, 21])
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=12))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.tick_params(labelsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#d9dde1", linewidth=0.7)
    handles = [
        Line2D([], [], marker=m, ls="", color=colors[k], label=t, ms=5.5)
        for m, k, t in [
            ("s", "ASAS", "ASAS-SN g"),
            ("D", "CV", "Hambsch CV"),
            ("o", "g", "ZTF g"),
            ("^", "r", "ZTF r"),
        ]
    ]
    handles.append(
        Line2D(
            [],
            [],
            marker="v",
            ls="",
            markerfacecolor="none",
            color="#424851",
            label="Fainter than limit",
        )
    )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.055, 0.924),
        frameon=False,
        ncol=3,
        columnspacing=1.35,
        handletextpad=0.25,
        fontsize=9,
    )
    axes[1].set_xlabel("Observation date in 2025 (UTC)", fontsize=10)
    fig.savefig(HERE / "results/light-curve.png", dpi=170, facecolor="white")


if __name__ == "__main__":
    main()
