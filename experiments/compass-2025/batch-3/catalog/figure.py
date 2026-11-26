"""Plot Hambsch's two quality-selected time series and the frozen injection results."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from timing import aavso, catalog, selected

HERE = Path(__file__).resolve().parent


def main():
    rows = selected(aavso())
    report = json.loads((HERE / "timing.json").read_text())
    fig, axes = plt.subplots(3, 1, figsize=(5.4, 6.3), layout="constrained",
                             gridspec_kw={"height_ratios": [1.2, 1.2, 0.95]})
    for ax, (day, night) in zip(axes[:2], report["nights"].items()):
        series = [r for r in rows if str(catalog.pilot.day_of(float(r["JD"]))) == day]
        x = [(float(r["JD"])-night["first_jd"])*1440 for r in series]
        y = [float(r["mag"]) for r in series]
        errors = [float(r["uncert"]) for r in series]
        ax.errorbar(x, y, yerr=errors, ls="none", marker="o", markersize=3,
                    color="#134866", ecolor="#a4b5c1", elinewidth=1, capsize=0)
        ax.axhline(night["median_mag"], color="#bd531b", lw=1, ls="--")
        ax.invert_yaxis()
        ax.set_xlim(-2, 120)
        ax.set_ylabel("CV magnitude")
        ax.set_title(f"April {int(day[-2:])}  ·  {len(y)} exposures", loc="left", fontsize=11, weight="bold")
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.grid(axis="x", color="#ddd", linewidth=0.5)
    axes[1].set_xlabel("Minutes from each night's first exposure")
    ax = axes[-1]
    draws = report["injections"]
    heights = [v["transfer_pass"] for v in draws.values()]
    ax.bar(range(3), heights, width=0.45, color="#527e72")
    for i, n in enumerate(heights):
        ax.text(i, n + 0.6, f"{n}/60", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(range(3), ["0", "0.03", "0.10"])
    ax.set_ylim(0, 15)
    ax.set_ylabel("Tests passed")
    ax.set_xlabel("Injected semi-amplitude (mag), 90-minute period")
    ax.set_title("Artificial signals passing the overnight test", loc="left", fontsize=11, weight="bold")
    for ax in axes:
        ax.spines[["right", "top"]].set_visible(False)
    fig.savefig(HERE / "timing.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
