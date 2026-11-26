"""Plot the six-day rise and same-image comparison measurements."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    field = json.loads((HERE / "field.json").read_text())
    fig, axes = plt.subplots(2, 1, figsize=(5.8, 5.5), height_ratios=[1.55, 1], layout="constrained")
    fig.suptitle("MGAB-V488 in December 2013", fontsize=16)
    colors = {"r": "#b34324", "i": "#234e91"}
    ax, bottom = axes
    for b, shift in [("r", -0.12), ("i", 0.12)]:
        subset = [r for r in field["target_exposures"] if r["band"] == b]
        for r in subset:
            same_date = [row for row in subset if int(row["mjd"]) == int(r["mjd"])]
            spread = (same_date.index(r) - (len(same_date) - 1) / 2) * .13
            x = int(r["mjd"]) - 56632 + shift + spread
            marker = "s" if r["defect_flag"] else "o"
            face = colors[b] if r["passes_strict"] and not r["defect_flag"] else "white"
            ax.errorbar(
                x,
                r["mag"],
                yerr=r["magerr"],
                fmt=marker,
                ms=6.5,
                color=colors[b],
                markerfacecolor=face,
                capsize=2,
            )
        ax.plot([], [], "o", color=colors[b], label=f"{b} band")
    ax.set(
        xlim=(-0.8, 7.0),
        ylim=(21.2, 17.8),
        ylabel="Pan-STARRS PSF magnitude",
        xticks=[0, 6],
        xticklabels=["Dec 6", "Dec 12"],
    )
    ax.legend(loc="lower right", ncol=2, frameon=False)
    ax.text(
        0.025,
        0.95,
        "Open symbols have partial masks or a defect flag.\nThe square marks the flagged r exposure.",
        transform=ax.transAxes,
        va="top",
        fontsize=8.5,
    )
    for i, b in enumerate(["r", "i"]):
        entry = field["same_image_changes"][b]
        ds = list(entry["control_changes_mag"].values())
        bottom.scatter(
            np.array(ds),
            i + np.linspace(-0.08, 0.08, len(ds)),
            s=38,
            marker="o",
            color="#64748b",
            label="Field stars" if i == 0 else None,
        )
        bottom.scatter(
            entry["target_change_mag"],
            i,
            s=65,
            marker="D",
            color=colors[b],
            label="MGAB-V488" if i == 0 else None,
        )
    bottom.axvline(0, color="#94a3b8", lw=1, zorder=0)
    bottom.set(
        yticks=[0, 1],
        yticklabels=["r", "i"],
        ylim=(-0.5, 1.55),
        xlim=(-3.1, 0.45),
        xlabel="Magnitude change from Dec 6 to 12",
    )
    bottom.legend(loc="upper center", ncol=2, frameon=False, fontsize=8.5)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#e2e8f0", zorder=0, linewidth=0.7)
        axis.set_axisbelow(True)
    fig.savefig(HERE / "december-2013.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
