"""Plot individual date-pair vectors from the selected sign-dependent astrometry packet."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = Path(__file__).resolve().parent


def draw():
    packet = json.loads((HERE / "packet.json").read_text())
    case = packet["case"]
    band = case["best_band"]
    stats = case["bands"][band]
    fig, ax = plt.subplots(figsize=(4.9, 4.6), layout="constrained")
    ax.axvline(0, color="#c7cbd0", lw=1, zorder=0)
    ax.axhline(0, color="#c7cbd0", lw=1, zorder=0)
    for key, prefix, color, marker, label in (
        ("training", "training", "#536b8e", "o", "Earlier"),
        ("held_out", "test", "#bb531b", "^", "Later"),
    ):
        pairs = packet["pairs"][key]
        ax.scatter(
            [p["east"] for p in pairs],
            [p["north"] for p in pairs],
            c=color,
            marker=marker,
            s=27,
            alpha=0.76,
            label=f"{label} pairs ({len(pairs)})",
        )
        dx, dy = stats[prefix + "_east_north_arcsec"]
        ax.annotate(
            "",
            xy=(dx, dy),
            xytext=(0, 0),
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "lw": 2.3,
                "mutation_scale": 15,
            },
        )
    ax.set(
        xlabel="Eastward displacement (arcsec)",
        ylabel="Northward displacement (arcsec)",
        title=f"{case['id']} · {band}-band centroid shifts",
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.xaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", framealpha=0.95, fontsize=9)
    fig.savefig(HERE / "shifts.png", dpi=180)


if __name__ == "__main__":
    draw()
