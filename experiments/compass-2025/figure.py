"""Show how the astrometric shift depends on the reduction of the same exposures."""

import json
from pathlib import Path
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def draw():
    data = json.loads((HERE / "astrometry/reference.json").read_text())[
        "same_exposure_positions"
    ]["held_out"]
    pairs = data["paired_projections"]
    fig, ax = plt.subplots(figsize=(5.7, 3.8), layout="constrained")
    for i, r in enumerate(pairs):
        ax.plot(
            [i, i],
            [r["direct_arcsec"], r["difference_arcsec"]],
            color="#a5a8ad",
            lw=1,
            zorder=1,
        )
    for key, color, marker, label in [
        ("difference_arcsec", "#a54026", "o", "Difference images"),
        ("direct_arcsec", "#176778", "s", "Direct images"),
    ]:
        ax.scatter(
            range(len(pairs)),
            [r[key] for r in pairs],
            s=42,
            c=color,
            label=label,
            marker=marker,
            zorder=2,
        )
    ax.axhline(0, lw=1, color="#606468", ls="--")
    ax.set_xlabel("Held-out date pair, in time order")
    ax.set_ylabel("Shift along predicted direction (arcsec)")
    ax.set_xticks(range(len(pairs)), range(1, len(pairs) + 1))
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_title("Same exposures, two position measurements", loc="left", fontsize=13)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    fig.savefig(HERE / "results/position-check.png", dpi=180)


if __name__ == "__main__":
    draw()
