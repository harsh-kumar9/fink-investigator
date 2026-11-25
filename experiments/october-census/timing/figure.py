"""Plot the later red nights and the green nights missing from both red series."""

import json
import matplotlib.pyplot as plt
import numpy as np

from analyze import CENSUS, HERE, nights
from epochs import Ephemeris, NAME
from green import check as green


def draw():
    fit = Ephemeris()
    dips = json.loads((HERE / "dips.json").read_text())
    box = dips["fits"][dips["training_winner"]]
    red = fit.data[fit.n :]
    cases = json.loads((CENSUS / "results/direct.json").read_text())["exploration"]
    case = next(c for c in cases if c["id"] == NAME)
    red_dates = [np.floor(fit.data[:, 0])]
    for oid in dips["other_series_disjoint_dates"]:
        red_dates.append(np.floor(nights({**case, "direct_source_id": oid})[:, 0]))
    _, stars = green.select(False)
    data = green.green_data(stars, False)[stars[0]["oid"]]
    blue = data[~np.isin(np.floor(data[:, 0]), np.concatenate(red_dates))]
    fig, axes = plt.subplots(
        2, 1, figsize=(5.6, 5.7), sharex=True, layout="constrained"
    )
    h = box["duration_days"] / (2 * box["period_days"])
    for ax, values, color, title, label in zip(
        axes,
        [red, blue],
        ["#943a30", "#197062"],
        [
            f"Red · {len(red)} later nights",
            f"Green · {len(blue)} nights without red data",
        ],
        ["r", "g"],
    ):
        phase = (
            (values[:, 1] - box["center_hjd"]) / box["period_days"] + 0.25
        ) % 1 - 0.25
        ax.axvspan(-h, h, color="#e9d8ae", alpha=0.75, label="Predicted deep fade")
        ax.axvspan(
            0.5 - h, 0.5 + h, color="#e3e8ef", alpha=0.9, label="Half a cycle later"
        )
        ax.scatter(phase, values[:, 2], s=16, color=color, alpha=0.78, linewidths=0)
        ax.set(ylabel=f"{label} magnitude", xlim=(-0.25, 0.75))
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=13)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(loc="lower right", fontsize=10, framealpha=0.95)
    axes[1].set_xticks([-0.25, 0, 0.25, 0.5, 0.75])
    axes[1].set_xlabel("Phase of the 11.62-day cycle")
    fig.savefig(HERE / "eclipse.png", dpi=180)


if __name__ == "__main__":
    draw()
