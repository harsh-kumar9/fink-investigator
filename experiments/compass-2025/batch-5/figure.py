"""Draw the measured evolution of SN 2025xza, including the late image checks."""

from datetime import datetime, timedelta
import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
COLORS = {"g": "#14766e", "r": "#a24229", "i": "#54468e"}


def main():
    curve = json.loads((HERE / "results/light-curve.json").read_text())
    images = json.loads((HERE / "results/field-calibrated.json").read_text())
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(6.3, 5.7),
        gridspec_kw={"height_ratios": [1.15, 1]},
        layout="constrained",
    )
    a, b = axes
    for band, fid in [("g", 1), ("r", 2)]:
        rows = [row for row in curve["nights"] if row["fid"] == fid]
        a.errorbar(
            [datetime(1858, 11, 17) + timedelta(days=r["mjd"]) for r in rows],
            [r["flux"] for r in rows],
            yerr=[r["error"] for r in rows],
            fmt="o" if band == "g" else "s",
            ms=3.7,
            color=COLORS[band],
            alpha=0.88,
            elinewidth=0.8,
            label=f"ZTF {band}",
        )
    a.set_title("SN 2025xza reddened as it faded", loc="left", weight="bold")
    a.legend(frameon=False, ncol=2, loc="upper right")
    a.set_ylabel("Difference flux (µJy)")
    a.xaxis.set_major_locator(MonthLocator(interval=1))
    a.xaxis.set_major_formatter(DateFormatter("%b"))
    for band, win, label, x in [
        ("r", "dec", "Dec\n2025", 0),
        ("r", "winter", "Jan\n2026", 1),
        ("i", "march", "Mar", 2),
        ("i", "summer", "Jul–Aug", 3),
    ]:
        w = images[band]["comparisons"][win]
        if "delta_ujy" in w["matched"]:
            y, err = w["matched"]["delta_ujy"], w["matched"]["error_ujy"]
        else:
            val = w["nightly"][0]["differences"]["source"]
            y, err = val["excess"], val["error"]
        b.errorbar(x, y, yerr=err, fmt="o", ms=7, capsize=4, color=COLORS[band])
        b.annotate(
            band,
            (x, y),
            xytext=(10, 0),
            textcoords="offset points",
            va="center",
            color=COLORS[band],
            weight="bold",
        )
    b.set_title("Checks on the original survey images", loc="left", weight="bold")
    b.set_ylabel("Flux above pre-event\nbaseline (µJy)")
    b.set_xticks(range(4), ["Dec\n2025", "Jan\n2026", "Mar", "Jul–Aug"])
    b.set_xlim(-0.45, 3.5)
    for ax in axes:
        ax.axhline(0, color="#626b70", lw=0.9, ls="--", zorder=0)
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["right", "top"]].set_visible(False)
        ax.yaxis.set_major_locator(MaxNLocator(5))
    fig.savefig(HERE / "results/supernova.png", dpi=165)
    plt.close(fig)


if __name__ == "__main__":
    main()
