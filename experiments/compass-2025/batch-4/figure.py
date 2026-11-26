"""Plot the transient's signed flux and its measured color change against the control sample."""

from datetime import datetime, timedelta
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"


def main():
    light = json.loads((OUT / "onset.json").read_text())
    cmp = json.loads((OUT / "comparison.json").read_text())
    fig, (lc, shape) = plt.subplots(2, 1, figsize=(5.75, 7.8), height_ratios=(1, 1))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.94, bottom=0.16, hspace=0.50)
    lc.set_title(
        "AT 2025iac stayed bright longer in red", fontsize=14, loc="left", pad=10
    )
    for f, label, marker, color in [
        (1, "ZTF g", "o", "#008274"),
        (2, "ZTF r", "^", "#b3471a"),
    ]:
        v = [r for r in light["nightly_forced"] if r["fid"] == f]
        times = [datetime(1858, 11, 17) + timedelta(days=r["mjd"]) for r in v]
        lc.errorbar(
            times,
            [r["flux"] for r in v],
            yerr=[r["error"] for r in v],
            fmt=marker,
            color=color,
            markersize=4,
            capsize=1.4,
            elinewidth=0.8,
            label=label,
        )
    lc.xaxis.set_major_locator(mdates.WeekdayLocator(interval=3))
    lc.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    lc.set_ylabel("Difference flux (µJy)")
    lc.axhline(0, color="#919499", lw=0.8, zorder=0)
    lc.legend(frameon=False, loc="upper left", ncol=2, fontsize=9)
    lc.set_xlabel("Observation date in 2025 (UTC)")
    shape.set_title(
        "A comparison with classified transients", fontsize=14, loc="left", pad=10
    )
    points = {}
    for panel in cmp["panels"].values():
        for r in panel["controls"]:
            q = r["all"]
            if not q:
                continue
            family = r["family"]
            if family in ("IIb", "Ib/c"):
                family = "IIb or Ib/c"
            if family.startswith("SLSN"):
                family = "Superluminous"
            points.setdefault(family, []).append((q["decline"], q["reddening"]))
    palette = [
        ("Ia", "#5076b2", "o"),
        ("II", "#c05a38", "s"),
        ("IIn", "#765d3e", "v"),
        ("IIb or Ib/c", "#818487", "^"),
        ("TDE", "#7d589d", "D"),
        ("Superluminous", "#007d79", "P"),
    ]
    for name, col, m in palette:
        pts = points.get(name)
        if pts:
            x, y = zip(*pts)
            shape.scatter(
                x, y, c=col, marker=m, label=f"{name} ({len(pts)})", s=32, alpha=0.78
            )
    target = cmp["target"]
    x = target["decline"]
    y = target["reddening"]
    shape.scatter([x], [y], c="#111820", marker="*", s=180, zorder=5)
    shape.annotate(
        "AT 2025iac",
        (x, y),
        xytext=(0.06, 0.94),
        textcoords="axes fraction",
        fontsize=9,
        weight="bold",
        arrowprops=dict(arrowstyle="-", lw=0.9, color="#243039"),
    )
    shape.set(
        xlabel="Red fading 25–40 days after red peak (mag)",
        ylabel="Change in g − r color (mag)",
    )
    shape.margins(x=0.09, y=0.16)
    shape.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=3,
        frameon=False,
        fontsize=8.7,
        columnspacing=1.1,
        handletextpad=0.15,
    )
    for ax in (lc, shape):
        ax.tick_params(labelsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", linewidth=0.6, color="#dadde1")
        ax.set_axisbelow(True)
    fig.savefig(OUT / "transient.png", dpi=165)
    plt.close(fig)


if __name__ == "__main__":
    main()
