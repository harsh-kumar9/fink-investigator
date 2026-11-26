"""Plot measured total brightness near the two J0722 events."""

from datetime import datetime, timezone, timedelta
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from onset import NAME, HERE, NMW, direct

MJD = datetime(1858, 11, 17, tzinfo=timezone.utc)
COLORS = {"g": "#08766e", "r": "#a53b25"}


def draw():
    results = json.loads((HERE / "results/onset.json").read_text())
    broker = json.loads((HERE / "results/broker.json").read_text())
    ledger = json.loads((HERE / "results/ledger.json").read_text())
    (case,) = [r for r in ledger["objects"] if r["id"] == NAME]
    _, rows, _ = direct.direct_data(case, HERE / "direct" / NAME)
    fig, axes = plt.subplots(2, 1, figsize=(7, 6.7), layout="constrained")
    ranges = [
        ("2025-01-09", "2025-01-22", "Before the discovery"),
        ("2025-04-01", "2025-05-05", "A late brightening during the decline"),
    ]
    for ax, (lo, hi, title) in zip(axes, ranges):
        start, stop = [
            datetime.fromisoformat(d).replace(tzinfo=timezone.utc) for d in (lo, hi)
        ]
        for b, color in COLORS.items():
            ticks = []
            for row in rows:
                t = MJD + timedelta(days=float(row["mjd"]))
                if (
                    row["filtercode"] == "z" + b
                    and not int(row["catflags"])
                    and 0 < float(row["magerr"]) <= 0.2
                    and start <= t < stop
                ):
                    ticks.append((t, float(row["mag"]), float(row["magerr"])))
            x, y, err = zip(*ticks)
            ax.errorbar(
                x,
                y,
                yerr=err,
                fmt="o" if b == "g" else "^",
                ms=6,
                capsize=2,
                color=color,
                alpha=0.95,
            )
        ax.invert_yaxis()
        ax.set_xlim(start, stop)
        ax.xaxis.set_major_locator(
            mdates.DayLocator(interval=3 if ax == axes[0] else 7)
        )
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.set_title(title, loc="left", weight="bold", fontsize=13)
        ax.set_ylabel("Apparent magnitude")
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["right", "top"]].set_visible(False)
    ax = axes[0]
    for b, color in COLORS.items():
        ax.axhline(
            results["bands"][b]["raw"]["january"]["median"],
            c=color,
            alpha=0.5,
            ls=":",
            lw=1.5,
        )
        for d in broker["dates"].values():
            a = d[b]
            ax.scatter(
                MJD + timedelta(days=a["mjd"]),
                a["forced_mag"],
                marker="s",
                s=67,
                facecolors="none",
                edgecolors=color,
                lw=1.3,
                zorder=5,
            )
    ax.axvline(MJD + timedelta(days=NMW - 2400000.5), color="#444", ls="--", lw=1)
    ax.text(
        MJD + timedelta(days=NMW - 2400000.5 - 0.13),
        0.96,
        "Discovery",
        ha="right",
        va="top",
        transform=ax.get_xaxis_transform(),
        fontsize=9,
    )
    axes[1].annotate(
        "April 16",
        xy=(datetime(2025, 4, 16, 4, tzinfo=timezone.utc), 16.42),
        xytext=(12, 1),
        textcoords="offset points",
        fontsize=10,
    )
    keys = [
        Line2D([0], [0], marker=m, label=n, color=c, ls="none")
        for n, m, c in [("ZTF g", "o", COLORS["g"]), ("ZTF r", "^", COLORS["r"])]
    ]
    keys.extend(
        [
            Line2D(
                [0],
                [0],
                marker="s",
                label="Forced flux",
                mfc="none",
                color="#555",
                ls="none",
            ),
            Line2D([0], [0], ls=":", c="#666", label="Pre-event median"),
        ]
    )
    fig.legend(
        handles=keys, loc="outside lower center", ncol=2, frameon=False, fontsize=10
    )
    fig.savefig(HERE / "events.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    draw()
