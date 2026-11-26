"""Plot measured eruption strengths and the contemporaneous TESS aperture series."""

import argparse
from datetime import date
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from analyze_vars import HERE, load, pilot
from cadence import NAME
from tess_data import SELECTION, bin6, cube, sectors
from tess_events import corrected


def main(out):
    case, pre, _ = load(NAME)
    opts = sectors(case, False)
    groups = json.loads((HERE / "alternating.json").read_text())
    weak = groups["2019_short"]
    strong = groups["2019_intervening"]
    full = groups["all_2019_groups"]
    start = date(2019, 7, 15)
    end = date(2019, 10, 31)
    f, (ax, bx) = plt.subplots(
        2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.5, 1]}
    )
    f.subplots_adjust(left=0.13, right=0.98, top=0.81, bottom=0.11, hspace=0.18)
    f.text(
        0.13,
        0.965,
        "Two kinds of peak in MGAB-V960",
        ha="left",
        va="top",
        fontsize=17,
        weight="bold",
    )
    f.text(
        0.13,
        0.89,
        "ZTF and TESS observed the same 2019 eruption sequence.",
        ha="left",
        fontsize=10.5,
    )
    for b, color in [("g", "#33886e"), ("r", "#a54832")]:
        x = [d for d in pre[b] if start <= d <= end]
        ax.scatter(
            x,
            [pre[b][d][1] for d in x],
            s=12,
            color=color,
            label=f"ZTF {b}",
            alpha=0.78,
            linewidths=0,
        )
    for seq, color, mark, label in [
        (weak, "#644d87", "o", "Short peak"),
        (strong, "#233e65", "D", "Intervening peak"),
    ]:
        ax.scatter(
            [date.fromisoformat(g["r"]["date"]) for g in seq],
            [g["r"]["mag"] for g in seq],
            color="none",
            edgecolor=color,
            marker=mark,
            lw=1.3,
            s=47,
            label=label,
            zorder=5,
        )
    ax.invert_yaxis()
    ax.set_ylabel("Magnitude")
    ax.legend(loc="lower left", frameon=True, fontsize=8, ncols=2, facecolor="white")
    peaks = [
        pilot.utc_jd(date.fromisoformat(g["r"]["date"] if g["r"] else g["first"])) + 0.5
        for g in full
    ]
    for sector, color in zip(SELECTION[NAME], ["#416782", "#46766d", "#75628c"]):
        m, t, p = cube(case, sector, *opts[sector], False)
        arr = bin6(t, p)
        tt = arr[:, 0]
        train = np.ones(len(tt), dtype=bool)
        for d in peaks:
            train &= abs(tt - d) > 1.5
        px = corrected(arr, train)[:, 1:].reshape(-1, 9, 9)
        flux = px[:, 3:6, 3:6].sum(axis=(1, 2))
        flux -= np.median(flux)
        epoch = mdates.date2num(date(1970, 1, 1)) + tt - 2440587.5
        bx.scatter(epoch, flux, s=4, alpha=0.35, lw=0, color=color)
        day = np.floor(epoch)
        unique = np.unique(day)
        bx.plot(
            unique + 0.5,
            [np.median(flux[day == i]) for i in unique],
            lw=1.5,
            color=color,
            label=f"Sector {sector}",
        )
    bx.set_ylabel("TESS aperture\nflux offset (e⁻/s)")
    bx.legend(fontsize=8, ncols=3, frameon=False, loc="upper right")
    bx.set_xlim(start, end)
    bx.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1, 15]))
    bx.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    for a in (ax, bx):
        a.spines[["top", "right"]].set_visible(False)
        a.grid(axis="y", color="#dddddd", lw=0.5)
        a.set_axisbelow(True)
    f.savefig(out, dpi=170)
    plt.close(f)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "mgab-v960.png")
    main(p.parse_args().out)
