"""Plot the transferred historical waveform amplitude by phase and temporal sample."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def main(out):
    d = json.loads((HERE / "template.json").read_text())
    primary = d["full_pre_anchor"][0]
    f, (ax, bx) = plt.subplots(
        2, 1, figsize=(6.5, 6), gridspec_kw={"height_ratios": [1.4, 1]}
    )
    f.subplots_adjust(left=0.15, right=0.97, top=0.83, bottom=0.11, hspace=0.64)
    f.text(0.15, 0.965, "MGAB-V488", fontsize=20, weight="bold", ha="left", va="top")
    f.text(
        0.15,
        0.902,
        "Does the published half-day waveform repeat in ZTF?",
        fontsize=11,
        ha="left",
        va="top",
    )
    xx = primary["phase_shifts"] + [1]
    y = primary["amplitude_at_phase"]
    lo, hi = primary["pointwise_bootstrap_p025_p975"]
    blue = "#115f86"
    ax.fill_between(
        xx,
        lo + [lo[0]],
        hi + [hi[0]],
        color="#b4d8e8",
        alpha=0.85,
        lw=0,
        label="95% across resampled blocks",
    )
    ax.plot(xx, y + [y[0]], c=blue, lw=2.3, label="Fitted waveform scale")
    ax.axhline(0, color="#666666", lw=0.6)
    ax.axhline(1, color="#ad481b", ls="--", lw=1.5)
    ax.text(0.03, 1.05, "Published amplitude", color="#943e19", fontsize=9)
    ax.set(
        xlim=(0, 1),
        ylim=(-0.75, 1.27),
        xlabel="Shift from published phase (cycles)",
        ylabel="Amplitude multiplier",
    )
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    widths = [7, 14, 28]
    xx = np.arange(3)
    for key, label, offset, color, marker in [
        ("full_pre_anchor", "All pre-anchor dates", -0.07, blue, "o"),
        ("later_only", "Later dates only", 0.07, "#8c4e89", "s"),
    ]:
        rr = {r["local_bin_width_days"]: r for r in d[key]}
        bx.scatter(
            xx + offset,
            [rr[w]["bootstrap_p95_maximum"] for w in widths],
            s=40,
            c=color,
            marker=marker,
            label=label,
            zorder=3,
        )
    bx.axhline(1, color="#ad481b", ls="--", lw=1.5)
    bx.set(
        xticks=xx,
        xticklabels=[str(w) for w in widths],
        xlim=(-0.47, 2.47),
        ylim=(0, 1.4),
        xlabel="Days per local brightness offset",
        ylabel="95th percentile of\nmaximum scale",
    )
    bx.legend(fontsize=8, loc="lower center", ncols=2, frameon=False)
    for axis in (ax, bx):
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", lw=0.5, color="#e0e0e0")
        axis.set_axisbelow(True)
    f.savefig(out, dpi=175)
    plt.close(f)


def state_figure(out):
    data = json.loads((HERE / "pixels.json").read_text())["apertures_pixels"]["2.2"]
    f, (ax, bx) = plt.subplots(
        2, 1, figsize=(6.5, 5.7), gridspec_kw={"height_ratios": [1.25, 1]}
    )
    f.subplots_adjust(left=0.14, right=0.96, top=0.81, bottom=0.10, hspace=0.70)
    f.text(0.14, 0.96, "ZTF17aabutiy", fontsize=19, weight="bold", ha="left", va="top")
    f.text(
        0.14,
        0.90,
        "What changes when the missing visits are measured?",
        fontsize=11,
        ha="left",
        va="top",
    )
    samples = [
        ("catalogued_flux", "Catalog flux"),
        ("forced_catalogued", "Aperture,\ndetected dates"),
        ("forced_all", "Aperture,\nall dates"),
    ]
    colors = {"r": "#a54c26", "g": "#23765c"}
    for b, shift in [("r", -0.16), ("g", 0.16)]:
        yy = [
            data[b]["matched_bin_medians"][k]["equivalent_dimming_mag"]
            for k, label in samples
        ]
        ax.bar(
            np.arange(len(samples)) + shift,
            yy,
            width=0.30,
            color=colors[b],
            label=f"ZTF {b}",
        )
        for i, y in enumerate(yy):
            ax.text(
                i + shift,
                y + 0.045,
                f"{y:.2f}",
                ha="center",
                fontsize=9,
                color=colors[b],
            )
        groups = [data[b]["windows"][k] for k in ("before", "after")]
        fb = [v["bright_fraction"] * 100 for v in groups]
        bx.plot([0, 1], fb, marker="o", color=colors[b], lw=1.5)
        for i, (frac, v) in enumerate(zip(fb, groups)):
            bx.text(
                i,
                frac + (1.4 if b == "r" else -3),
                f"{round(v['n']*frac/100)}/{v['n']} visits",
                color=colors[b],
                ha="center",
                fontsize=9,
            )
    ax.legend(fontsize=9, frameon=False, loc="upper left", ncols=2)
    ax.set(
        xticks=np.arange(len(samples)),
        xticklabels=[label for k, label in samples],
        ylim=(0, 2.42),
        ylabel="Binned dimming (mag)",
    )
    bx.set(
        xticks=[0, 1],
        xticklabels=["Earlier 180 days", "Later 180 days"],
        ylabel="Bright visits (%)",
        xlim=(-0.23, 1.23),
        ylim=(8, 34),
    )
    for axis in (ax, bx):
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", lw=0.5, color="#e0e0e0")
        axis.set_axisbelow(True)
    f.savefig(out, dpi=175)
    plt.close(f)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--state",
        action="store_true",
        help="Plot the fixed-calendar ZTF17aabutiy measurement",
    )
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    if args.state:
        state_figure(args.out or HERE / "field-visits.png")
    else:
        main(args.out or HERE / "waveform.png")
