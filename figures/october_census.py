"""Render the two archived census figures. Requires Matplotlib and Pillow."""

import json
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import census
import census_photometry
import pilot

BG = "#faf8f2"
INK = "#193840"
GREEN = "#006b64"
ORANGE = "#aa4328"


def images():
    folders = sorted((census.CENSUS / "review").glob("*/images"))
    fig, axes = plt.subplots(
        len(folders), 3, figsize=(5.25, 1.7 * len(folders) + 0.65), facecolor=BG
    )
    fig.subplots_adjust(
        left=0.06, right=0.96, top=0.92, bottom=0.014, hspace=0.75, wspace=0.08
    )
    fig.suptitle(
        "Repeat observations of each alert position",
        fontsize=15.4,
        weight="bold",
        color=INK,
        y=0.99,
    )
    for row, folder in zip(axes, folders):
        meta = json.loads((folder / "manifest.json").read_text())
        for a in row:
            a.axis("off")
        for a, e in zip(row, meta["epochs"]):
            item = next(i for i in e["images"] if i["kind"] == "Difference")
            with Image.open(folder / item["file"]) as im:
                a.imshow(im, cmap="gray", vmin=0, vmax=255, interpolation="nearest")
            a.set_title(
                e["date_utc"] + " " + ("+" if e["isdiffpos"] in ("t", "1") else "−"),
                fontsize=8.6,
                color=INK,
                pad=3,
            )
        row[0].text(
            0,
            1.30,
            meta["id"],
            size=10,
            weight="bold",
            color=INK,
            transform=row[0].transAxes,
        )
    fig.savefig(census.CENSUS / "review/differences.png", dpi=140)
    plt.close(fig)


def photometry():
    previous, cases = census.assess(False)
    _, anchors, _ = census.locked_cohort()
    _, raw = pilot.snapshot(census.BASE / "snapshot")
    old = pilot.histories(
        [r for n, v in raw.items() if n.startswith("history-") for r in v]
    )
    new = pilot.histories(census.queries("fink", anchors, False))
    control = next(r for r in previous["objects"] if r["id"] == "ZTF18achfhai")
    target = next(r for r in cases if r["id"] == "ZTF17aacydyv")
    fig, axes = plt.subplots(2, 1, figsize=(5.3, 6.1), facecolor=BG)
    fig.subplots_adjust(left=0.165, right=0.945, top=0.86, bottom=0.135, hspace=0.62)
    fig.text(
        0.08,
        0.968,
        "Two ways to measure brightness",
        va="top",
        fontsize=15.6,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.08,
        0.915,
        "One dot per matched UTC date in ZTF’s red filter",
        fontsize=11.1,
        color=INK,
    )
    for a, case, h, tag in zip(
        axes,
        [control, target],
        [old, new],
        ["Known RR Lyrae control", "Unclassified ZTF position"],
    ):
        data = census_photometry.compare(case, h, False, include_pairs=True)
        pts = data["paired_medians"]
        x = [p[0] for p in pts]
        y = [20 - 2.5 * math.log10(p[1]) for p in pts]
        a.scatter(
            x,
            y,
            c=[GREEN if n > 0 else ORANGE for n in x],
            s=21,
            alpha=0.77,
            linewidths=0,
            zorder=3,
        )
        a.axvline(0, color="#819296", lw=1)
        a.invert_yaxis()
        a.grid(axis="y", color="#dce1df")
        a.set_axisbelow(True)
        a.set_facecolor(BG)
        for side in ["right", "top"]:
            a.spines[side].set_visible(False)
        a.set_ylabel("Direct magnitude\n(brighter above)", fontsize=10.6)
        a.set_title(
            f"{tag} · {case['id']}",
            loc="left",
            fontsize=12.5,
            color=INK,
            weight="bold",
            pad=7,
        )
        a.text(
            0.025,
            0.055,
            f"ρ = {data['spearman']:.2f}  ·  {len(pts)} dates",
            transform=a.transAxes,
            fontsize=10,
            bbox={"facecolor": BG, "edgecolor": "none", "alpha": 0.9},
        )
    axes[-1].set_xlabel("Signed difference flux (relative units)", fontsize=11)
    fig.text(
        0.08,
        0.035,
        "Both reductions use the same telescope exposures.",
        fontsize=10,
        color=INK,
    )
    fig.savefig(census.CENSUS / "results/photometry.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    images()
    photometry()
