"""Regenerate the public VB 10 figure from the archived measurements and ZTF PNGs.

Optional plotting environment: Python 3.10+, matplotlib and Pillow.
Run from the repository root with python figures/vb10.py.
"""

import json
from datetime import date
import gzip
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pilot

INK = "#17343d"
SECONDARY = "#455c60"
TEAL = "#14786d"
BURNT = "#a6411e"
PAPER = "#faf8f2"


def inputs():
    case = next(
        r
        for r in pilot.evaluate(ROOT / "data/pilot")["objects"]
        if r["id"] == "ZTF18abgobyb"
    )
    _, raw = pilot.snapshot(ROOT / "data/pilot")
    source_id = case["gaia_matches"][0]["source_id"]
    gaia = next(
        s
        for key, rows in raw.items()
        if key.startswith("gaia-")
        for s in rows
        if s["source_id"] == source_id
    )
    x, y = pilot.offsets(case["ra"], case["dec"], float(gaia["ra"]), float(gaia["dec"]))
    vx, vy = float(gaia["pmra"]) / 1000, float(gaia["pmdec"]) / 1000
    epoch = float(gaia["ref_epoch"]) + (x * vx + y * vy) / (vx * vx + vy * vy)
    crossing = pilot.day_of(2451545.0 + (epoch - 2000) * 365.25)
    references = json.loads(
        gzip.decompress((ROOT / "results/moving-star/reference.json.gz").read_bytes())
    )
    reference = next(r for r in references if str(r["i:candid"]) == case["candid"])
    start, end = [pilot.day_of(reference[k]) for k in ("i:jdstartref", "i:jdendref")]
    return case, crossing, reference["i:nframesref"], start, end


def draw():
    case, crossing, count, first, last = inputs()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "text.color": INK,
            "savefig.facecolor": PAPER,
            "figure.facecolor": PAPER,
        }
    )
    # Explicit page coordinates keep labels readable when the README shrinks to phone width.
    fig = plt.figure(figsize=(6.3, 6.85), dpi=160)

    def label(x, y, value, size, **kwargs):
        fig.text(x, y, value, fontsize=size, va="top", **kwargs)

    label(0.055, 0.966, "An alert at VB 10’s old position", 21.5, weight="bold")
    label(
        0.055,
        0.899,
        "ZTF compared a new exposure with an older image\nof the same patch of sky.",
        14.5,
        linespacing=1.3,
        color=SECONDARY,
    )
    columns = [
        ("template", "Reference", "Spring 2018"),
        ("science", "New image", "12 Oct 2024"),
        ("difference", "Difference", "Paired residual"),
    ]
    for i, (kind, title, caption) in enumerate(columns):
        x = 0.055 + i * 0.302
        label(x, 0.783, title, 14.4, weight="bold")
        axis = fig.add_axes([x, 0.482, 0.286, 0.267])
        stamp = ROOT / "review" / case["id"] / (case["candid"] + "-" + kind + ".png")
        with Image.open(stamp) as im:
            axis.imshow(im, cmap="gray", vmin=0, vmax=255, interpolation="nearest")
        axis.set_axis_off()
        label(x, 0.473, caption, 12.8, color=SECONDARY)
    fig.add_artist(
        plt.Line2D(
            [0.055, 0.945],
            [0.407, 0.407],
            transform=fig.transFigure,
            color="#c0cdca",
            lw=1.1,
        )
    )
    label(0.055, 0.376, "When did VB 10 pass this position?", 17.5, weight="bold")
    # Calendar scale deliberately uses the full interval between the archive endpoints.
    axis = fig.add_axes([0.095, 0.158, 0.805, 0.14], facecolor=PAPER)
    axis.set_axis_off()
    axis.set_xlim(date(2018, 3, 1), date(2018, 7, 15))
    axis.set_ylim(-1, 1)
    axis.plot([date(2018, 3, 1), date(2018, 7, 15)], [0, 0], lw=1.5, color="#9daaa9")
    axis.plot([first, last], [0, 0], color=TEAL, lw=8, solid_capstyle="butt")
    for d in (first, last):
        axis.plot([d, d], [-0.10, 0.10], color=TEAL, lw=2.2)
        axis.text(
            d,
            -0.24,
            d.strftime("%d %b").lstrip("0"),
            ha="center",
            va="top",
            size=13,
            color=TEAL,
        )
    axis.plot(
        [crossing],
        [0],
        marker="D",
        markersize=8,
        markerfacecolor=BURNT,
        markeredgecolor=PAPER,
        markeredgewidth=1.3,
    )
    axis.text(
        crossing,
        0.28,
        crossing.strftime("%d %b") + " estimate",
        ha="center",
        va="bottom",
        size=14,
        color=BURNT,
        weight="bold",
    )
    label(
        0.055,
        0.129,
        f"The inferred crossing falls within the dates of the\n{count} exposures used to build the reference image.",
        14,
        linespacing=1.3,
    )
    label(
        0.055,
        0.044,
        "ZTF images via Fink  ·  Gaia DR3 astrometry  ·  Dates in UTC",
        10.8,
        color=SECONDARY,
    )
    fig.savefig(ROOT / "figures/vb10.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    draw()
