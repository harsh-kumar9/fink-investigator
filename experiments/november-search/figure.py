"""Draw annual ZTF brightness distributions from the replayed measurements."""

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
result = json.loads((HERE / "results/brightening.json").read_text())
fig, ax = plt.subplots(figsize=(4.6, 3.3), layout="constrained")
for label, series, dx, color, marker in (
    ("Red", result["exploratory_red"]["annual_medians"], -0.09, "#ac3323", "o"),
    (
        "Green, field-corrected",
        result["annual_corrected_medians"],
        0.09,
        "#146956",
        "s",
    ),
):
    # A positive offset in brightness corresponds to a smaller apparent magnitude.
    baseline = series["2021"]["mag"]
    positions, values, lower, upper = [], [], [], []
    for year, s in series.items():
        q1, q3 = s["interquartile_mag"]
        positions.append(int(year) + dx)
        values.append(baseline - s["mag"])
        lower.append(q3 - s["mag"])
        upper.append(s["mag"] - q1)
    ax.errorbar(
        positions,
        values,
        yerr=(lower, upper),
        fmt=marker,
        color=color,
        markersize=5,
        capsize=3.5,
        linewidth=1.3,
        label=label,
    )
ax.axhline(0, color="#8b8b8b", linewidth=0.8, zorder=0)
ax.set(
    xlim=(2020.62, 2025.38),
    xticks=range(2021, 2026),
    ylabel="Brightening from 2021 median (mag)",
    xlabel="Year of observation",
)
ax.set_title("ZTF18abqclig", loc="left", weight="bold")
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(labelsize=9)
fig.savefig(HERE / "results/brightening.png", dpi=185)
