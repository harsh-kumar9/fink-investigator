"""Arrange archived Fink cutouts without changing their respective image stretches."""

import json
from pathlib import Path
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def panels(source, epochs):
    n = len(epochs)
    fig, axes = plt.subplots(
        n, 3, figsize=(4.45, 1.45 * n + 0.5), layout="constrained", squeeze=False
    )
    for i, (folder, e) in enumerate(epochs):
        for j, item in enumerate(e["images"]):
            ax = axes[i, j]
            img = plt.imread(folder / item["file"])
            ax.imshow(img, cmap="gray", interpolation="nearest", vmin=0, vmax=1)
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_title(item["kind"], fontsize=11)
            if j == 0:
                sign = "+" if str(e["isdiffpos"]).lower() in {"t", "1", "true"} else "−"
                ax.set_ylabel(f"{e['date_utc']}  {sign}", fontsize=9)
    fig.suptitle(source, fontsize=14)
    fig.savefig(HERE / source / "montage.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    for case in json.loads((HERE / "results.json").read_text()):
        source = case["id"]
        epochs = []
        for folder in (HERE / source / "epochs", HERE / source / "negative"):
            if folder.exists():
                epochs.extend(
                    (folder, e)
                    for e in json.loads((folder / "manifest.json").read_text())[
                        "epochs"
                    ]
                )
        panels(source, sorted(epochs, key=lambda item: item[1]["jd"]))
    extra = json.loads((HERE / "offset-results.json").read_text())
    source = extra["id"]
    epochs = []
    for receipt in (HERE / source).glob("*/manifest.json"):
        epochs.extend(
            (receipt.parent, e) for e in json.loads(receipt.read_text())["epochs"]
        )
    panels(source, sorted(epochs, key=lambda item: item[1]["jd"]))
