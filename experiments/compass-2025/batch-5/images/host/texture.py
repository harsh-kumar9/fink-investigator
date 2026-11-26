"""Compare a weak pre-event arm knot with nearby positions of similar brightness."""

import math
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.spatial import cKDTree


def highpass(images, width):
    peaks = {}
    for layer, (data, wcs, pixel) in images.items():
        for band, im in zip("grz", data):
            smooth = gaussian_filter(im, 0.22 / pixel)
            background = gaussian_filter(im, width / pixel)
            residual = smooth - background
            maxima = (
                residual == maximum_filter(residual, 2 * round(0.4 / pixel) + 1)
            ) & (residual > 0)
            y, x = np.nonzero(maxima)
            peaks[(layer, band)] = (
                cKDTree(np.column_stack((x, y))),
                wcs,
                pixel,
                residual,
                background,
            )
    return peaks


def distances(peaks, center):
    values = {}
    for (layer, band), (tree, wcs, pixel, _, _) in peaks.items():
        x, y = map(float, wcs.world_to_pixel_values(*center))
        d, _ = tree.query([x, y])
        values[layer + "_" + band] = float(d * pixel)
    return values


def test(images, center, catalog, offset, distance, plot_path=None):
    output = {}
    nominal = None
    for width in (0.9, 1.5, 2.5):
        peaks = highpass(images, width)
        output[str(width)] = {
            k: round(v, 3) for k, v in distances(peaks, center).items()
        }
        if width == 1.5:
            nominal = peaks
    tree, wcs, pixel, residual, background = nominal[("hsc-dr3", "r")]

    def surface(p):
        x, y = map(float, wcs.world_to_pixel_values(*p))
        return background[round(y), round(x)]

    target_level = surface(center)
    if target_level <= 0:
        raise ValueError("Missing disk light")
    n = hits = 0
    for dx in range(-20, 21, 2):
        for dy in range(-20, 21, 2):
            if not 3 <= math.hypot(dx, dy) <= 20:
                continue
            pos = offset(*center, dx, dy)
            ratio = surface(pos) / target_level
            if not 2 / 3 <= ratio <= 1.5:
                continue
            if min(distance(*pos, r["ra"], r["dec"]) for r in catalog) <= 1:
                continue
            n += 1
            ds = distances(nominal, pos)
            hits += all(v <= 0.5 for key, v in ds.items() if key != "ls-dr10_z")
    if n < 5:
        raise ValueError("Too few local arm controls")
    if plot_path:
        import matplotlib.pyplot as plt

        f, axes = plt.subplots(2, 3, figsize=(6, 4.2), layout="constrained")
        for row, layer in zip(axes, ["hsc-dr3", "ls-dr10"]):
            for ax, b in zip(row, "grz"):
                tree, wcs, pixel, residual, _ = nominal[(layer, b)]
                x, y = map(float, wcs.world_to_pixel_values(*center))
                k = round(5 / pixel)
                xx, yy = round(x), round(y)
                v = residual[yy - k : yy + k + 1, xx - k : xx + k + 1]
                amp = max(abs(np.percentile(v, [3, 97])))
                ax.imshow(v, origin="lower", cmap="RdBu_r", vmin=-amp, vmax=amp)
                ax.plot(k + x - xx, k + y - yy, "k+", ms=7, mew=1)
                ax.set(
                    xticks=[],
                    yticks=[],
                    title=("HSC " if layer == "hsc-dr3" else "DECam ") + b,
                )
        f.suptitle("Local contrast before the transient", fontsize=13)
        f.savefig(plot_path, dpi=130)
        plt.close(f)
    return {
        "nearest_positive_peak_arcsec_by_broad_sigma_arcsec": output,
        "comparison_positions": n,
        "five_band_alignments_at_comparison_positions": hits,
        "target_five_band_alignment": all(
            v <= 0.5 for k, v in distances(nominal, center).items() if k != "ls-dr10_z"
        ),
    }
