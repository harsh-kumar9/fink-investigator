"""Read plotted marker centers from Tarasenkov et al.'s public Figure 1 EPS source."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from statistics import median
import tarfile

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
YEAR = BATCH.parent
sys.path.insert(0, str(YEAR.parents[1]))
import compass_archive
import pilot

URL = "https://arxiv.org/src/2505.20842v1"
SHA = "78c54a6c244bba3f179aff6bc50ded13df8618cd5af4a763f9886a7a71aa2664"
T0 = 2460696.4416
TARGET = "ZTF25aacfjde"


def ticks(text, horizontal):
    # Labels' baselines at y=960 for x ticks or x=930 for y ticks.
    command = r"^(\d+) 960 M" if horizontal else r"^930 (\d+) M"
    labels = re.findall(command + r"\n\[ \[\(Times-Roman\)[^\n]*?\(\s*(-?\d+)\)\]", text, re.M)
    points = [(int(coord), int(label)) for coord, label in labels]
    if len(points) < 2:
        raise ValueError("Tick labels missing")
    (x0, y0), (x1, y1) = points[0], points[-1]
    scale = (y1 - y0) / (x1 - x0)
    intercept = y0 - x0 * scale
    if any(abs((coord * scale + intercept) - label) > abs(scale) for coord, label in points):
        raise ValueError("Inconsistent tick spacing")
    return scale, intercept


def panel(archive, filename):
    eps = archive.extractfile(filename).read()
    body = eps.decode("latin-1").split("%%Page: 1 1", 1)[1]
    sx, ix = ticks(body, True)
    sy, iy = ticks(body, False)
    layers = {}
    for idx, content in re.findall(r"% Begin plot #(\d+)\n(.*?)% End plot #\1\n", body, re.S):
        points = re.findall(r"^(\d+) (\d+) (CircleF|TriD)$", content, re.M)
        # Keyed layers draw their final circle at the legend's sample and matching y label.
        legend = re.search(r"(\d+) (\d+) M\n\[ \[\(Times-Roman\).*?\]\n\] -100.0 MRshow", content)
        if legend:
            x, y, kind = points.pop()
            if kind != "CircleF" or y != legend[2] or int(x) - int(legend[1]) != 571:
                raise ValueError("Unexpected legend symbol")
        layers[int(idx)] = [(sx * int(x) + ix, sy * int(y) + iy, kind) for x, y, kind in points]
    return layers, {"filename": filename, "sha256": hashlib.sha256(eps).hexdigest(),
                    "days_per_horizontal_unit": round(sx, 6),
                    "magnitudes_per_vertical_unit": round(abs(sy), 6),
                    "latest_marker_days": round(max(day for layer in layers.values()
                                                     for day, _, _ in layer), 3)}


def atlas_brackets(top, middle, exposure_days):
    # Both panels use cyan for layer 4 and orange for layer 5. Prefer the zoom when covered.
    result = {}
    for band, layer in [("c", 4), ("o", 5)]:
        for side in ("before", "after"):
            for panel_name, panel_data in [("middle", middle), ("full", top)]:
                data = [(day, mag) for day, mag, marker in panel_data[layer] if marker == "CircleF"
                        and (day < min(exposure_days) if side == "before" else day > max(exposure_days))]
                if data:
                    break
            nearest = (max if side == "before" else min)(d for d, _ in data)
            group = [(d, m) for d, m in data if abs(d - nearest) < 0.5]
            days, mags = zip(*group)
            dec = 2 if panel_name == "middle" else 1
            result[band + "_" + side] = {"panel": panel_name, "markers_within_half_day": len(group),
                                         "days_range": [round(min(days), dec), round(max(days), dec)],
                                         "magnitude_range": [round(min(mags), 2), round(max(mags), 2)],
                                         "median_mag": round(median(mags), 2)}
    return result


def check(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA:
        raise ValueError("arXiv source checksum differs")
    with tarfile.open(path) as tar:
        top, top_info = panel(tar, "j0722_combined_lightcurve_all.eps")
        middle, middle_info = panel(tar, "j0722_combined_lightcurve.eps")
    if sorted(top) != list(range(1, 13)) or sorted(middle) != list(range(1, 11)):
        raise ValueError("Unexpected number of series")
    _, histories, _ = compass_archive.load(YEAR / "snapshot")
    direct = json.loads((BATCH / "images/rebrightenings.json").read_text())
    january = []
    overlaps = {}
    for band, layer, info in [("g", middle[6], middle_info), ("r", top[6], top_info)]:
        fid = 1 if band == "g" else 2
        alerts = [r for r in histories[TARGET] if r["i:fid"] == fid]
        matched = []
        for day, mag, kind in layer:
            # The coordinate quantization and the rounded axis ticks each contribute.
            # A unique single-filter time match within one horizontal grid unit is required.
            rows = [r for r in alerts if abs(r["i:jd"] - T0 - day) <= info["days_per_horizontal_unit"]]
            if kind == "CircleF" and len(rows) == 1 and abs(mag - rows[0]["i:magpsf"]) <= 0.005:
                row = rows[0]
                matched.append(str(row["i:candid"]))
                if str(pilot.day_of(row["i:jd"])) == "2025-01-19":
                    pos = next(c for c in direct if c["candid"] == str(row["i:candid"]))
                    total = pos["direct"]["magnitude"]
                    january.append({"band": band, "candid": str(row["i:candid"]),
                                    "figure_day": round(day, 3), "exposure_day": round(row["i:jd"]-T0, 4),
                                    "figure_magnitude": round(mag, 3),
                                    "difference_magnitude": row["i:magpsf"],
                                    "direct_magnitude": total,
                                    "figure_minus_direct_mag": round(mag - total, 3)})
        if len(set(matched)) != len(matched):
            raise ValueError("Reused alert exposure")
        overlaps[band] = {"plotted_markers": len(layer), "matching_alerts": len(matched)}
    if len(january) != 2:
        raise ValueError("Both January 19 filters expected")
    limits = [m for d,m,k in middle[8] if k == "TriD" and -1.942 <= d < -0.942]
    return {"arxiv_source": URL, "source_sha256": SHA, "discovery_jd": T0,
            "panels": [top_info, middle_info], "january19": january, "ztf_series_comparison": overlaps,
            "atlas_bracketing_january19": atlas_brackets(top, middle, [r["exposure_day"] for r in january]),
            "january19_ASASSN_upper_limits": {"count": len(limits), "magnitude_range":
                                                [round(min(limits), 2), round(max(limits), 2)]}}


def verify_journal(path, result):
    import pymupdf

    if hashlib.sha256(path.read_bytes()).hexdigest() != "557d6c98c3261aeb065e13f3615ec2fa3b08e771c95b673ddc03be48b40d6c79":
        raise ValueError("Publisher PDF checksum differs")
    with pymupdf.open(path) as doc:
        page = next(p for p in doc if "upper limits (triangles)" in p.get_text())
        drawings = page.get_drawings()
    for r in result["january19"]:
        # Independently checked axis ticks in final journal page 4, in PDF points from top left.
        if r["band"] == "g":
            xmin, xmax, y12, color = -20, 60, 229.266, (0.071, 0.161, 0.290)
        else:
            xmin, xmax, y12, color = -3500, 500, 82.302, (0.553, 0.0, 0.0)
        x = 115.529 + (r["figure_day"] - xmin) / (xmax - xmin) * (514.411 - 115.529)
        y = y12 + (r["figure_magnitude"] - 12) * 8.588
        matches = []
        for shape in drawings:
            if shape["type"] != "f" or tuple(round(c, 3) for c in shape["fill"]) != color:
                continue
            b = shape["rect"]
            if abs((b.x0 + b.x1) / 2 - x) < 0.1 and abs((b.y0 + b.y1) / 2 - y) < 0.1:
                matches.append(b)
        if len(matches) != 1:
            raise ValueError("Journal marker missing or ambiguous")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source_archive", type=Path)
    p.add_argument("--record", action="store_true")
    p.add_argument("--publisher-pdf", type=Path, help="Also verify final journal markers with PyMuPDF")
    args = p.parse_args()
    result = check(args.source_archive)
    file = HERE / "results.json"
    if args.record:
        pilot.json_write(file, result)
    elif json.loads(file.read_text()) != result:
        raise ValueError("Figure or exposure comparison changed")
    if args.publisher_pdf:
        verify_journal(args.publisher_pdf, result)
    print(json.dumps(result, indent=2))
