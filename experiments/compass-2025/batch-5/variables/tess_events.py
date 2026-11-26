"""Suppress local TESS common trends and measure source-position and timing evidence."""

import argparse
from datetime import date
import json, math
from pathlib import Path

import numpy as np

from analyze_vars import HERE, NAMES, load, pilot
from tess_data import SELECTION, bin6, cube, sectors, previous_tess
from alternating import group_peaks

CORNER = previous_tess.APERTURES


def corrected(arr, training):
    p = arr[:, 1:]
    yy, xx = np.indices((9, 9))
    mask = (np.maximum(abs(yy - 4), abs(xx - 4)) >= 3).ravel()
    comp = p[:, mask]
    center = np.median(comp[training], axis=0)
    scale = 1.4826 * np.median(abs(comp[training] - center), axis=0)
    keep = scale > 0
    norm = (comp[:, keep] - center[keep]) / scale[keep]
    pc = np.linalg.svd(norm[training], full_matrices=False)[2][:4]
    x = np.column_stack([np.ones(len(arr)), norm @ pc.T])
    coef = np.linalg.lstsq(x[training], p[training], rcond=None)[0]
    return np.column_stack([arr[:, 0], p - x @ coef])


def contrast(px, bright, faint, xy):
    if not (bright.any() and faint.any()):
        return None
    dif = np.median(px[bright], axis=0) - np.median(px[faint], axis=0)
    yy, xx = np.unravel_index(dif.argmax(), dif.shape)
    x, y = xy
    ap = {
        k: round(float(dif[r : r + 3, c : c + 3].sum()), 6)
        for k, (r, c) in CORNER.items()
    }
    return {
        "brightest_pixel": [int(xx), int(yy)],
        "distance_pixels": round(math.hypot(x - xx, y - yy), 4),
        "central_pixel": round(float(dif[round(y), round(x)]), 6),
        "anywhere_peak": round(float(dif.max()), 6),
        "apertures": ap,
        "difference_image": np.round(dif, 5).tolist(),
    }


def plateaus(t, ap):
    mid = np.median(ap)
    scale = 1.4826 * np.median(abs(ap - mid))
    bright = ap > mid + 1.5 * scale
    runs = []
    group = []
    for k in np.flatnonzero(bright):
        if group and t[k] - t[group[-1]] > 0.5001:
            runs.append(group)
            group = []
        group.append(k)
    if group:
        runs.append(group)
    return [
        {
            "first_bjd": round(float(t[g[0]]), 5),
            "last_bjd": round(float(t[g[-1]]), 5),
            "sampled_days": round(float(t[g[-1]] - t[g[0]]), 4),
            "bins": len(g),
            "accepted": bool(
                t[g[-1]] - t[g[0]] >= 5
                and len(g) >= 0.7 * (1 + 4 * (t[g[-1]] - t[g[0]]))
            ),
        }
        for g in runs
    ]


def run():
    groups = [g for g in group_peaks() if "2019-07-15" <= g["first"] <= "2019-10-31"]
    out = {}
    for n in NAMES:
        case, pre, _ = load(n)
        options = sectors(case, False)
        selection = groups if n == "ZTF18abcckqd" else []
        if n == "ZTF17aaaoxxi":
            allg = json.loads((HERE / "duration.json").read_text())[n]["tests"][
                "primary"
            ]["groups"]
            selection = [
                {
                    **g,
                    "r": {
                        "date": min(
                            (str(d) for d in pre["r"] if str(d) in g["dates"]),
                            key=lambda d: pre["r"][date.fromisoformat(d)][1],
                            default=g["first"],
                        )
                    },
                }
                for g in allg
            ]
        out[n] = []
        for sector in SELECTION[n]:
            meta, raw_t, raw = cube(case, sector, *options[sector], False)
            original = bin6(raw_t, raw)
            t = original[:, 0]
            peaks = [
                pilot.utc_jd(
                    date.fromisoformat(g["r"]["date"] if g["r"] else g["first"])
                )
                + 0.5
                for g in selection
            ]
            train = np.ones(len(t), dtype=bool)
            for p in peaks:
                train &= abs(t - p) > 1.5
            arr = corrected(original, train)
            met = previous_tess.compare(meta, arr, pre)
            bands = met["bands"]
            passed = any(
                v.get("matched_nights", 0) >= 6
                and v.get("spearman", {}).get("center", 0) >= 0.65
                and v.get("peak_distance_pixels", 999) <= 1
                for v in bands.values()
            )
            per = {
                "sector": sector,
                "n_six_hour_bins": len(t),
                "training_bins": int(train.sum()),
                "source_localized": passed,
                "corrected_ztf_matched": {
                    b: {k: v for k, v in x.items() if k != "matches"}
                    for b, x in bands.items()
                },
                "events": [],
            }
            pixs = {
                "raw": original[:, 1:].reshape(-1, 9, 9),
                "corrected": arr[:, 1:].reshape(-1, 9, 9),
            }
            for g, p in zip(selection, peaks):
                if p < t[0] - 1.5 or p > t[-1] + 1.5:
                    continue
                bright = abs(t - p) <= 1.5
                faint = np.zeros(len(t), dtype=bool)
                brackets = []
                for k in ("previous_faint", "next_faint"):
                    if g.get(k):
                        dd = pilot.utc_jd(date.fromisoformat(g[k])) + 0.5
                        match = abs(t - dd) <= 1.5
                        if abs(dd - p) <= 10 and match.any():
                            brackets.append(g[k])
                            faint |= match
                overlap = bright & faint
                bright &= ~overlap
                faint &= ~overlap
                per["events"].append(
                    {
                        "group": g["first"],
                        "ztf_red_peak": g["r"]["date"] if g["r"] else None,
                        "brackets": brackets,
                        "bright_bins": int(bright.sum()),
                        "faint_bins": int(faint.sum()),
                        **{
                            lab: contrast(px, bright, faint, meta["target_xy"])
                            for lab, px in pixs.items()
                        },
                    }
                )
            a = pixs["corrected"][:, 3:6, 3:6].sum(axis=(1, 2))
            if n == "ZTF17aaaoxxi":
                cand = plateaus(t, a)
                per["longest_candidate_days"] = max(
                    (v["sampled_days"] for v in cand), default=None
                )
                per["eligible_plateaus"] = [v for v in cand if v["accepted"] and passed]
            out[n].append(per)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=HERE / "tess-events.json")
    pilot.json_write(p.parse_args().out, run())
