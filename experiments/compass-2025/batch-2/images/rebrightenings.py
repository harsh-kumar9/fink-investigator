"""Check the positions in two colors at the start and during the late brightening."""

import argparse
import json
import math

import replay
import compass_direct

HERE = replay.HERE
NAME = "ZTF25aacfjde"
BANDS = {"g": 1, "r": 2}
SAMPLE = {
    "2025-01-19": ["g", "r"],
    "2025-02-10": ["g"],
    "2025-02-20": ["g"],
    "2025-04-16": ["g", "r"],
}


def cutout_folder(day, band):
    if day == "2025-01-19" and band == "g":
        return HERE / NAME / "epochs"
    return HERE / NAME / (day if band == "g" else day + "-r")


def run(fetch=False):
    im = replay.images
    im.HERE = HERE
    anchors, histories, _ = im.compass_archive.load(replay.YEAR / "snapshot")
    a = next(a for a in anchors if im.pilot.oid(a) == NAME)
    history = histories[NAME]
    stars = im.gaia(a, False)
    candidates = sorted(stars, key=lambda s: im.sep(a, s, True))
    if not candidates or im.sep(a, candidates[0], True) > 5:
        raise ValueError("Source assignment changed")
    star = candidates[0]
    case = {"ra": a["i:ra"], "dec": a["i:dec"]}
    _, photometry, _ = compass_direct.direct_data(case, replay.BATCH / "direct" / NAME)
    result = []
    for day, bands in SAMPLE.items():
        for band in bands:
            rows = [r for r in history if str(im.pilot.day_of(r["i:jd"])) == day
                    and r["i:fid"] == BANDS[band] and im.quality(r) and im.positive(r) is True]
            if len(rows) != 1:
                raise ValueError(f"Unexpected exposure count for {day} {band}")
            row = rows[0]
            cand = str(row["i:candid"])
            folder = cutout_folder(day, band)
            if folder.name == "epochs":
                if cand not in im.stamps(history, a, folder, False):
                    raise ValueError("Expected g image in original packet")
            else:
                im.stamps(rows, row, folder, fetch)
            diff_offset = im.relative_position(row, star, True)
            # ZTF DR24 exposure IDs match the observing timestamp prefix in the alert candid.
            matches = [d for d in photometry if cand.startswith(d["expid"])
                       and d["filtercode"] == "z" + band and int(d["catflags"]) == 0
                       and abs((float(d["mjd"]) + 2400000.5) - row["i:jd"]) < 60 / 86400]
            if len(matches) > 1:
                raise ValueError("Multiple direct positions for exposure")
            record = {
                "date": day, "band": band, "candid": cand, "gaia_id": star["source_id"],
                "difference_mag": row["i:magpsf"], "difference_magerr": row["i:sigmapsf"],
                "difference_xy_arcsec": [round(v, 3) for v in diff_offset],
                "difference_separation_arcsec": round(math.hypot(*diff_offset), 3),
                "next_gaia_separation_arcsec": round(min(im.sep(row,s,True) for s in candidates[1:]), 3),
                "direct": None,
            }
            if matches:
                d = matches[0]
                direct_pos = {"i:jd": float(d["mjd"]) + 2400000.5,
                              "i:ra": float(d["ra"]), "i:dec": float(d["dec"])}
                offset = im.relative_position(direct_pos, star, True)
                record["direct"] = {
                    "expid": d["expid"], "magnitude": float(d["mag"]), "magerr": float(d["magerr"]),
                    "xy_arcsec": [round(v, 3) for v in offset],
                    "separation_arcsec": round(math.hypot(*offset), 3),
                    "difference_to_direct_arcsec": round(math.dist(diff_offset, offset), 3),
                }
            result.append(record)
    return result


def plot(results):
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(len(results), 3, figsize=(4.5, 1.4 * len(results) + 0.45),
                             layout="constrained", squeeze=False)
    for i, r in enumerate(results):
        folder = cutout_folder(r["date"], r["band"])
        for j, kind in enumerate(("Science", "Template", "Difference")):
            ax = axes[i,j]
            image = plt.imread(folder / (r["candid"] + "-" + kind.lower() + ".png"))
            ax.imshow(image, cmap="gray", interpolation="nearest", vmin=0, vmax=1)
            ax.set(xticks=[], yticks=[])
            if i == 0:
                ax.set_title(kind, fontsize=11)
            if j == 0:
                ax.set_ylabel(r["date"] + "  " + r["band"], fontsize=10)
    fig.suptitle(NAME, fontsize=14)
    fig.savefig(HERE / NAME / "followup.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()
    results = run(args.fetch)
    file = HERE / "rebrightenings.json"
    if args.fetch:
        replay.images.pilot.json_write(file, results)
    elif json.loads(file.read_text()) != results:
        raise ValueError("Source positions changed")
    if args.plot:
        plot(results)
    print("Verified", len(results), "exposures")
