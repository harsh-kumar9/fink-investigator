"""Retrieve DR24 standards farther from the disk on the same quadrant."""

from pathlib import Path
from urllib.parse import urlencode
import math
import numpy as np
from standards import HERE, POS, late, pilot, capture


def catalog(band, fetch):
    fid = {"r": 2, "i": 3}[band]
    cols = "oid,ra,dec,field,ccdid,qid,fid,nobs,meanmag"
    sql = f"SELECT TOP 6000 {cols} FROM ztf_objects_dr24 WHERE field=401 AND ccdid=9 AND qid=3 AND fid={fid} AND CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{POS[0]},{POS[1]},0.08333334))=1"
    url = (
        pilot.IRSA
        + "?"
        + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=sql))
    )
    rows = capture(
        HERE / "standards-sources", "wide-" + band, url, set(cols.split(",")), fetch
    )
    if len(rows) >= 6000 or len({r["oid"] for r in rows}) != len(rows):
        raise ValueError("Truncated wide cone")
    return rows


def standards(rows, band):
    meta = dict(field="401", ccdid="9", qid="3", filtercode="z" + band)
    field, xy, inds = late.phot.star_list(rows, meta, *POS, radius=280)
    selected = []
    for i in inds:
        distance = math.hypot(*pilot.offsets(*POS, *xy[i]))
        if distance > 40:
            selected.append((distance, field[i]))
    return [
        r for d, r in sorted(selected, key=lambda item: (item[0], item[1]["oid"]))[:8]
    ]


if __name__ == "__main__":
    for band in "ri":
        rows = catalog(band, True)
        v = standards(rows, band)
        print(band, len(rows), [(r["oid"], r["meanmag"]) for r in v])


def tile_counts(row, refs, central_folder, fetch):
    import gzip
    import hashlib
    import io
    import json
    import subprocess
    from datetime import datetime, timezone
    from astropy.io import fits
    from astropy.wcs import WCS

    directory = HERE / "standards-sources" / "tiles"
    directory.mkdir(exist_ok=True)
    pid = row["pid"]
    path, receipt = directory / (pid + ".fits.gz"), directory / (pid + ".json")
    url = late.phot.cutout_url(row, *POS, size=600)
    if fetch and not receipt.exists():
        errors = []
        for size in (600, 480, 510):
            requested = late.phot.cutout_url(row, *POS, size=size)
            try:
                raw = pilot.curl(requested, [])
                method = {"method": "ibe", "url": requested}
                break
            except subprocess.CalledProcessError as e:
                errors.append([requested, e.returncode])
            try:
                raw, method = late.image_fetch.from_alternate_cutout(
                    row, POS, size=size
                )
                break
            except subprocess.CalledProcessError as e:
                errors.append(["alternate-" + str(size), e.returncode])
        else:
            try:
                raw, method = late.image_fetch.from_original(row, POS, size=480)
            except subprocess.CalledProcessError as e:
                errors.append(["original", e.returncode])
                pilot.json_write(
                    receipt,
                    dict(
                        url=url,
                        pid=pid,
                        expid=row["expid"],
                        oids=[r[0] for r in refs],
                        method="unavailable",
                        attempts=errors,
                        archived_utc=datetime.now(timezone.utc).isoformat(
                            timespec="seconds"
                        ),
                    ),
                )
                return {}
        with fits.open(io.BytesIO(raw)) as hdus, fits.open(
            io.BytesIO(late.image(row, False, POS, central_folder))
        ) as cen:
            h, data = hdus[0].header, hdus[0].data
            for k in ("EXPID", "FIELDID", "QID", "CCD_ID", "FILTER"):
                if h[k] != cen[0].header[k]:
                    raise ValueError("Cutouts are from different exposures")
            w = WCS(h)
            x, y = map(round, w.all_world2pix([POS], 0)[0])
            coord = w.all_pix2world([[x, y]], 0)
            cx, cy = map(round, WCS(cen[0].header).all_world2pix(coord, 0)[0])
            core = data[y - 5 : y + 6, x - 5 : x + 6]
            if not np.array_equal(
                core, cen[0].data[cy - 5 : cy + 6, cx - 5 : cx + 6], equal_nan=True
            ):
                raise ValueError("Same-exposure central pixels disagree")
            output = [fits.PrimaryHDU()]
            for oid, pos, _ in refs:
                fx, fy = w.all_world2pix([pos], 0)[0]
                ix, iy = round(fx), round(fy)
                tile = data[iy - 13 : iy + 14, ix - 13 : ix + 14]
                if tile.shape != (27, 27):
                    raise ValueError("Partial stellar tile")
                ext = fits.ImageHDU(tile.copy(), name=oid)
                ext.header["XPOS"] = fx - ix + 13
                ext.header["YPOS"] = fy - iy + 13
                ext.header["RA"] = pos[0]
                ext.header["DEC"] = pos[1]
                output.append(ext)
            payload = io.BytesIO()
            fits.HDUList(output).writeto(payload)
        data = payload.getvalue()
        path.write_bytes(gzip.compress(data, mtime=0))
        pilot.json_write(
            receipt,
            dict(
                url=url,
                request=method,
                pid=pid,
                expid=row["expid"],
                oids=[r[0] for r in refs],
                cutout_sha256=hashlib.sha256(raw).hexdigest(),
                tiles_sha256=hashlib.sha256(data).hexdigest(),
                archived_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ),
        )
    meta = json.loads(receipt.read_text())
    if meta.get("method") == "unavailable":
        if (
            meta["url"] != url
            or meta["pid"] != pid
            or meta["oids"] != [r[0] for r in refs]
            or path.exists()
        ):
            raise ValueError("Invalid unavailable receipt")
        return {}
    data = gzip.decompress(path.read_bytes())
    if (
        meta["url"] != url
        or meta["pid"] != pid
        or meta["expid"] != row["expid"]
        or meta["oids"] != [r[0] for r in refs]
        or meta["tiles_sha256"] != hashlib.sha256(data).hexdigest()
    ):
        raise ValueError("Stellar tiles changed")
    with fits.open(io.BytesIO(data)) as hdus:
        result = {}
        for (oid, pos, _), ext in zip(refs, hdus[1:], strict=True):
            h = ext.header
            if (
                h["EXTNAME"] != oid
                or abs(h["RA"] - pos[0]) + abs(h["DEC"] - pos[1]) > 1.0e-9
            ):
                raise ValueError("Wrong calibration star")
            result[oid] = late.phot.aperture(ext.data, h["XPOS"], h["YPOS"])
        return result
