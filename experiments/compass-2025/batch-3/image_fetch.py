"""Extract a small ZTF science cutout from its public original when IBE cutouts fail."""

from datetime import datetime, timezone
import gzip
import hashlib
import io
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from compass_aperture import cutout_url

HEADER_RANGE = 65536
MAX_BYTES = 40_000_000


def request_range(url, start, end, budget):
    """Read one bounded response, checking whether the server honored the range."""
    if budget < end - start + 1:
        raise ValueError("Range exceeds remaining download budget")
    with TemporaryDirectory() as directory:
        headers = Path(directory) / "headers"
        body = Path(directory) / "body"
        cmd = [
            "curl",
            "--fail",
            "--silent",
            "--show-error",
            "--location",
            "--max-time",
            "120",
            "--max-filesize",
            str(budget),
            "--range",
            f"{start}-{end}",
            "--dump-header",
            str(headers),
            "--output",
            str(body),
            "--write-out",
            "%{http_code}",
            url,
        ]
        status = int(subprocess.check_output(cmd, text=True))
        if body.stat().st_size > budget:
            raise ValueError("Download exceeded byte budget")
        raw = body.read_bytes()
        hdr = headers.read_text()
    if status == 206:
        ranges = re.findall(
            r"^content-range:\s*bytes (\d+)-(\d+)/(\d+)\s*$", hdr, re.I | re.M
        )
        if not ranges:
            raise ValueError("Missing Content-Range")
        first, last, total = map(int, ranges[-1])
        if (first, last) != (start, min(end, total - 1)) or len(
            raw
        ) != last - first + 1:
            raise ValueError("Incomplete or mismatched range response")
    elif status == 200:
        total = len(raw)
    else:
        raise ValueError(f"Unexpected original-image response {status}")
    meta = {
        "requested_bytes": [start, end],
        "http_status": status,
        "received_bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }
    return raw, total, meta


def read_header(data):
    for i in range(0, len(data) - 79, 80):
        if data[i : i + 80] == b"END".ljust(80):
            offset = (i + 80 + 2879) // 2880 * 2880
            return (
                fits.Header.fromstring(data[: i + 80].decode("ascii"), sep=""),
                offset,
            )
    raise ValueError("No FITS primary header END within header range")


def from_original(
    row, position, size=165, max_bytes=MAX_BYTES, *, read_range=request_range
):
    """Return (losslessly cropped FITS bytes, source receipt). Bound total response bytes."""
    if not (0 < max_bytes <= MAX_BYTES and np.isfinite(size) and 0 < size <= 600):
        raise ValueError("Invalid transfer ceiling or cutout size")
    url = cutout_url(row, *position, size=size).split("?", 1)[0]
    first, total, meta = read_range(url, 0, HEADER_RANGE - 1, max_bytes)
    requests = [meta]
    h, offset = read_header(first)
    if (
        not h.get("SIMPLE")
        or h.get("NAXIS") != 2
        or h.get("BITPIX") != -32
        or h.get("BSCALE", 1) != 1
        or h.get("BZERO", 0) != 0
        or not (
            h.get("EXPID") == int(row["expid"])
            and h.get("CCDID") == int(row["ccdid"])
            and h.get("QID") == int(row["qid"])
            and h.get("FILTER") == "ZTF_" + row["filtercode"][1:]
        )
    ):
        raise ValueError(
            "FITS type, exposure, detector or filter does not match inventory"
        )
    width, height = h["NAXIS1"], h["NAXIS2"]
    if min(width, height) <= 0 or total < offset + 4 * width * height:
        raise ValueError("Truncated original image")
    wcs = WCS(h)
    center = wcs.all_world2pix([position], 0)[0]
    if not (
        np.isfinite(center).all()
        and np.all(center >= -0.5)
        and np.all(center < np.array([width, height]) - 0.5)
    ):
        raise ValueError("Cutout center lies outside the original quadrant")
    half = 0.5 * size / (3600 * proj_plane_pixel_scales(wcs))
    start = np.maximum(np.floor(center - half + 0.5).astype(int), 0)
    stop = np.minimum(np.floor(center + half + 0.5).astype(int) + 1, [width, height])
    x0, y0 = map(int, start)
    x1, y1 = map(int, stop)
    p0 = offset + 4 * (y0 * width + x0)
    p1 = offset + 4 * ((y1 - 1) * width + x1) - 1
    if meta["http_status"] == 200:
        data, pos = first, p0
    else:
        data, second_total, meta = read_range(url, p0, p1, max_bytes - len(first))
        requests.append(meta)
        if second_total != total:
            raise ValueError("Original changed during transfer")
        pos = p0 if meta["http_status"] == 200 else 0
    pixels = np.ndarray(
        (y1 - y0, x1 - x0), dtype=">f4", buffer=data, offset=pos, strides=(4 * width, 4)
    ).copy()
    h["CRPIX1"] -= x0
    h["CRPIX2"] -= y0
    for key in ["DATASUM", "CHECKSUM"]:
        h.remove(key, ignore_missing=True)
    f = io.BytesIO()
    fits.PrimaryHDU(pixels, header=h).writeto(f)
    raw = f.getvalue()
    return raw, {
        "method": "local_crop",
        "original_url": url,
        "original_bytes": total,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "transfers": requests,
        "slice_xy_zero_based_half_open": [x0, y0, x1, y1],
        "source_header_bytes": offset,
        "size_arcsec": size,
        "sha256_cutout": hashlib.sha256(raw).hexdigest(),
    }


def from_alternate_cutout(row, position, size=165):
    """Try the documented literal-coordinate URL with default and explicit gzip."""
    if not (0 < size <= 600 and np.isfinite(position).all()):
        raise ValueError("Invalid center or cutout size")
    original_url = cutout_url(row, *position, size=size).split("?", 1)[0]
    base = f"{original_url}?center={position[0]},{position[1]}&size={size}arcsec"
    attempts = []
    bound = 2_000_000
    with TemporaryDirectory() as directory:
        path = Path(directory) / "cutout"
        for url in (base, base + "&gzip=true"):
            command = [
                "curl",
                "--silent",
                "--show-error",
                "--location",
                "--fail",
                "--max-time",
                "40",
                "--max-filesize",
                str(bound),
                "--output",
                str(path),
                "--write-out",
                "%{http_code}",
                url,
            ]
            response = subprocess.run(command, stdout=subprocess.PIPE)
            attempt = {
                "url": url,
                "http_status": int(response.stdout or b"0"),
                "curl_status": response.returncode,
            }
            attempts.append(attempt)
            if response.returncode in (22, 28):
                continue
            response.check_returncode()
            if path.stat().st_size > bound or attempt["http_status"] != 200:
                raise ValueError("Invalid alternate cutout response")
            payload = path.read_bytes()
            if payload.startswith(b"\x1f\x8b"):
                with gzip.GzipFile(fileobj=io.BytesIO(payload)) as gz:
                    raw = gz.read(bound + 1)
            else:
                raw = payload
            if len(raw) > bound or not raw.startswith(b"SIMPLE  ="):
                raise ValueError("Malformed or excessive alternate FITS cutout")
            with fits.open(io.BytesIO(raw)) as f:
                h = f[0].header
                if (
                    h["EXPID"] != int(row["expid"])
                    or h["CCDID"] != int(row["ccdid"])
                    or h["QID"] != int(row["qid"])
                    or h["FILTER"] != "ZTF_" + row["filtercode"][1:]
                ):
                    raise ValueError("Alternate cutout identity mismatch")
                wcs = WCS(h)
                xy = wcs.all_world2pix([position], 0)[0]
                dims = np.array(f[0].data.shape[::-1])
                if (
                    not np.isfinite(xy).all()
                    or not np.all(xy >= -0.5)
                    or not np.all(xy < dims - 0.5)
                    or np.any(
                        dims > np.ceil(size / (3600 * proj_plane_pixel_scales(wcs))) + 3
                    )
                ):
                    raise ValueError(
                        "Alternate cutout size or target position mismatch"
                    )
            return raw, {
                "method": "alternate_ibe_cutout",
                "original_url": original_url,
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(
                    timespec="seconds"
                ),
                "attempts": attempts,
                "received_bytes": len(payload),
                "size_arcsec": size,
                "sha256_response": hashlib.sha256(payload).hexdigest(),
                "sha256_cutout": hashlib.sha256(raw).hexdigest(),
            }
    response.check_returncode()
