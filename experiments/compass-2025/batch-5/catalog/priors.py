"""Resolve the five-survey supernova report and published dwarf-nova classifications."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import html
import io
import json
import re
import tempfile
import zipfile
from pathlib import Path
from urllib.parse import urlencode
from xml.etree import ElementTree
import importlib.util

from check import HERE, archive

STARS = {
    "MGAB-V960": (307.1167854, 52.1020668, "1499309", "mgab-v951-v1000"),
    "MGAB-V306": (102.2581637, 8.995694, "1498454", "mgab-v301-v350"),
}
spec = importlib.util.spec_from_file_location(
    "batch_four_priors", HERE.parents[1] / "batch-4/catalog/priors.py"
)
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)
previous.HERE = HERE
previous.archive = archive
previous.POSITIONS = {n: v[:2] for n, v in STARS.items()}


def row_archive(key, url, raw, selected, extra=None):
    (HERE / (key + ".gz")).write_bytes(gzip.compress(selected, mtime=0))
    metadata = {
        "url": url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256_full_response": previous.digest(raw),
        "sha256_uncompressed": previous.digest(selected),
    }
    archive.pilot.json_write(HERE / (key + ".json"), metadata | (extra or {}))


def tns(fetch):
    url = previous.TNS
    if fetch and not (HERE / "tns.csv.json").exists():
        raw = archive.pilot.curl(url, ["-A", previous.GUEST])
        z = zipfile.ZipFile(io.BytesIO(raw))
        lines = z.read("tns_public_objects.csv").splitlines(keepends=True)
        (match,) = [x for x in lines[2:] if b'"2025xza"' in x and b"ZTF25absffwy" in x]
        row_archive("tns.csv", url, raw, lines[1] + match, {"staged_utc": lines[0].decode().strip()})
    (r,) = csv.DictReader(io.StringIO(archive.source("tns.csv", url, False).decode()))
    if r["name"] != "2025xza" or "ZTF25absffwy" not in r["internal_names"]:
        raise ValueError("Wrong TNS object")
    return r


def discovery(star, fetch):
    page = "https://sites.google.com/view/mgab-astronomy/" + STARS[star][3]
    key = star + "-discovery"
    if fetch and not (HERE / (key + ".json")).exists():
        url = (
            page
            + "?"
            + urlencode({"view": "published", "date": datetime.now(timezone.utc).strftime("%Y%m%d")})
        )
        raw = archive.pilot.curl(url, ["--location"])
        sections = re.findall(r"<section\b[^>]*>.*?</section>", raw.decode(), re.S)
        def heading(section):
            match = re.search(r'<h2\b[^>]*>(.*?)</h2>', section, re.S)
            if match:
                return "".join(html.unescape(re.sub(r"<[^>]+>", "", match.group(1))).split())
        (section,) = [s for s in sections if heading(s) == star]
        text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", section)).split())
        images = re.findall(r'<img\b[^>]*\bsrc="([^"]+)"', section)
        chart_hashes = []
        for index, img in enumerate(images):
            graph = archive.pilot.curl(html.unescape(img), ["--location", "-e", page])
            if not graph.startswith(b"\x89PNG"):
                raise ValueError("Unexpected discovery chart image")
            chart_hashes.append(previous.digest(graph))
            # Transient local copy for visual inspection of the discoverer's axes and legend.
            (Path(tempfile.gettempdir()) / f"fink-b5-{star}-{index}.png").write_bytes(graph)
        row_archive(
            key, page, raw, (text + "\n").encode(), {"chart_sha256": chart_hashes, "request_url": url}
        )
    return archive.source(key, page, False).decode().strip()


def public_row(key, url, fetch, selector):
    if fetch and not (HERE / (key + ".json")).exists():
        raw = archive.pilot.curl(url, ["--location"])
        selected = selector(raw.decode())
        row_archive(key, url, raw, (selected + "\n").encode())
    return archive.source(key, url, False).decode().strip()


def cats_row(page):
    (row,) = [r for r in re.findall(r"<tr\b[^>]*>.*?</tr>", page, re.S) if "ZTF25absffwy" in r]
    return " | ".join(" ".join(html.unescape(re.sub(r"<[^>]+>", "", el)).split())
                      for el in re.findall(r"<td\b[^>]*>(.*?)</td>", row, re.S))


def observer_dates(key, name, jd_start, jd_end, fetch):
    params = {"view": "api.delim", "ident": name, "fromjd": jd_start, "tojd": jd_end, "delimiter": "@@@"}
    url = "https://vsx.aavso.org/index.php?" + urlencode(params)
    lines = archive.source(key, url, fetch).decode("utf-8-sig").splitlines()
    fields = lines[0].split("@@@")
    if not {"obsID", "starName", "JD", "mag"} <= set(fields):
        raise ValueError("Unexpected observer table")
    return [dict(zip(fields, line.split("@@@"), strict=True)) for line in lines[1:]]


def run(fetch):
    result = {"tns": tns(fetch), "cats150": public_row(
        "cats150", "https://cats150.caltech.edu/sources", fetch, cats_row
    ), "dwarf_novae": {}}
    for n, (_, _, oid, _) in STARS.items():
        url = "https://vsx.aavso.org/index.php?" + urlencode({"view": "api.object", "ident": n})
        row = {el.tag: el.text for el in ElementTree.fromstring(archive.source(n + "-vsx.xml", url, fetch))}
        if row["OID"] != oid:
            raise ValueError("VSX identity mismatch")
        tables = [
            "J/AJ/159/198/table1",
            "J/AJ/159/198/table2",
            "J/AJ/159/198/table3",
            "J/AJ/162/94/table1",
            "J/AJ/162/94/table2",
            "J/MNRAS/456/4441/catalog",
            "J/AJ/153/204/table5",
        ]
        candidates = previous.literature(n, tables, fetch)
        obs = observer_dates(n + "-aid.txt", n, 2458119.5, 2461041.5, fetch)
        result["dwarf_novae"][n] = {
            "vsx": row,
            "published_rows": candidates,
            "discovery": discovery(n, fetch),
            "aavso_observations": len(obs),
            "aavso_dates": sorted({str(archive.pilot.day_of(float(v["JD"]))) for v in obs}),
        }
    result["dwarf_novae"]["MGAB-V306"]["aavso_2010_2014_count"] = len(
        observer_dates("MGAB-V306-early-aid.txt", "MGAB-V306", 2455197.5, 2457023.5, fetch)
    )
    def vsnet_row(page):
        (line,) = [s for s in page.splitlines() if s.lstrip().startswith("MGAB-V306")]
        return " ".join(html.unescape(line).split())
    result["dwarf_novae"]["MGAB-V306"]["vsnet_january_2025"] = public_row(
        "vsnet-january", "http://ooruri.kusastro.kyoto-u.ac.jp/mailarchive/vsnet-cvcirc/8794",
        fetch, vsnet_row
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
