"""Check source-specific publications, the TNS discovery record and public spectra."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import io
import json
import re
import zipfile
from urllib.parse import urlencode
from xml.etree import ElementTree

from check import HERE, catalog as previous

archive = previous.catalog
TNS = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
GUEST = 'tns_marker{"tns_id":0,"type":"guest","name":"Codex astronomy archive query"}'
STARS = {
    "MGAB-V488": (44.4238904, 48.862231, "1498636"),
    "ZTF17aabutiy": (329.5663886, 51.8389644, "1542955"),
}
POSITIONS = {name: value[:2] for name, value in STARS.items()} | {"AT2025iac": (250.0278614, 29.763246)}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def tns(fetch):
    # Retain the exact header and named row of the large, daily revised public registry.
    path, receipt = HERE / "tns.csv.gz", HERE / "tns.csv.json"
    if fetch and not receipt.exists():
        raw = archive.pilot.curl(TNS, ["-A", GUEST])
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            data = z.read("tns_public_objects.csv").splitlines(keepends=True)
        rows = [line for line in data[2:] if b'"2025iac"' in line and b"ZTF25aanlbls" in line]
        (target,) = rows
        subset = data[1] + target
        path.write_bytes(gzip.compress(subset, mtime=0))
        archive.pilot.json_write(
            receipt,
            {
                "url": TNS,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "staged_utc": data[0].decode().strip(),
                "sha256_full_zip": digest(raw),
                "sha256_uncompressed": digest(subset),
            },
        )
    (row,) = list(csv.DictReader(io.StringIO(archive.source("tns.csv", TNS, False).decode())))
    if row["name"] != "2025iac" or row["internal_names"] != "ZTF25aanlbls":
        raise ValueError("Unexpected registry source")
    return row


def literature(star, tables, fetch, label="literature"):
    ra, dec = POSITIONS[star]
    params = {"-source": ",".join(tables), "-c": f"{ra} {dec:+}", "-c.rs": 5, "-out.all": 1, "-out.max": 20}
    raw = archive.source(f"{star}-{label}.xml", archive.vsx.SERVICE + "?" + urlencode(params), fetch)
    result = {}
    for tab in ElementTree.fromstring(raw).findall(".//{*}TABLE"):
        fields = [f.get("name") for f in tab.findall("{*}FIELD")]
        rows = [
            dict(zip(fields, [v.text.strip() if v.text else None for v in tr], strict=True))
            for tr in tab.findall(".//{*}TR")
        ]
        if len(rows) >= 20 or any(float(row["_r"]) > 5.01 for row in rows):
            raise ValueError("Unexpected published cone rows")
        result[tab.get("name")] = rows
    if set(result) != set(tables):
        raise ValueError("Incomplete literature query")
    return result


def spectra(fetch):
    sdss = {}
    for name, (ra, dec) in POSITIONS.items():
        query = (
            "SELECT TOP 20 s.specObjID, s.ra, s.dec, s.class,s.z,s.zWarning FROM "
            f"dbo.fGetNearbySpecObjEq({ra},{dec},0.166667) AS n JOIN SpecObj AS s ON s.specObjID=n.specObjID"
        )
        url = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?" + urlencode(
            {"cmd": query, "format": "csv"}
        )
        raw = archive.source(name + "-sdss.csv", url, fetch)
        rows = csv.DictReader(line for line in raw.decode().splitlines() if not line.startswith("#"))
        if rows.fieldnames != ["specObjID", "ra", "dec", "class", "z", "zWarning"]:
            raise ValueError("Unexpected SDSS primary spectrum schema")
        sdss[name] = list(rows)
        if len(sdss[name]) >= 20:
            raise ValueError("Truncated SDSS query")
    query = " union all ".join(
        f"select '{name}' as src,targetid,mean_fiber_ra,mean_fiber_dec,spectype,z,zwarn from desi_dr1.zpix where q3c_radial_query(mean_fiber_ra,mean_fiber_dec,{ra},{dec},0.002777778)"
        for name, (ra, dec) in POSITIONS.items()
    )
    url = "https://datalab.noirlab.edu/query/query?" + urlencode(
        {"sql": query, "ofmt": "csv", "async": "False"}
    )
    receipt = HERE / "desi.csv.json"
    if fetch and not receipt.exists():
        raw = archive.pilot.curl(url, ["-H", "X-DL-AuthToken: anonymous.0.0.anon_access"])
        (HERE / "desi.csv.gz").write_bytes(gzip.compress(raw, mtime=0))
        archive.pilot.json_write(
            receipt,
            {
                "url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "sha256_uncompressed": digest(raw),
            },
        )
    d = csv.DictReader(io.StringIO(archive.source("desi.csv", url, False).decode()))
    if d.fieldnames != ["src", "targetid", "mean_fiber_ra", "mean_fiber_dec", "spectype", "z", "zwarn"]:
        raise ValueError("Unexpected DESI spectrum schema")
    return {"sdss_dr17_primary": sdss, "desi_dr1_zpix": list(d)}


def historical_notices(fetch):
    pages = {
        "vsnet": "http://ooruri.kusastro.kyoto-u.ac.jp/mailarchive/vsnet-cvcirc/7258",
        "murawski": "https://sites.google.com/view/mgab-astronomy/mgab-v451-v500",
    }
    notices = {}
    for name, url in pages.items():
        receipt = HERE / (name + ".json")
        if fetch and not receipt.exists():
            raw = archive.pilot.curl(url, ["--location"])
            page = raw.decode("utf-8", errors="replace")
            if name == "vsnet":
                labels = ("MGAB-V488 ", "ZTF17aabutiy ", "ZAD ZTF Alert DC mag")
                lines = [line.strip() for line in page.splitlines() if line.startswith(labels)]
                if len(lines) != 3:
                    raise ValueError("Missing star or source key in VSNET circular")
            else:
                sections = re.findall(r"<section\b[^>]*>.*?</section>", page, re.S)
                (section,) = [s for s in sections if ">MGAB-V488</div>" in s]
                prose = html.unescape(re.sub(r"<[^>]+>", " ", section))
                lines = [" ".join(prose.split())]
                if not all(s in lines[0] for s in ["MGAB-V488", "44.423839", "UG", "18.4 - 21.4 r"]):
                    raise ValueError("Unexpected discovery section")
            subset = ("\n".join(lines) + "\n").encode()
            (HERE / (name + ".gz")).write_bytes(gzip.compress(subset, mtime=0))
            archive.pilot.json_write(
                receipt,
                {
                    "url": url,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "sha256_full_response": digest(raw),
                    "sha256_uncompressed": digest(subset),
                },
            )
        notices[name] = archive.source(name, url, False).decode().splitlines()
    return notices


def discovery_chart(verify=False):
    """The discoverer's public page provides a short-lived URL for the chart image."""
    receipt = HERE / "discovery-chart.json"
    if verify or not receipt.exists():
        page_url = "https://sites.google.com/view/mgab-astronomy/mgab-v451-v500"
        date = datetime.now(timezone.utc).strftime("%Y%m%d")
        url = page_url + "?" + urlencode({"view": "published", "date": date})
        page = archive.pilot.curl(url, ["--location"])
        sections = re.findall(r"<section\b[^>]*>.*?</section>", page.decode(), re.S)
        (section,) = [s for s in sections if ">MGAB-V488</div>" in s]
        (image_url,) = re.findall(r'<img\b[^>]*\bsrc="([^"]+)"', section)
        raw = archive.pilot.curl(html.unescape(image_url), ["--location", "-e", page_url])
        if not raw.startswith(b"\x89PNG\r\n\x1a\n"):
            raise ValueError("The discovery chart is unavailable")
        if receipt.exists():
            if digest(raw) != json.loads(receipt.read_text())["sha256_chart"]:
                raise ValueError("The discoverer has changed the published chart")
        else:
            archive.pilot.json_write(
                receipt,
                {
                    "page_url": page_url,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "sha256_page": digest(page),
                    "sha256_chart": digest(raw),
                    "bytes_chart": len(raw),
                    "visually_read_chart_legend": ["g-band (ZTF)", "r-band (ZTF)"],
                    "visually_read_axis_HJD_minus_2400000": [58150, 59050],
                },
            )
    return json.loads(receipt.read_text())


def run(fetch):
    result = {"tns": tns(fetch), "primary_rows": {}, "aavso": {}, "vsx": {}}
    for name, (_, _, oid) in STARS.items():
        url = "https://vsx.aavso.org/index.php?" + urlencode({"view": "api.object", "ident": name})
        row = {v.tag: v.text for v in ElementTree.fromstring(archive.source(name + "-vsx.xml", url, fetch))}
        if row["Name"] != name or row["OID"] != oid:
            raise ValueError("Wrong VSX source")
        result["vsx"][name] = row
        aavso = name + "-aid.txt"
        url = "https://vsx.aavso.org/index.php?" + urlencode(
            {"view": "api.delim", "ident": name, "fromjd": 2458119.5, "tojd": 2461041.5, "delimiter": "@@@"}
        )
        lines = archive.source(aavso, url, fetch).decode("utf-8-sig").splitlines()
        header = lines[0].split("@@@")
        if not {"starName", "obsID", "JD", "band", "mag"} <= set(header):
            raise ValueError("Unexpected AAVSO observer archive")
        rows = [dict(zip(header, l.split("@@@"), strict=True)) for l in lines[1:]]
        result["aavso"][name] = {
            "rows": len(rows),
            "names": sorted({r["starName"] for r in rows}),
            "dates": sorted({str(archive.pilot.day_of(float(r["JD"]))) for r in rows}),
        }
    result["primary_rows"]["MGAB-V488"] = literature("MGAB-V488", ["J/AJ/153/204/table5"], fetch)
    result["primary_rows"]["ZTF17aabutiy"] = literature(
        "ZTF17aabutiy", ["J/AJ/159/198/table1", "J/AJ/159/198/table2", "J/AJ/159/198/table3"], fetch
    )
    result["historical_notices"] = historical_notices(fetch)
    if fetch or (HERE / "discovery-chart.json").exists():
        result["discovery_chart"] = discovery_chart()
    for name in STARS:
        result["primary_rows"][name].update(
            literature(name, ["J/MNRAS/456/4441/catalog"], fetch, "crts-outburst")
        )
    result["spectra_within_10arcsec"] = spectra(fetch)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument(
        "--verify-chart", action="store_true", help="Check the live chart against the observed checksum"
    )
    args = p.parse_args()
    if args.verify_chart:
        discovery_chart(verify=True)
    print(json.dumps(run(args.fetch), indent=2))
