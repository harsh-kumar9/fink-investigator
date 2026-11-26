"""Check earlier variable classifications, Gaia astrometry and independent observer coverage."""

import argparse
from datetime import datetime, timezone
import email
from email import policy
import gzip
import hashlib
import json
import math
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

from check import HERE, archive


def october_notices(fetch):
    addr = "vsnet-outburst@ooruri.kusastro.kyoto-u.ac.jp"
    url = (f"http://ooruri.kusastro.kyoto-u.ac.jp/mailman3/hyperkitty/list/{addr}/export/"
           f"{addr}-2024-10.mbox.gz?start=2024-10-10&end=2024-10-14")
    receipt = HERE / "october-notices.json"
    if fetch and not receipt.exists():
        raw = archive.pilot.curl(url, [])
        blocks = re.split(br"(?m)^From .+?\n", gzip.decompress(raw))[1:]
        notices = []
        for block in blocks:
            msg = email.message_from_bytes(block, policy=policy.default)
            content = msg.get_body(preferencelist=("plain",)).get_content()
            section = re.search(r"(?ms)^\*\* ASASSN-14iy[^\n]*\n(.*?)(?=^\*\*|\Z)", content)
            if section:
                rows = []
                for date, mag, source in re.findall(r"(?m)^\s*(\d{8}\.\d+)\s+(\S+)\s+\(([^\n]+)\)", section[1]):
                    rows.append([date, mag, source.split(" (")[0]])
                if not rows:
                    raise ValueError("Unrecognized VSNET observation format")
                notices.append({"subject": str(msg["Subject"]), "date": str(msg["Date"]), "measurements": rows})
        if not notices:
            raise ValueError("Missing expected VSNET bulletin")
        excerpt = json.dumps(notices, indent=2).encode()
        (HERE / "october-notices.gz").write_bytes(gzip.compress(excerpt, mtime=0))
        archive.pilot.json_write(receipt, {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "sha256_original_mbox_gz": hashlib.sha256(raw).hexdigest(), "messages_in_archive": len(blocks),
            "sha256_uncompressed": hashlib.sha256(excerpt).hexdigest()})
    return json.loads(archive.source("october-notices", url, False))


def vizier(key, name, ra, dec, fields, fetch):
    params = {"-source": name, "-c": f"{ra} {dec:+}", "-c.rs": 5,
              "-out": fields, "-out.add": "_r", "-out.max": 100}
    url = archive.vsx.SERVICE + "?" + urlencode(params)
    root = ElementTree.fromstring(archive.source(key + ".xml", url, fetch))
    if any(v.get("name", "").lower() == "error" or v.get("value", "").upper() in ("ERROR", "OVERFLOW")
           for v in root.findall(".//{*}INFO")):
        raise ValueError("Incomplete VizieR cone")
    (table,) = root.findall(".//{*}TABLE")
    fields = [el.get("name") for el in table.findall("{*}FIELD")]
    return [dict(zip(fields, [td.text for td in tr], strict=True)) for tr in table.findall(".//{*}TR")]


def run(fetch):
    vsx = {}
    for name in ["ASASSN-14iy", "MGAB-V1042", "ASASSN-16le"]:
        url = "https://vsx.aavso.org/index.php?" + urlencode({"view": "api.object", "ident": name})
        root = ElementTree.fromstring(archive.source(name + "-vsx.xml", url, fetch))
        vsx[name] = {el.tag: el.text for el in root}
        if vsx[name]["Name"] != name:
            raise ValueError("Unexpected VSX name")
    ra, dec = 15.413625, 75.1437291
    gaia = vizier("gaia-onset", "I/355/gaiadr3", ra, dec, "Source,RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag,Plx,RUWE", fetch)
    for row in gaia:
        if row["pmRA"] and row["pmDE"]:
            dt_years = (2460596.8183333 - 2457389.0) / 365.25  # Gaia DR3 J2016.0
            a, d = float(row["RA_ICRS"]), float(row["DE_ICRS"])
            x = a + dt_years * float(row["pmRA"]) / (3600000 * math.cos(math.radians(d)))
            y = d + dt_years * float(row["pmDE"]) / 3600000
            row["predicted_separation_at_alert_arcsec"] = round(math.hypot(*archive.pilot.offsets(x, y, ra, dec)), 4)
    ps1 = vizier("sesar-fading", "J/AJ/153/204/table5", 330.183135, 63.716245,
                 "recno,RAJ2000,DEJ2000,S3ab,S3c,Per,phi0,gAmp,rAmp,iAmp,zAmp,Tg,Tr,Ti,Tz,gmag,rmag,imag,zmag", fetch)
    params = {"view": "api.delim", "ident": "ASASSN-14iy", "fromjd": 2460554.5,
              "tojd": 2460624.5, "delimiter": "@@@"}
    url = "https://vsx.aavso.org/index.php?" + urlencode(params)
    lines = archive.source("onset-observers.txt", url, fetch).decode("utf-8-sig").splitlines()
    fields = lines[0].split("@@@")
    if not {"JD", "mag", "obsID"} <= set(fields):
        raise ValueError("Unknown AAVSO response")
    obs = [dict(zip(fields, line.split("@@@"), strict=True)) for line in lines[1:]]
    return {"vsx": vsx, "gaia_dr3_5arcsec": gaia, "sesar_5arcsec": ps1, "vsnet": october_notices(fetch),
            "aavso_2024_onset_count": len(obs),
            "aavso_2024_onset_dates": sorted({str(archive.pilot.day_of(float(row["JD"]))) for row in obs})}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    print(json.dumps(run(parser.parse_args().fetch), indent=2))
