"""Replay the archived VizieR cone for the October timing lead."""

import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
RA, DEC = 280.305911, 9.7241705
GAIA = "I/355/gaiadr3"
ATLAS = "J/AJ/156/241/table4"
ASASSN = "II/366/catv2021"
LAV = "J/A+A/648/A44/tabled1"
URL = "https://vizier.cfa.harvard.edu/viz-bin/votable?" + urlencode({
    "-source": ",".join([GAIA, ATLAS, ASASSN, LAV]),
    "-c": f"{RA} +{DEC}", "-c.rs": 15, "-out.all": 1, "-out.add": "_r", "-out.max": 50,
})


def separation(ra, dec):
    x = math.radians(ra - RA)
    y = math.radians(dec - DEC)
    a = math.sin(y/2)**2 + math.cos(math.radians(DEC)) * math.cos(math.radians(dec)) * math.sin(x/2)**2
    return math.degrees(2*math.asin(math.sqrt(a))) * 3600


def run():
    raw = gzip.decompress((HERE / "cone.xml.gz").read_bytes())
    receipt = json.loads((HERE / "source.json").read_text())
    if receipt["url"] != URL or hashlib.sha256(raw).hexdigest() != receipt["sha256_uncompressed"]:
        raise ValueError("VizieR source or request mismatch")
    root = ET.fromstring(raw)
    if any(i.get("name", "").lower() == "error" or i.get("value", "").upper() in {"ERROR","OVERFLOW"}
           for i in root.findall(".//{*}INFO")):
        raise ValueError("VizieR reported incomplete results")
    data = {}
    for table in root.findall(".//{*}TABLE"):
        fields = [f.attrib["name"] for f in table.findall("{*}FIELD")]
        values = [[(td.text or "").strip() for td in tr] for tr in table.findall(".//{*}TR")]
        data[table.attrib["name"]] = [dict(zip(fields, row, strict=True)) for row in values]
    if set(data) != {GAIA, ATLAS, ASASSN, LAV} or any(len(rows) >= 50 for rows in data.values()):
        raise ValueError("Missing or truncated catalog table")
    atlas, = data[ATLAS]
    offset = separation(float(atlas["RAJ2000"]), float(atlas["DEJ2000"]))
    if abs(offset - float(atlas["_r"])) > 0.05:
        raise ValueError("ATLAS coordinate disagrees with reported separation")
    timing = json.loads((HERE.parent / "results.json").read_text())
    target, = [r for r in timing["period_search"] if r["id"] == "ZTF18abjtfpx"]
    gaia_fields = ("Source", "RA_ICRS", "DE_ICRS", "_r", "Gmag", "pmRA", "pmDE", "VarFlag")
    atlas_fields = ("ATOID", "RAJ2000", "DEJ2000", "_r", "Class", "fp-c-pts", "fp-o-pts", "fp-LSper",
                    "fp-lngfitper", "fp-domper-c", "fp-domper-o", "ls-Pday", "ddc", "prox",
                    "ps-dist", "ps-dist0", "objID")
    return {
        "anchor_deg": [RA, DEC], "radius_arcsec": 15, "rows": {k:len(v) for k,v in data.items()},
        "atlas": {k:atlas[k] for k in atlas_fields},
        "computed_atlas_offset_arcsec": round(offset, 3),
        "relative_period_difference_to_atlas_ls_percent": round(100 * abs(target["period_days"] /
                                                                   float(atlas["fp-LSper"]) - 1), 3),
        "gaia_within_4arcsec_j2016": [
            {k:g[k] for k in gaia_fields} for g in sorted(data[GAIA], key=lambda r: float(r["_r"]))
            if float(g["_r"]) <= 4
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
