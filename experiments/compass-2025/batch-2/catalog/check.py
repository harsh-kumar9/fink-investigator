"""Resolve SIMBAD and variable catalog neighbors against the three frozen batch anchors."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
from xml.etree import ElementTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]))
import pilot
import vsx

TABLES = ["B/vsx/vsx", "J/AJ/156/241/table4", "II/366/catv2021", "I/358/vrrlyr"]
GAIA_VARIABLES = {
    "ZTF18aagstwh": ("I/358/vrrlyr", "1403797832926731008", ["PF", "e_PF", "P1O", "Class"]),
    "ZTF18abmfrjj": ("I/358/vlpv", "4239802740358274688", ["Freq", "e_Freq", "Amp"]),
}


def source(name, url, fetch):
    path = HERE / (name + ".gz")
    receipt = HERE / (name + ".json")
    if fetch and not receipt.exists():
        raw = pilot.curl(url, [])
        path.write_bytes(gzip.compress(raw, mtime=0))
        pilot.json_write(receipt, {
            "url": url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
        })
    meta = json.loads(receipt.read_text())
    raw = gzip.decompress(path.read_bytes())
    if meta["url"] != url or hashlib.sha256(raw).hexdigest() != meta["sha256_uncompressed"]:
        raise ValueError(f"Wrong provenance for {name}")
    return raw


def vizier(case, fetch):
    params = {
        "-source": ",".join(TABLES),
        "-c": f"{case['ra']} {case['dec']:+}",
        "-c.rs": 5, "-out.max": 100, "-out.add": "_r",
        "-out": "OID,Name,ATOID,ASASSN-V,Source,RAJ2000,DEJ2000,RA_ICRS,DE_ICRS,Type,Class,Period,Per,fp-LSper,ls-Pday,P1,Pfund,P1O,RRtype",
    }
    raw = source(case["id"] + ".xml", vsx.SERVICE + "?" + urlencode(params), fetch)
    root = ElementTree.fromstring(raw)
    for info in root.findall(".//{*}INFO"):
        if info.get("name", "").lower() == "error" or info.get("value", "").upper() in ("ERROR", "OVERFLOW"):
            raise ValueError(f"VizieR error for {case['id']}")
    result = {table: None for table in TABLES}
    for tab in root.findall(".//{*}TABLE"):
        cols = [f.attrib["name"] for f in tab.findall("{*}FIELD")]
        data = []
        for tr in tab.findall(".//{*}TR"):
            row = dict(zip(cols, [(td.text or "").strip() for td in tr], strict=True))
            if float(row["_r"]) > 5.01:
                raise ValueError("Position outside query")
            data.append(row)
        if len(data) >= 100:
            raise ValueError("Truncated VizieR cone")
        result[tab.attrib["name"]] = data
    return result


def gaia_period(case, fetch):
    table, source_id, measures = GAIA_VARIABLES[case["id"]]
    params = {"-source": table, "Source": source_id, "-out.all": 1, "-out.max": 5}
    raw = source(case["id"] + "-gaia.xml", vsx.SERVICE + "?" + urlencode(params), fetch)
    root = ElementTree.fromstring(raw)
    (tab,) = root.findall(".//{*}TABLE")
    (tr,) = tab.findall(".//{*}TR")
    cols = [f.attrib["name"] for f in tab.findall("{*}FIELD")]
    row = dict(zip(cols, [(td.text or "").strip() for td in tr], strict=True))
    distance = math.hypot(*pilot.offsets(case["ra"], case["dec"], float(row["RA_ICRS"]), float(row["DE_ICRS"])))
    if tab.attrib["name"] != table or row["Source"] != source_id or distance > 1:
        raise ValueError("Unexpected Gaia variability counterpart")
    data = {k: row[k] for k in ["Source", "RA_ICRS", "DE_ICRS", *measures]}
    if "Freq" in data:
        frequency, error = float(data["Freq"]), float(data["e_Freq"])
        data.update(period_days=1 / frequency, linearized_period_error_days=error / frequency**2)
    return data


def run(fetch):
    cases = [r for r in json.loads((HERE.parent / "results/ledger.json").read_text())["objects"] if r["selected_for"]]
    clauses = [f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{a['ra']},{a['dec']},{10/3600}))=1" for a in cases]
    query = "SELECT TOP 100 main_id,otype,ra,dec,pmra,pmdec FROM basic WHERE " + " OR ".join(clauses)
    url = pilot.SIMBAD + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    raw = source("simbad.csv", url, fetch)
    reader = csv.DictReader(io.StringIO(raw.decode().lstrip()))
    if "main_id" not in (reader.fieldnames or []):
        raise ValueError("SIMBAD query failed")
    neighbors = list(reader)
    if len(neighbors) >= 100:
        raise ValueError("Truncated SIMBAD cone")
    result = {}
    for a in cases:
        local = []
        for row in neighbors:
            distance = math.hypot(*pilot.offsets(a["ra"], a["dec"], float(row["ra"]), float(row["dec"])))
            if distance <= 10.01:
                local.append({**row, "separation_arcsec": round(distance, 3)})
        result[a["id"]] = {"simbad_10arcsec": local, "catalogs_5arcsec": vizier(a, fetch)}
        if a["id"] in GAIA_VARIABLES:
            result[a["id"]]["gaia_variability"] = gaia_period(a, fetch)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    args = p.parse_args()
    print(json.dumps(run(args.fetch), indent=2))
