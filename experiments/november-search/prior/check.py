"""Freeze and summarize complementary public recognition and registry checks."""

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parent
CASES = {"ZTF18abqclig": (92.6614871, 20.6345684),
         "ZTF20aaothpc": (130.6710409, -2.381935),
         "ZTF18abjtfpx": (280.305911, 9.7241705)}


def capture(name, url, fetch, payload=None):
    receipt = BASE / (name + ".meta.json")
    path = BASE / (name + ".gz")
    if fetch and not receipt.exists():
        with NamedTemporaryFile() as tmp:
            args = ["curl", "-sSL", "--max-time", "45", "--max-filesize", "2000000",
                    "-w", "%{http_code}", "-o", tmp.name]
            if payload is not None:
                args += ["-H", "Content-Type: application/json", "-d", json.dumps(payload)]
            status = int(subprocess.check_output([*args, url], text=True))
            raw = tmp.read()
        meta = {"url": url, "json_payload": payload, "http_status": status,
                "sha256": hashlib.sha256(raw).hexdigest(),
                "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        path.write_bytes(gzip.compress(raw, mtime=0))
        receipt.write_text(json.dumps(meta, indent=2) + "\n")
    raw = gzip.decompress(path.read_bytes())
    meta = json.loads(receipt.read_text())
    if meta["url"] != url or meta["json_payload"] != payload or meta["sha256"] != hashlib.sha256(raw).hexdigest():
        raise ValueError(f"Provenance mismatch for {name}")
    return meta["http_status"], raw


def sep(a, b):
    ra1, d1, ra2, d2 = map(math.radians, (*a, *b))
    h = math.sin((d1-d2)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((ra1-ra2)/2)**2
    return round(7200*math.degrees(math.asin(math.sqrt(h))), 3)


def catalogs(case, pos, fetch):
    tables = ["II/360", "J/MNRAS/524/4867/tablea1"]
    if case == "ZTF20aaothpc":
        tables += ["J/A+A/675/A195/candztfg", "J/A+A/675/A195/candztfr"]
    params = {"-source": ",".join(tables), "-c": f"{pos[0]} {pos[1]:+}",
              "-c.rs": "5", "-out.add": "_r", "-out.max": "100",
              "-out": "GaiaDR2,RA_ICRS,DE_ICRS,R,LY,SY,LE,SE,LM,SM,LEG,SEG,SDSS,Name,RAJ2000,DEJ2000,ObjectId,RAdeg,DEdeg,Nepochs,PredIClass,PredIClassProb,PredClass,PredClassProb,Period"}
    status, raw = capture(case + "-vizier", "https://vizier.cfa.harvard.edu/viz-bin/votable?" + urlencode(params), fetch)
    if status != 200:
        return {"http_status": status}
    root = ET.fromstring(raw)
    if any(i.get("name", "").upper() == "ERROR" or i.get("value", "").upper() in ("ERROR", "OVERFLOW")
           for i in root.findall(".//{*}INFO")):
        raise ValueError(raw[-500:])
    result = {}
    for table in root.findall(".//{*}TABLE"):
        name = table.attrib["name"]
        fields = [f.attrib["name"] for f in table.findall("{*}FIELD")]
        rows = []
        for tr in table.findall(".//{*}TR"):
            row = dict(zip(fields, [(c.text or "").strip() for c in tr.findall("{*}TD")], strict=True))
            key = next((pair for pair in [("RA_ICRS", "DE_ICRS"), ("RAdeg", "DEdeg"), ("RAJ2000", "DEJ2000")]
                        if all(row.get(k) for k in pair)), None)
            if key:
                row["calculated_arcsec"] = sep(pos, tuple(float(row[k]) for k in key))
                if abs(float(row["_r"]) - row["calculated_arcsec"]) > 0.105:
                    raise ValueError(f"VizieR separation mismatch: {row}")
            rows.append(row)
        result[name] = rows
    expected = {"II/360/catalog", *tables[1:]}
    missing = expected - result.keys()
    if result.keys() - expected or missing - {"J/A+A/675/A195/candztfg", "J/A+A/675/A195/candztfr"}:
        raise ValueError(f"Unexpected VizieR tables: {result.keys()}")
    for name in sorted(missing):
        result[name] = {"status": "requested table absent from VizieR response"}
    return result


def run(fetch):
    results = {}
    for case, pos in CASES.items():
        item = results[case] = {}
        status, _ = capture(case + "-zvar", "https://api.zvar.caltech.edu/imdiff/cone_search",
                            fetch, {"ra": pos[0], "dec": pos[1], "radius": 5})
        item["zvar_http"] = status
        if status == 200:
            raise ValueError("Successful ZVAR response requires period adjudication")
        if case != "ZTF18abjtfpx":
            item["vizier"] = catalogs(case, pos, fetch)
        api = "https://api.alerce.online/ztf/v1/objects/" + case
        status, raw = capture(case + "-alerce-object", api, fetch)
        if status == 200:
            obj = json.loads(raw)
            if obj["oid"] != case:
                raise ValueError("Wrong ALeRCE object")
            item["alerce_separation_arcsec"] = sep(pos, (obj["meanra"], obj["meandec"]))
        else:
            item["alerce_object_http"] = status
        status, raw = capture(case + "-alerce-class", api + "/probabilities", fetch)
        if status == 200:
            item["alerce_rank1"] = [c for c in json.loads(raw) if c["ranking"] == 1]
        else:
            item["alerce_class_http"] = status
        if case == "ZTF18abjtfpx":
            status, raw = capture(case + "-alerce-features", api + "/features", fetch)
            if status != 200:
                raise ValueError("Timing feature response failed")
            item["alerce_periods"] = [r for r in json.loads(raw) if r["name"] in ("Multiband_period", "Period_band")]
            continue
        query = urlencode({"ra": pos[0], "decl": pos[1], "radius": 10, "coords_unit": "arcsec", "format": "csv"})
        status, raw = capture(case + "-tns", "https://www.wis-tns.org/search?" + query, fetch)
        item["tns_http"] = status
        if status == 200:
            raise ValueError("Successful TNS response requires source adjudication")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    print(json.dumps(run(p.parse_args().fetch), indent=2))
