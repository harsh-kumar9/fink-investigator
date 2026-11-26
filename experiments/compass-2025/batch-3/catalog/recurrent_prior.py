"""Check the photometric QSO candidate label against its source and two spectroscopy archives."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from urllib.parse import urlencode
from xml.etree import ElementTree

from check import HERE, catalog

CATALOGS = ["VII/283/catalog", "VII/290/catalog", "VII/294/catalog", "J/ApJS/219/39/cand"]
PLOT = "https://vsx.aavso.org/index.php?view=docs.open&fid=25415"


def run(fetch=False):
    position = "99.526908 +45.2491446"
    params = {"-source": ",".join(CATALOGS), "-c": position, "-c.rs": 10, "-out.all": 1, "-out.max": 20}
    raw = catalog.source("quasar-priors.xml", catalog.vsx.SERVICE + "?" + urlencode(params), fetch)
    root = ElementTree.fromstring(raw)
    catalog_rows = {}
    for t in root.findall(".//{*}TABLE"):
        name = t.get("name")
        cols = [f.get("name") for f in t.findall("{*}FIELD")]
        fields = (["_r", "Name", "Cl+Ass", "z", "rName", "rz", "Qpct", "XName", "RName"] if name.startswith("VII/")
                  else ["_r", "RAJ2000", "DEJ2000", "Class", "zsp", "zphO", "b_zphO", "B_zphO"])
        rows = [dict(zip(cols, [v.text for v in tr], strict=True)) for tr in t.findall(".//{*}TR")]
        if len(rows) >= 20 or any(float(r["_r"]) > 10 for r in rows):
            raise ValueError("Quasar cone exceeds limits")
        catalog_rows[name] = [{k: row[k] for k in fields} for row in rows]
    if set(catalog_rows) != set(CATALOGS):
        raise ValueError("A quasar catalog table was unavailable")
    query = ("SELECT TOP 20 s.specObjID, s.ra, s.dec, s.class,s.z,s.zWarning FROM "
             "dbo.fGetNearbySpecObjEq(99.526908,45.2491446,0.166667) AS n JOIN SpecObj AS s ON s.specObjID=n.specObjID")
    url = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?" + urlencode({"cmd": query, "format": "csv"})
    sdss = catalog.source("sdss-dr17.csv", url, fetch).decode().splitlines()
    reader = csv.DictReader(line for line in sdss if not line.startswith("#"))
    if reader.fieldnames != ["specObjID", "ra", "dec", "class", "z", "zWarning"]:
        raise ValueError("SDSS did not return a spectrum table")
    spectra = {"sdss_dr17_primary": list(reader)}
    query = ("select targetid,mean_fiber_ra,mean_fiber_dec,spectype,z,zwarn from desi_dr1.zpix where "
             "q3c_radial_query(mean_fiber_ra,mean_fiber_dec,99.526908,45.2491446,0.002777778) limit 20")
    url = "https://datalab.noirlab.edu/query/query?" + urlencode({"sql": query, "ofmt": "csv", "async": "False"})
    receipt = HERE / "desi-dr1.csv.json"
    if fetch and not receipt.exists():
        # The NOIRLab query manager documents this fixed token for anonymous access.
        raw = catalog.pilot.curl(url, ["-H", "X-DL-AuthToken: anonymous.0.0.anon_access"])
        (HERE / "desi-dr1.csv.gz").write_bytes(gzip.compress(raw, mtime=0))
        catalog.pilot.json_write(receipt, {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                         "sha256_uncompressed": hashlib.sha256(raw).hexdigest()})
    raw = catalog.source("desi-dr1.csv", url, False)
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if reader.fieldnames != ["targetid", "mean_fiber_ra", "mean_fiber_dec", "spectype", "z", "zwarn"]:
        raise ValueError("DESI did not return a spectrum table")
    spectra["desi_dr1_zpix"] = list(reader)
    if any(len(v) >= 20 for v in spectra.values()):
        raise ValueError("Spectroscopic cone truncated")
    return {"quasar_catalogs": catalog_rows, "spectra_within_10arcsec": spectra, "discovery_plot": PLOT}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    print(json.dumps(run(p.parse_args().fetch), indent=2))
