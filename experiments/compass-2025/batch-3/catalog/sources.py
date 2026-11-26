"""Archive the discovery survey row, VSX identifications, Chen catalog and AAVSO series."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

from check import HERE, catalog

VSX = "https://vsx.aavso.org/index.php?"
NAME = "ASASSN-25ba"
ASAS = "https://www.astronomy.ohio-state.edu/asassn/transients.html"


def sky_patrol(fetch=False, read_table=False):
    columns = ["asas_sn_id", "ra_deg", "dec_deg", "catalog_sources"]
    url = "http://asassn-lb01.ifa.hawaii.edu:9006/lookup_cone/radius0.0016666666666666668_ra309.8640412_dec-9.5602805"
    post = {"catalog": "master_list", "cols": columns, "format": "arrow", "download": True}
    receipt = HERE / "sky-patrol.json"
    path = HERE / "sky-patrol-cone.parquet"
    if fetch and not receipt.exists():
        import pyarrow.parquet as pq
        raw = catalog.pilot.curl(url, ["-H", "Content-Type: application/json", "--data", json.dumps(post)])
        table = pq.read_table(io.BytesIO(raw))
        if not set(columns).issubset(table.column_names):
            raise ValueError("Unexpected Sky Patrol cone fields")
        path.write_bytes(raw)
        catalog.pilot.json_write(receipt, {
            "url": url, "post": post, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256_response": hashlib.sha256(raw).hexdigest(), "radius_arcsec": 6,
            "returned_columns": table.column_names, "rows": table.num_rows,
        })
    result = json.loads(receipt.read_text())
    raw = path.read_bytes()
    if (result["url"] != url or result["post"] != post
            or hashlib.sha256(raw).hexdigest() != result["sha256_response"]):
        raise ValueError("Mismatched Sky Patrol cone provenance")
    if read_table:
        import pyarrow.parquet as pq
        table = pq.read_table(io.BytesIO(raw))
        if result["rows"] != table.num_rows or result["returned_columns"] != table.column_names:
            raise ValueError("Sky Patrol archive and receipt disagree")
    return result


def aavso(fetch=False):
    params = dict(view="api.delim", ident=NAME, fromjd=2460735.5, tojd=2460857.5, delimiter="@@@")
    raw = catalog.source("onset.aid", VSX + urlencode(params), fetch)
    lines = raw.decode("utf-8-sig").splitlines()
    fields = lines[0].split("@@@")
    rows = [dict(zip(fields, row.split("@@@"), strict=True)) for row in lines[1:]]
    if len({r["obsID"] for r in rows}) != len(rows) or any(
        r["starName"].upper() != NAME.upper() or not params["fromjd"] <= float(r["JD"]) <= params["tojd"]
        for r in rows
    ):
        raise ValueError("Mismatched AAVSO source, repeated observation or date outside requested window")
    return rows


def discovery(fetch):
    # The continuously growing full survey table includes thousands of unrelated records.
    path = HERE / "discovery.json"
    saved_row = HERE / "discovery-row.html.gz"
    if fetch and not path.exists():
        raw = catalog.pilot.curl(ASAS, [])
        rows = re.findall(r"<tr\b[^>]*>(.*?)</tr>", raw.decode("utf-8", errors="replace"), re.S | re.I)
        matched = [r for r in rows if re.search(r"<td>\s*ASASSN-25ba\s*</td>", r)]
        (record,) = matched
        saved_row.write_bytes(gzip.compress(record.encode("utf-8"), mtime=0))
        catalog.pilot.json_write(path, {
            "url": ASAS, "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256_full_response": hashlib.sha256(raw).hexdigest(),
            "sha256_row": hashlib.sha256(record.encode("utf-8")).hexdigest(),
        })
    data = json.loads(path.read_text())
    record = gzip.decompress(saved_row.read_bytes())
    cells = [re.sub(r"<[^>]+>", "", v).strip() for v in re.findall(r"<td>\s*(.*?)\s*</td>", record.decode(), re.S)]
    if data["url"] != ASAS or data["sha256_row"] != hashlib.sha256(record).hexdigest() or cells[0] != NAME:
        raise ValueError("Unexpected ASAS-SN discovery row")
    return {**data, "name": cells[0], "other_name": cells[1], "ra": cells[3], "dec": cells[4],
            "discovery_utc": cells[5], "discovery_mag": float(cells[6]), "survey_note": cells[-1]}


def run(fetch=False, verify_parquet=False):
    names = {}
    for star in [NAME, "MGAB-V243", "ZTF J192701.28+073352.1"]:
        url = VSX + urlencode(dict(view="api.object", ident=star))
        raw = catalog.source("vsx-" + star.split()[0], url, fetch)
        obj = ElementTree.fromstring(raw)
        if obj.findtext("Name") != star:
            raise ValueError("Wrong VSX identification")
        names[star] = {e.tag: e.text for e in obj}
    params = {
        "-source": "J/ApJS/249/18/table2,J/ApJS/249/18/table3",
        "-c": "291.7553476 +7.5644761", "-c.rs": 5, "-out.all": 1, "-out.max": 5,
    }
    raw = catalog.source("chen.xml", catalog.vsx.SERVICE + "?" + urlencode(params), fetch)
    (tr,) = ElementTree.fromstring(raw).findall(".//{*}TR")
    tb = next(tab for tab in ElementTree.fromstring(raw).findall(".//{*}TABLE") if tab.findall(".//{*}TR"))
    row = dict(zip([f.get("name") for f in tb.findall("{*}FIELD")], [v.text for v in tr], strict=True))
    if row["ID"] != "ZTFJ192701.28+073352.1" or float(row["_r"]) > 1:
        raise ValueError("Unexpected source in Chen catalog")
    counts = Counter(str(catalog.pilot.day_of(float(r["JD"]))) for r in aavso(fetch))
    return {"discovery": discovery(fetch), "vsx": names, "chen": row, "aavso_obs_per_date": dict(sorted(counts.items())),
            "sky_patrol_6arcsec": sky_patrol(fetch, verify_parquet)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--verify-parquet", action="store_true", help="Read the Sky Patrol cone table with PyArrow")
    args = p.parse_args()
    print(json.dumps(run(args.fetch, args.verify_parquet), indent=2))
