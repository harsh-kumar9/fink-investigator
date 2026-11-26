"""Audit observer archive coverage of TCP J07222683+6220548. Only inventory VSOLJ data."""

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

from check import HERE, pilot, source

NAME = "TCP J07222683+6220548"
AUID = "000-BQC-247"
VSX = "https://vsx.aavso.org/index.php?"
APRIL = (2460766.5, 2460796.5)


def aid(star, start, stop, name, fetch):
    params = dict(view="api.delim", ident=star, fromjd=start, tojd=stop, delimiter="|")
    raw = source(name + ".aid", VSX + urlencode(params), fetch)
    reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")), delimiter="|", quoting=csv.QUOTE_NONE)
    header = next(reader)
    if not {"starName", "JD", "obsID", "comment"}.issubset(header):
        raise ValueError("AAVSO did not return its data table")
    rows = []
    for fields in reader:
        # The archive uses unquoted pipes inside some observers' comments.
        comment = header.index("comment")
        after_comment = comment + 1 + len(fields) - len(header)
        cells = fields[:comment] + ["|".join(fields[comment:after_comment])] + fields[after_comment:]
        rows.append(dict(zip(header, cells, strict=True)))
    if len({r["obsID"] for r in rows if r["obsID"].isdigit()}) != len(rows):
        raise ValueError("Malformed or repeated AAVSO observation IDs")
    if any(not start <= float(r["JD"]) <= stop for r in rows):
        raise ValueError("AAVSO returned an out-of-window observation")
    if any(r["starName"].upper() not in {NAME.upper(), "ASASSN-25AJ"} for r in rows):
        raise ValueError("Unexpected target in AAVSO response")
    return {
        "observations": len(rows),
        "dates": sorted({str(pilot.day_of(float(r["JD"]))) for r in rows}),
        "returned_names": sorted({r["starName"] for r in rows}),
    }


def vsolj_inventory(start, end, name, fetch, verify_live):
    # The public discovery index asks users to contact VSOLJ before research use
    # of its raw photometry. Keep just the availability of dated records.
    params = dict(mode="view", object_name=NAME.replace(" ", ""), sdate=start,
                  edate=end, observer="", tz="UT-VSNET", asfile="Y")
    url = "http://www.cetus-net.org/cgi-bin/obs_search.cgi?" + urlencode(params)
    path = HERE / (name + ".json")
    if verify_live or (fetch and not path.exists()):
        from datetime import datetime, timezone
        raw = pilot.curl(url, [])
        body = raw.decode("ascii")
        counts = Counter()
        for line in body.splitlines():
            if line.startswith(params["object_name"] + " "):
                fields = line.split()
                date = fields[1][:8]
                if not (len(fields) == 4 and start <= date <= end and fields[3].endswith(".VSOLJ")):
                    raise ValueError("Unexpected VSOLJ observation format")
                counts[(date, fields[3])] += 1
        total = re.search(r"Total:\s*(\d+) obs\.", body)
        if total is None or int(total[1]) != sum(counts.values()):
            raise ValueError("Incomplete VSOLJ coverage extraction")
        receipt = {
            "url": url, "sha256_response": hashlib.sha256(raw).hexdigest(),
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "total_rows": sum(counts.values()),
            "coverage": [{"date_utc": d, "observer_code": o, "rows": n}
                         for (d, o), n in sorted(counts.items())],
        }
        if verify_live:
            saved = json.loads(path.read_text())
            if any(receipt[k] != saved[k] for k in ("url", "sha256_response", "coverage", "total_rows")):
                raise ValueError(f"VSOLJ source changed for {name}")
        else:
            pilot.json_write(path, receipt)
    data = json.loads(path.read_text())
    if data["url"] != url or sum(i["rows"] for i in data["coverage"]) != data["total_rows"]:
        raise ValueError("Unexpected VSOLJ query or counts")
    return data


def run(fetch=False, verify_live=False):
    obj = ElementTree.fromstring(source("vsx.xml", VSX + urlencode(dict(view="api.object", ident=NAME)), fetch))
    if (obj.findtext("Name"), obj.findtext("AUID"), obj.findtext("OID")) != (NAME, AUID, "8278004"):
        raise ValueError("AAVSO identity differs")
    aids = {
        "april_name": aid(NAME, *APRIL, "april-name", fetch),
        "april_auid": aid(AUID, *APRIL, "april-auid", fetch),
        "april_alias": aid("ASASSN-25aj", *APRIL, "april-alias", fetch),
        "active_comparison": aid(NAME, 2460700.5, 2460725.5, "active", fetch),
    }
    coverage = {
        label: vsolj_inventory(start, stop, "vsolj-" + label, fetch, verify_live)
        for label, start, stop in [("january", "20250117", "20250121"),
                                   ("april", "20250401", "20250430")]
    }
    return {"aavso": aids, "vsolj": coverage}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--verify-live", action="store_true", help="Refetch and verify both VSOLJ response hashes and coverage")
    args = p.parse_args()
    print(json.dumps(run(args.fetch, args.verify_live), indent=2))
