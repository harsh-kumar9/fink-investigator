"""Identify archived May bulletins and VSOLJ observing dates for ASASSN-25ba."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import email
from email import policy
from email.utils import parsedate_to_datetime
import gzip
import hashlib
import json
import re
from urllib.parse import urlencode

from check import HERE, catalog

LISTS = ["vsnet-alert", "vsnet-campaign-dn", "vsnet-obs", "vsnet-outburst", "vsnet-cvcirc"]
ALIASES = ["asassn25ba", "ztf25aaohzdc", "j203927375093336899"]


def bulletins(raw):
    data = gzip.decompress(raw)
    blocks = re.split(br"(?m)^From .+?\n", data)[1:]
    subjects = []
    dates = []
    for block in blocks:
        msg = email.message_from_bytes(block, policy=policy.default)
        part = msg.get_body(preferencelist=("plain",))
        if part is None:
            raise ValueError("Unexpected message format in VSNET mailbox")
        date = parsedate_to_datetime(str(msg["Date"])).astimezone(timezone.utc).date().isoformat()
        dates.append(date)
        body = re.sub(r"[^a-z0-9]", "", (str(msg["Subject"]) + " " + part.get_content()).lower())
        if any(name in body for name in ALIASES):
            subjects.append({"date_sent_utc": date, "subject": str(msg["Subject"])})
    return {
        "sha256_decompressed": hashlib.sha256(data).hexdigest(),
        "messages": len(blocks), "from": min(dates), "through": max(dates), "matching_subjects": subjects,
    }


def coverage(raw):
    rows = [(r[1][:8], r[3]) for line in raw.decode("ascii").splitlines()
            if (r := line.split()) and r[0].lower() == "asassn-25ba"]
    total = re.search(br"Total:\s*(\d+) obs\.", raw)
    if total is None or int(total[1]) != len(rows):
        raise ValueError("VSOLJ inventory does not match response")
    return {"sha256_response": hashlib.sha256(raw).hexdigest(),
            "dates": [{"utc": d, "observer_code": o, "observations": n} for (d, o), n in sorted(Counter(rows).items())]}


def run(fetch=False, verify_live=False):
    path = HERE / "coverage.json"
    if path.exists() and not verify_live:
        return json.loads(path.read_text())
    if not fetch and not verify_live:
        raise FileNotFoundError(path)
    result = {}
    for name in LISTS:
        addr = name + "@ooruri.kusastro.kyoto-u.ac.jp"
        url = f"http://ooruri.kusastro.kyoto-u.ac.jp/mailman3/hyperkitty/list/{addr}/export/{addr}-2025-06.mbox.gz?start=2025-05-01&end=2025-06-01"
        result[name] = {"url": url, **bulletins(catalog.pilot.curl(url, []))}
    query = dict(mode="view", object_name="ASASSN-25ba", sdate="20250301", edate="20250630",
                 observer="", tz="UT-VSNET", asfile="Y")
    url = "http://www.cetus-net.org/cgi-bin/obs_search.cgi?" + urlencode(query)
    result["vsolj"] = {"url": url, **coverage(catalog.pilot.curl(url, []))}
    if verify_live:
        if json.loads(path.read_text())["sources"] != result:
            raise ValueError("Observer coverage has changed")
    else:
        catalog.pilot.json_write(path, {"retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "sources": result})
    return json.loads(path.read_text())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--verify-live", action="store_true")
    args = p.parse_args()
    print(json.dumps(run(args.fetch, args.verify_live), indent=2))
