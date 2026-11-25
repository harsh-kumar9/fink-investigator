"""Replay the AllWISE counterpart and NEOWISE visit comparison from archived IRSA responses."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from statistics import median
import subprocess
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parent
ANCHOR = (87.2258635, 26.8631435)
BASE = "https://irsa.ipac.caltech.edu/TAP/sync"


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def separation(ra, dec, center):
    r1, d1, r2, d2 = map(math.radians,(ra,dec,*center))
    s = math.sin((d1-d2)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r1-r2)/2)**2
    return 2*math.degrees(math.asin(math.sqrt(s)))*3600


def pick():
    target, = read(ROOT / "target-cone.csv")
    pool = []
    for star in read(ROOT / "field.csv"):
        distance = separation(float(star["ra"]),float(star["dec"]),ANCHOR)
        if (star["w1mpro"] and abs(float(star["w1mpro"])-float(target["w1mpro"]))<=0.5
                and float(star["w1snr"] or 0)>=20 and star["cc_flags"][:2]=="00"
                and star["nb"]=="1" and star["ext_flg"]=="0" and distance>=15):
            pool.append((distance,star["designation"],star))
    if len(pool)<2:
        raise ValueError("Insufficient field controls")
    return [target]+[s for _,_,s in sorted(pool)[:2]]


def query(star):
    params = {"catalog":"neowiser_p1bs_psd","spatial":"cone","objstr":f"{star['ra']} {star['dec']}",
              "radius":"2.5","radunits":"arcsec","outfmt":"1",
              "selcols":"mjd,ra,dec,w1mpro,w1sigmpro,w1snr,cc_flags,moon_masked,qual_frame,nb"}
    return "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?"+urlencode(params)


def read_ipac(path):
    source = path.read_bytes()
    if path.suffix == ".gz":
        source = gzip.decompress(source)
    lines = source.decode().splitlines()
    header = next(line for line in lines if line.startswith("|"))
    boundaries = [i for i,c in enumerate(header) if c == "|"]
    columns = [header[a+1:b].strip() for a,b in zip(boundaries,boundaries[1:])]
    data = []
    for line in lines:
        if not line.strip() or line.startswith((chr(92),"|")):
            continue
        values = [line[a+1:b+1].strip() for a,b in zip(boundaries,boundaries[1:])]
        data.append(dict(zip(columns,["" if v.lower()=="null" else v for v in values])))
    count = next(line for line in lines if line.startswith(chr(92)+"RowsRetrieved"))
    if len(data) != int(count.split("=")[1]):
        raise ValueError(f"Incomplete exposure rows in {path}")
    return data


def quantile(values, frac):
    vals = sorted(values)
    idx = (len(vals)-1)*frac
    gap = idx-math.floor(idx)
    return vals[math.floor(idx)]*(1-gap)+vals[math.ceil(idx)]*gap


def summarize(star, role, fetch):
    path = ROOT/(star["designation"]+".tbl.gz")
    if fetch and not path.exists():
        tmp = path.with_suffix(".tmp")
        subprocess.run(["curl","-fsSL","--max-time","90","--retry","1",query(star),"-o",str(tmp)],check=True)
        if not read_ipac(tmp):
            raise ValueError("Unexpected archive response")
        path.write_bytes(gzip.compress(tmp.read_bytes(),mtime=0))
        tmp.unlink()
    raw = read_ipac(path)
    center = (float(star["ra"]), float(star["dec"]))
    if any(separation(float(point["ra"]), float(point["dec"]), center) > 2.51 for point in raw):
        raise ValueError(f"Exposure outside the requested cone of {star['designation']}")
    good = []
    for p in raw:
        if (p["w1mpro"] and p["w1sigmpro"] and 0<float(p["w1sigmpro"])<=0.15
                and float(p["w1snr"] or 0)>=10 and p["cc_flags"][:1]=="0"
                and p["moon_masked"][:1]=="0" and int(p["qual_frame"])>0 and p["nb"]=="1"):
            good.append((float(p["mjd"]),float(p["w1mpro"])))
    visits = []
    for day,mag in sorted(good):
        if not visits or day-visits[-1][-1][0]>30:
            visits.append([])
        visits[-1].append((day,mag))
    retained = [{"median_mjd":round(median([t for t,_ in v]),3),"n":len(v),
                 "median_w1":median([mag for _,mag in v])} for v in visits if len(v)>=5]
    mags = [v["median_w1"] for v in retained]
    metrics = {}
    if len(mags)>=2:
        metrics = {"range_mag":round(max(mags)-min(mags),3),
                   "p90_minus_p10_mag":round(quantile(mags,.9)-quantile(mags,.1),3)}
    return {"role":role,"designation":star["designation"],"query":query(star),"w1_catalog":star["w1mpro"],
            "separation_arcsec":round(separation(float(star["ra"]),float(star["dec"]),ANCHOR),3),
            "rows":len(raw),"good_rows":len(good),"visits":retained,**metrics}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch",action="store_true")
    args = parser.parse_args()
    result = [summarize(star,"target" if i==0 else "field control",args.fetch) for i,star in enumerate(pick())]
    print(json.dumps(result,indent=2))
