"""Plot the 21 visit medians from the archived target and comparison light curves."""

import json
from datetime import datetime, timedelta
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter, YearLocator

ROOT = Path(__file__).resolve().parent
series = json.loads((ROOT/"results.json").read_text())
fig, ax = plt.subplots(figsize=(8,4.5))
fig.set_facecolor("#faf8f2")
ax.set_facecolor("#faf8f2")
for row, name, color, mark in zip(series, ["Target position","Nearby star A","Nearby star B"],
                                 ["#00736c","#815d91","#a96436"],["o","s","^"]):
    visits = row["visits"]
    baseline = median(v["median_w1"] for v in visits)
    times = [datetime(1858,11,17)+timedelta(days=v["median_mjd"]) for v in visits]
    flux = [baseline-v["median_w1"] for v in visits]
    ax.plot(times,flux,linewidth=1.5,linestyle="-" if name=="Target position" else ":",
            marker=mark,color=color,markersize=5,label=name)
ax.axhline(0,color="#83918d",linewidth=.7,zorder=0)
ax.set_title("Infrared brightness varied at the alert position",loc="left",fontsize=16,pad=15)
ax.set_ylabel("Brighter than each source’s median (mag)")
ax.xaxis.set_major_locator(YearLocator(2))
ax.xaxis.set_major_formatter(DateFormatter("%Y"))
ax.spines[["right","top"]].set_visible(False)
ax.legend(loc="upper left",frameon=False,ncol=3,fontsize=9)
ax.margins(x=.025,y=.27)
fig.tight_layout()
fig.savefig(ROOT/"visits.png",dpi=145,facecolor=fig.get_facecolor())
