"""Select the first sufficiently covered month in the preregistered order."""

from collections import Counter
from datetime import date
import json
from pathlib import Path

import pilot

PERIODS = [
    (date(2026, 10, 1), date(2026, 11, 1)),
    (date(2025, 10, 1), date(2025, 11, 1)),
    (date(2024, 11, 1), date(2024, 12, 1)),
]


def coverage(directory):
    directory.mkdir(parents=True, exist_ok=True)
    output = []
    for start, stop in PERIODS:
        folder = directory / start.strftime("%Y-%m")
        folder.mkdir(exist_ok=True)
        receipt = folder / "manifest.json"
        manifest = (
            json.loads(receipt.read_text()) if receipt.exists() else {"sources": {}}
        )
        n = 10 * ((stop - start).days + 2) + 1
        data = pilot.store_raw(
            folder,
            manifest,
            "anomaly.json.gz",
            pilot.FINK + "/anomaly",
            {
                "start_date": str(start),
                "stop_date": str(stop),
                "n": n,
                "columns": "i:objectId,i:jd",
                "output-format": "json",
            },
        )
        if len(data) >= n:
            raise ValueError("archive coverage query saturated its cap")
        monthly = pilot.store_raw(
            folder,
            manifest,
            "statistics.json.gz",
            pilot.FINK + "/statistics",
            {
                "date": start.strftime("%Y%m"),
                "output-format": "json",
            },
        )
        rows = [
            r for r in data if pilot.utc_jd(start) <= r["i:jd"] < pilot.utc_jd(stop)
        ]
        by_day = Counter(str(pilot.day_of(r["i:jd"])) for r in rows)
        ids = {pilot.oid(r) for r in rows}
        passed = len(by_day) >= 10 and len(ids) >= 20
        stats = {r["key:key"]: int(r["basic:sci"]) for r in monthly}
        if any(not s.startswith("ztf_" + start.strftime("%Y%m")) for s in stats):
            raise ValueError("unexpected observing month from Fink statistics")
        output.append(
            {
                "start": str(start),
                "stop_exclusive": str(stop),
                "anomaly_packets": len(rows),
                "unique_objects": len(ids),
                "utc_dates_with_anomalies": len(by_day),
                "anomalies_per_utc_date": dict(sorted(by_day.items())),
                "observing_nights_in_statistics": len(stats),
                "science_alerts_in_statistics": sum(stats.values()),
                "eligible": passed,
            }
        )
        if passed:
            break
    pilot.json_write(directory / "coverage.json", output)
    return output


if __name__ == "__main__":
    for row in coverage(
        Path(__file__).resolve().parent / "experiments/replication/coverage"
    ):
        print(
            row["start"],
            row["anomaly_packets"],
            "packets,",
            row["unique_objects"],
            "objects,",
            row["utc_dates_with_anomalies"],
            "dates;",
            "eligible =",
            row["eligible"],
        )
