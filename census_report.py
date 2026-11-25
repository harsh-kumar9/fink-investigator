"""Reconcile the entire frozen month and publish evidence for discriminating follow-up."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import census
import census_followup
import census_photometry
import pilot
import vsx

# Visual observations of the archived earliest, middle and anchor stamps.
IMAGE_READINGS = {
    "ZTF17aaaaitr": (
        "A compact central residual recurs on all three dates. Two Gaia stars lie inside 5″.",
        "Deblend both Gaia positions on sharp images and measure which changes.",
    ),
    "ZTF17aacydyv": (
        "A central residual recurs in 2019, 2023 and 2026. A separate edge residual is visible.",
        "Measure the centroid at bright and faint phases with higher angular resolution and obtain multicolor or spectral data.",
    ),
    "ZTF18aaawbgm": (
        "A central residual recurs beside a second residual to the right. The selected anchor is positive after earlier negative detections.",
        "Register the images to Gaia and compare the centroids and signs of the two residuals across seeing conditions.",
    ),
    "ZTF18abjtfpx": (
        "Three Gaia stars lie inside 5″. The central signal is weaker in the earliest and latest stamps.",
        "Fit all three source positions jointly to test whether a neighbor accounts for the light curve.",
    ),
    "ZTF18ablrppd": (
        "Two Gaia stars lie inside 5″. Dark structure adjoins the central residual. The anchor's broker lookup failed and its DRB is low.",
        "Fit both stars on independent high-resolution images and check whether variability shifts toward the fainter neighbor.",
    ),
    "ZTF18abuizsi": (
        "The central negative signal has adjacent structure in two of the selected stamps.",
        "Test whether flux persists in higher-resolution photometry after control for seeing and reference selection.",
    ),
}

PUBLIC = "https://github.com/harsh-kumar9/fink-investigator/blob/main/experiments/october-census/"


def summarize(out):
    prev, more = census.assess(False)
    extended = {r["id"]: r for r in census_followup.run(False)}
    phot = census_photometry.run(False)
    if not prev["gate_passed"] or not phot["calibration_passed"]:
        raise ValueError("Review gate failed")
    measures = {r["id"]: r for r in phot["exploration"]}
    old = [dict(r) for r in prev["objects"]]
    vsx.annotate_cases(
        census.BASE / "crossmatch-audit",
        [
            r
            for r in old
            if r["status"] == "broker crossmatch failed" and not r["simbad_neighbors"]
        ],
    )
    groups = {"initial sample": old, "remainder": more}
    members = [r["id"] for rs in groups.values() for r in rs]
    coverage = next(
        r
        for r in json.loads((census.BASE.parent / "coverage/coverage.json").read_text())
        if r["start"] == prev["start"]
    )
    if len(set(members)) != len(members) or len(members) != coverage["unique_objects"]:
        raise ValueError("Month membership overlaps or disagrees with coverage")
    sn = json.loads((census.CENSUS / "sources/registry/sn-2026afgd.json").read_text())
    outcomes, ledger, packets = [], [], []
    for cohort, group in groups.items():
        for r in group:
            key = r["id"]
            matches = r.get("vsx_matches")
            simbad = r["simbad_neighbors"]
            gaia = (
                extended.get(key, {})
                .get("tables", {})
                .get(census_followup.TABLES[0], [])
            )
            gaia = sorted(
                [g for g in gaia if g["at_alert_epoch_arcsec"] is not None],
                key=lambda g: g["at_alert_epoch_arcsec"],
            )
            if key in measures:
                d = measures[key]
                star = gaia[0]
                gap = d["negative_minus_positive_median_direct_mag"]
                rho = d["spearman"]
                direct = f"{d['matched_dates']} matched r-band dates ({d['positive_dates']} positive, " f"{d['negative_dates']} negative), from {d['first_matched_date']} to {d['last_matched_date']}. " + (
                    f"Spearman rho {rho}. "
                    if rho is not None
                    else "Rank correlation unavailable. "
                ) + (
                    f"Negative minus positive median direct magnitude {gap} mag."
                    if gap is not None
                    else "Fewer than five dates of one sign, so the median gap is withheld."
                )

                status = (
                    "positional variability supported; class open"
                    if d["threshold_supported"]
                    else "position linked; variability test inconclusive"
                )
                observed, next_test = IMAGE_READINGS[key]
                evidence = [
                    {
                        "id": "gaia-position",
                        "source": extended[key]["query"],
                        "detail": f"Gaia DR3 {star['Source']} lies {star['at_alert_epoch_arcsec']} arcsec from the anchor at its epoch. "
                        f"{len(gaia)} Gaia sources with measured motion in the 5-arcsec cone. No variable-class entry in Gaia DR3 or the two Chen et al. tables.",
                    },
                    {"id": "direct-photometry", "source": d["query"], "detail": direct},
                    {
                        "id": "image-comparison",
                        "source": PUBLIC + f"review/{key}/images/gallery.md",
                        "detail": observed,
                    },
                ]
                if key == "ZTF17aacydyv":
                    evidence.extend(
                        [
                            {
                                "id": "independent-reading",
                                "source": PUBLIC
                                + f"review/{key}/independent-review.md",
                                "detail": "The independent gpt-6-luna review favored recurrent variability at the position and kept a blended neighbor open. This opinion used the same archived optical evidence.",
                            },
                            {
                                "id": "infrared",
                                "source": PUBLIC + "infrared/README.md",
                                "detail": "Independent NEOWISE W1 photometry spanned 0.295 mag between the 10th and 90th percentiles of 21 visit medians. Two comparison stars spanned 0.032 and 0.075 mag. A roughly 6-arcsec beam and unmatched optical/IR phases leave the varying component unresolved.",
                            },
                        ]
                    )
                alternatives = [
                    "A variable at the Gaia position",
                    "A blended neighbor or shared seeing/subtraction effects",
                ]
                p = pilot.case_packet(r)
                p["question"] = next_test
                p["evidence"][-1]["detail"] = "No VSX result within 5 arcsec."
                p["evidence"].extend(evidence)
                packets.append(p)
                ledger.append(
                    {
                        "id": key,
                        "method_outputs": evidence,
                        "competing_interpretations": alternatives,
                        "discriminating_evidence": "A matching signed flux change in direct images supports variability at the position. A centroid that follows another star with brightness, or a correlation with seeing, would favor a blend or instrumental origin.",
                        "status": status,
                        "next_discriminating_observation": next_test,
                    }
                )
            elif key == sn["ztf_id"]:
                status = "registered type Ia supernova"
                evidence = [
                    {
                        "id": "broker",
                        "source": pilot.PORTAL + key,
                        "detail": "No TNS name in the selected packet. SIMBAD and VSX returned empty cones.",
                    },
                    {
                        "id": "registry",
                        "source": sn["tns"]["url"],
                        "detail": "TNS reports SN Ia, with a P60/SEDM spectrum on 2026-10-14.",
                    },
                    {
                        "id": "identifier",
                        "source": sn["lions"]["url"],
                        "detail": "LIONS lists the exact ZTF identifier alongside SN 2026afgd.",
                    },
                ]
                ledger.append(
                    {
                        "id": key,
                        "method_outputs": evidence,
                        "competing_interpretations": [
                            "Previously unreported transient",
                            "Registered supernova with missing broker metadata",
                        ],
                        "status": status,
                        "discriminating_evidence": "An explicit ZTF identifier in the independent registry-linked catalog and a classified spectrum can distinguish earlier reporting from an unregistered transient.",
                        "next_discriminating_observation": "The exact external alias and reported spectrum resolve the registry question. Check updated TNS spectra for any proposed subtype revision.",
                    }
                )
            elif r["gaia_matches"]:
                status = "moving-star explanation supported"
                ledger.append(
                    {
                        "id": key,
                        "method_outputs": [
                            {
                                "id": "motion",
                                "source": PUBLIC + f"review/{key}/reference/report.md",
                                "detail": "G 143-8 crossed within 0.059 arcsec near 2018-07-21, inside the actual 15-image reference interval. All 325 usable detections are negative.",
                            },
                            {
                                "id": "images",
                                "source": PUBLIC + f"review/{key}/images/gallery.md",
                                "detail": "A paired residual appears across three epochs with increasing apparent separation.",
                            },
                        ],
                        "competing_interpretations": [
                            "Subtraction at an old stellar position",
                            "Superposed variable",
                        ],
                        "status": status,
                        "discriminating_evidence": "A crossing during the measured coadd interval and increasing distance between residual lobes are expected when the star moves away from its former position.",
                        "next_discriminating_observation": "Register multi-epoch sharp images and measure both lobes relative to the Gaia track.",
                    }
                )
            elif matches and vsx.closed_alias(r):
                status = "VSX counterpart"
            elif simbad:
                status = "SIMBAD counterpart"
            else:
                status = "insufficient quality or positive dates"
            outcomes.append(
                {
                    "id": key,
                    "cohort": cohort,
                    "candid": r["candid"],
                    "jd": r["jd"],
                    "anchor_ra_deg": r["ra"],
                    "anchor_dec_deg": r["dec"],
                    "score": r["score"],
                    "anchor_drb": r["drb"],
                    "residual": (
                        "positive" if r["isdiffpos"] in ("t", "1", 1) else "negative"
                    ),
                    "usable_dates": r["usable_nights"],
                    "positive_dates": r["positive_nights"],
                    "broker_cds": r["cdsxmatch"],
                    "packet_tns": r["tns"] or "not supplied",
                    "first_pass": r["status"],
                    "simbad_name": "|".join(m["name"] for m in simbad),
                    "simbad_separation_arcsec": "|".join(
                        str(m["separation_arcsec"]) for m in simbad
                    ),
                    "vsx_check": (
                        "not queried"
                        if matches is None
                        else ("matched" if matches else "empty")
                    ),
                    "vsx_oid": "|".join(m["oid"] for m in matches or []),
                    "vsx_name": "|".join(m["name"] for m in matches or []),
                    "vsx_type": "|".join(m["type"] for m in matches or []),
                    "vsx_separation_arcsec": "|".join(
                        str(m["separation_arcsec"]) for m in matches or []
                    ),
                    "gaia_extended_source": gaia[0]["Source"] if gaia else "",
                    "gaia_extended_separation_arcsec": (
                        gaia[0]["at_alert_epoch_arcsec"] if gaia else ""
                    ),
                    "motion_source_id": "|".join(
                        m["source_id"] for m in r["gaia_matches"]
                    ),
                    "track_offset_arcsec": "|".join(
                        str(m["distance_to_track_arcsec"]) for m in r["gaia_matches"]
                    ),
                    "disposition": status,
                    "object_url": pilot.PORTAL + key,
                }
            )
    out.mkdir(parents=True, exist_ok=True)
    with (out / "outcomes.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=outcomes[0], lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(outcomes, key=lambda r: r["id"]))
    counts = {
        "archive_packets": prev["anomaly_alerts"],
        "archive_ids": len(members),
        "anomaly_dates": coverage["utc_dates_with_anomalies"],
        "observing_nights": coverage["observing_nights_in_statistics"],
        "science_alerts": coverage["science_alerts_in_statistics"],
        "moving_star_screen_controls": {
            "true_positives": sum(
                r["expected"] and r["pass"] for r in prev["calibration"]
            ),
            "true_negatives": sum(
                not r["expected"] and r["pass"] for r in prev["calibration"]
            ),
            "total": len(prev["calibration"]),
        },
        "cohort_sizes": {k: len(v) for k, v in groups.items()},
        "first_pass_by_cohort": {
            k: dict(Counter(r["status"] for r in v)) for k, v in groups.items()
        },
        "first_pass": dict(Counter(r["first_pass"] for r in outcomes)),
        "final": dict(Counter(r["disposition"] for r in outcomes)),
        "by_first_pass": {
            s: dict(Counter(r["disposition"] for r in outcomes if r["first_pass"] == s))
            for s in sorted(set(r["first_pass"] for r in outcomes))
        },
    }
    pilot.json_write(out / "summary.json", counts)
    pilot.json_write(out / "discrepancies.json", ledger)
    # Keep two deterministic three-case review packets, preserving the agent review budget.
    for i in range(0, len(packets), 3):
        pilot.json_write(
            out / f"review_packet_{i//3+1}.json",
            {"gate_passed": True, "cases": packets[i : i + 3]},
        )
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=census.CENSUS / "results")
    print(json.dumps(summarize(parser.parse_args().out), indent=2))
