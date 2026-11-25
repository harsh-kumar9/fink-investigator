# Fink anomaly investigator

This is an experiment in iteratively building pipelines for open-world discovery with frontier language models. Agents propose explanations for unusual telescope alerts, check them against public measurements and catalogs, and record what ruled each explanation in or out. Each new search starts with a test on objects whose identities are already known.

Can an investigator recognize a known source of false alarms, then apply that check to historical telescope alerts? This pilot looks for **subtraction residuals from fast-moving foreground stars** in the Zwicky Transient Facility stream. Fink selects the anomaly alerts. This repository tests their interpretation.

The first run investigated 20 of the 144 object IDs in 210 October 2024 anomaly packets. It recovered the exact Gaia identifiers of four published moving-star controls and retained three known supernovae. In the exploration sample, ZTF18abgobyb sits on the past track of the star VB 10. Its predicted crossing date falls inside the exposure interval for the actual reference image. Three other initially unexplained sources all had subarcsecond VSX counterparts. [Cohort and case table](results/report.md) · [VB 10 reference-date check](results/moving-star/report.md) · [Three epochs of VB 10 stamps](review/ZTF18abgobyb/gallery.md) · [Variable-star stamps](review/ZTF17aaagxae/gallery.md) · [Independent image inspection](review/inspection.md)

## Replay the pilot

Run from this directory on Linux or macOS with Python 3.10 or later. The replay uses the standard library and the checked-in archive responses.

```sh
make verify
python3 pilot.py replay --out results
python3 motion.py ZTF18abgobyb --out results/moving-star
python3 agents/review.py plan --limit 3
```

The first command runs tests, reproduces the generated analysis, report and review packet, and compares them with the saved outputs. The second command exposes the same pipeline directly. The third checks the reference exposure dates and the SIMBAD identity against their separately archived responses. The fourth shows which cases agents can investigate. [Agent roles, budgets and cache](agents/README.md)

## Calibration and sampling

1. Freeze the labels in `data/controls.json`. Four distinct physical stars come from Table 4 of Pruzhinskaya et al. G 260-1 appears under two ZTF identifiers in that table and contributes one control. SN 2023mtp, SN 2023ixf and SN 2023bee are negative controls for moving-star contamination. The classifier receives a sky position, observation time and Gaia measurements. It does not receive a star name, label or expected Gaia ID. The answer key compares the result to the published identification and its SIMBAD Gaia DR3 alias. [Fink follow-up study](https://arxiv.org/html/2603.29511v1) · [SN 2023ixf](https://www.wis-tns.org/object/2023ixf) · [SN 2023bee](https://www.wis-tns.org/object/2023bee) · [SIMBAD identifier lookup](https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT%20ids.id%2Cgaia.id%20FROM%20ident%20AS%20ids%20JOIN%20ident%20AS%20gaia%20ON%20ids.oidref%3Dgaia.oidref%20WHERE%20ids.id%20IN%20(%27G%20260-1%27%2C%27G%20135-67%27%2C%27LP%20624-55%27%2C%27Wolf%20611%27)%20AND%20gaia.id%20LIKE%20%27Gaia%20DR3%20%25%27).
2. For each control use its most recent Fink detection before 2026-02-01 UTC. Query Gaia DR3 around the alert centroid within 60″. Require total proper motion of at least 150 mas/year and at least five times the quadrature sum of its two reported component errors. Propagate the Gaia reference position from its catalog epoch (J2016.0 here) to the alert observation time. Accept a position within 2.5″ of the intervening line segment. A star can leave a residual near its old or current position. The calculation uses the local tangent plane and catalog `pmra`, which includes cos δ. Parallax, perspective acceleration and positional covariance are omitted. The search cone also limits sensitivity to the fastest stars. [Gaia DR3 fields](https://gaia.aip.de/metadata/gaiadr3/gaia_source/).
3. Before examining the one-month sample, fix its UTC interval and budget. Keep the packet with the lowest Fink anomaly score per object in that interval, sort these packets by time and take 20 at evenly spaced indices. Lower scores are more anomalous. Time coverage follows available anomaly packets. The sample is purposive. Only history available through each selected packet enters its measurements. Current external catalog contents provide retrospective identification. [Fink anomaly API](https://doc.ztf.fink-broker.org/services/search/anomaly_detection/) · [Object histories](https://doc.ztf.fink-broker.org/services/search/objectid/).
4. Run the detector unchanged on the sampled targets. Exclude control IDs, their documented alias and matches to their recovered physical Gaia sources from independent exploration. A failed calibration closes the packet gate. All seven agreements are required for this pilot.
5. Apply the initial candidate filter. Triage checks include the foreground-star screen, per-packet CDS/TNS fields, a current SIMBAD cone of 3″, at least three quality-selected UTC dates, at least two dates with positive subtraction residuals, and anchor DRB ≥ 0.8. A history point requires DRB ≥ 0.8 and 0 < difference-magnitude uncertainty ≤ 0.2. Negative residuals can reflect a reference brighter than the science image. Missing TNS metadata is explicitly recorded. A failure of the broker's CDS service is recorded separately. [ZTF subtraction and filter fields](https://zwickytransientfacility.github.io/ztf-avro-alert/schema.html).
6. The first run exposed an important catalog gap. After independent review, a second-stage VSX query was added for **all** objects clearing the first filter. It archives each 5″ cone and checks the returned separations against its coordinates. One match within 1″ closes a new-source claim and makes the case eligible for an alias audit. Wider or multiple matches remain for inspection. Both the first-pass disposition and downstream catalog decision are preserved. Each packet has at most three cases. [AAVSO International Variable Star Index](https://vsx.aavso.org/) · [VSX catalog at VizieR](https://cdsarc.cds.unistra.fr/viz-bin/cat/B/vsx).

These seven literature cases test one contaminant screen. They are a retrospective positive/negative control check with known labels available to the developer. Sensitivity on a representative population, false discovery rate and readiness for astronomical confirmation remain unmeasured. A star track or near catalog neighbor suggests a counterpart. Registration, source confusion, reference subtraction and the quality of an underlying catalog classification can still affect the interpretation.

## Run another bounded window

For a fresh replication, register a discovery window and review budget before retrieving the new labels. This example tests the next disjoint calendar month.

```sh
python3 pilot.py fetch --start 2024-11-01 --stop 2024-12-01 --limit 20 --snapshot runs/2024-11
python3 pilot.py replay --snapshot runs/2024-11 --out runs/2024-11-results
python3 agents/review.py plan --packet runs/2024-11-results/review_packet.json
```

Fetch uses `curl`, sequential bounded HTTPS requests, per-request timeouts, limited retries, response-size caps and content checks. It supports 1–31 days and 1–20 objects. A partial first-stage download resumes when the same command is repeated. A completed primary snapshot is preserved. If the subsequent VSX download was interrupted, run `python3 pilot.py vsx --snapshot runs/2024-11`, then replay. Each raw response has its original request, UTC retrieval time and SHA-256 of the uncompressed bytes in its manifest. A frozen replay never silently refreshes a source. Historical alert selection and retrospective catalog snapshots can differ across refreshes. [Primary requests and checksums](data/pilot/manifest.json) · [VSX requests and checksums](data/pilot/vsx/manifest.json)

For the next moving-star case, compare the inferred crossing date with the coadd's earliest and latest exposure dates. Run `python3 motion.py OBJECT_ID --snapshot SNAPSHOT --out NEW_DIRECTORY --fetch` for the first capture and omit `--fetch` for offline replay. A control population of shifted sky coordinates and a second disjoint set of labeled stars would test chance alignments before scaling that screen. The number of closed false leads per expert review hour and the number of surviving cases with independent image evidence are the useful operating measures here.

This work uses public ZTF alert data served by Fink, ESA Gaia DR3 through NASA/IPAC IRSA, SIMBAD and VizieR at CDS, the AAVSO VSX and TNS. Individual rows retain their source identities. [Fink](https://fink-broker.org/) · [ZTF data access](https://irsa.ipac.caltech.edu/Missions/ztf.html) · [Gaia acknowledgements](https://www.cosmos.esa.int/web/gaia-users/credits) · [SIMBAD](https://simbad.cfa.harvard.edu/simbad/)

## Replication on October 2026

The next run used the last completed calendar month and a new set of 20 object IDs. The originally registered screen admitted four. All four had published VSX counterparts, so the sample yielded zero uncatalogued survivors. An exploratory check also found catalog counterparts for all five objects withheld after their historical broker crossmatches failed. [Registration](experiments/replication/protocol.md) · [Counts and object-level identifications](experiments/replication/2026-10/results/summary.md)

```sh
python3 coverage.py
python3 replication.py
make verify
```

These commands replay the checked-in archive on a clean checkout. The coverage script can retrieve missing source responses. The replication script reads the frozen cohort and its two additional VSX responses for historical service failures. The three alias-audit cases can be inspected with `python3 agents/review.py plan --limit 3 --packet experiments/replication/2026-10/results/review_packet.json`.
