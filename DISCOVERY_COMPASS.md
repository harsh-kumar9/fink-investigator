# Fink discovery compass

Start with a question whose answer could change after a second measurement. Fink's anomaly score supplies a small, unusual set of alerts. The compass samples that set across the calendar, uses the observations available at each alert to rank objects, and saves what happened to each one. The first wide run covered 974 ZTF identifiers and recovered a supernova, a dwarf-nova outburst and a position shift amplified by image subtraction. [2025 findings](experiments/compass-2025/README.md) · [Ranked queue](experiments/compass-2025/results/report.md) · [Fink's nightly selection](https://doc.ztf.fink-broker.org/services/search/anomaly_detection/)

## Choose the next measurement

| What drew attention | Question to ask | Measurement that can change the interpretation |
|---|---|---|
| A short history with positive detections in two filters | Did something new brighten here? | Locate the residual relative to the host or a foreground star. Search for a dated discovery and spectrum. Compare with pre-event photometry or upper limits. SN 2025oq had a spectrum and no matching DR24 static-source light curve. |
| Bright residuals on widely separated dates | Were there separate outbursts? | Measure the source directly inside both groups and between them. Count distinct episodes when there is an observed return to the usual level. At FT Cam, one of the two groups was fainter than the ordinary red level. |
| Repeated positive and negative residuals | How much did the source fade, and is it predictable? | Compare direct magnitudes by sign in each filter. Predict held-out seasons or dates with another filter. Try the double period when deep and shallow minima alternate. [11.62-day example](experiments/october-census/timing/README.md) |
| Bright and faint residuals land at different positions | Which source varied, and did the reduction move its centroid? | Learn a direction on early opposite-sign pairs and test later pairs. Propagate Gaia coordinates, separate detector/reference groups, and measure direct positions on the same exposures. [Worked diagnostic](experiments/compass-2025/README.md#where-did-the-position-shift-come-from) |
| A catalog service failed or a cataloged star looks familiar | Could this position still show unmeasured behavior? | Repeat the missing lookup and preserve its status. Investigate the behavior separately from the name. Matched field stars helped establish a later bright state in ZTF18abqclig. [Service-failure and field-star example](experiments/november-search/README.md) |

For a faint source present in the reference image, compare the signed difference flux, the reference-plus-difference total and the direct source on the very same exposure. A difference magnitude near the quiescent magnitude can mean the star nearly doubled in total flux. A published plot can already contain that measurement. Check the figure, measured baseline and earlier observer records before assigning novelty. [J0722 worked example](experiments/compass-2025/batch-2/README.md)

The first three rows have label-blind routing rules. The fourth has a separate centroid search. The last is a recovery procedure for candidates encountered anywhere in the pipeline. A route score sets review order. It is neither a probability of discovery nor an astrophysical classification.

## Resume the 2025 queue

The fifth round used the final eligible recent-onset source in the 2025 sample. **53 unique positions remain** on the recurrent and fading routes. The independent 2024 search added a fresh year with ten eligible positions left after its first packet. The two commands below advance the respective queues and preserve their original ranks. [Fifth-round findings](experiments/compass-2025/batch-5/README.md) · [2025 queue](experiments/compass-2025/batch-5/results/report.md) · [2024 findings and queue](experiments/compass-2024/README.md)

```sh
base=experiments/compass-2025
python3 compass.py --snapshot "$base/snapshot" --out runs/compass-2025-batch-6 \
  --reviewed-packet "$base/results/packet.json" --reviewed-packet "$base/astrometry/packet.json" \
  --reviewed-packet "$base/batch-2/results/packet.json" --reviewed-packet "$base/batch-3/results/packet.json" \
  --reviewed-packet "$base/batch-4/results/packet.json" --reviewed-packet "$base/batch-5/results/packet.json"

other=experiments/compass-2024
python3 compass.py --snapshot "$other/snapshot" --out runs/compass-2024-batch-2 \
  --exclude-ledger "$base/results/ledger.json" --reviewed-packet "$other/results/packet.json"
```

Three checks from the fifth round travel to other sources.

* **Separate a peak sequence from a lasting clock.** Measure intervening faint returns and compare peak strengths in the same filter. Try the same rule in independent seasons. For a coarse-pixel telescope, localize the variable signal on difference images and shift the time series as a control. [MGAB-V960 with ZTF and TESS](experiments/compass-2025/batch-5/variables/README.md)
* **Resolve the host before interpreting a faint pre-event source.** Center the event on a registered deep image, propagate reference-star positions to the image epoch and measure offset positions on equally bright host structure. State the aperture and local depth. [Subaru, DECam and Hubble at SN 2025xza](experiments/compass-2025/batch-5/images/host/hst/README.md)
* **Count how many standards and comparison objects actually survive.** A sparse field can require a wider stellar calibration on the same detector. A bright comparison panel can still have zero usable late baselines. Preserve those counts when reporting the result of a relaxed check. [Original image measurements and Type II panel](experiments/compass-2025/batch-5/README.md#four-months-of-red-light-from-sn-2025xza)

Two tests from the fourth round apply to other sources.

* **Test a historical classification at its own epochs.** Retrieve the original model, amplitudes and dates. A measured brightness change larger than the model’s full range can challenge it without knowing the absolute phase. Inspect same-exposure field stars and flagged pixels, then transfer the model to an independent telescope. [The six-day rise at MGAB-V488](experiments/compass-2025/batch-4/catalog/README.md)
* **Separate an observed light-curve resemblance from a physical identification.** Measure red fading and color evolution on matched dates, test known spectroscopic labels, and perturb the target fluxes to see how often the ranking changes. Check the host distance: one similar-looking light curve can imply very different luminosities. [AT 2025iac comparison](experiments/compass-2025/batch-4/comparison-protocol.md)

Earlier rounds added these checks.

* **Test eruption spacing across surveys.** Measure onset gaps where the sampling resolves both eruptions and an intervening fade. Freeze a range of delays. On unused dates from another telescope, compare brightness changes over that delay with those at shorter and longer delays starting on the same dates. Check nearby stars on the same images and time-dependent calibration. This can recover characteristic spacing even if a single clock fails on later seasons. [MGAB-V243 test and controls](experiments/compass-2025/batch-3/variables/recurrence-protocol.md)
* **Measure the quiet nights.** The broker can stop returning observations while the survey continues imaging the position. Use the exposure inventory, freeze image quality and selection criteria, and measure signed flux at the same position. Withhold a nearby star of known brightness to check sensitivity. For a proposed outburst duration, count how much of the observing window could have caught it. [ASASSN-25ba image protocol](experiments/compass-2025/batch-3/onset-protocol.md)
* **Challenge a single archival flare on its own exposure.** Ten nearby stars showed the same apparent excess as the candidate on a Pan-STARRS image. Inspect shared flags, local sky and comparison-star flux before assigning an earlier eruption. [Pan-STARRS field comparison](experiments/compass-2025/batch-3/images/ps1/README.md)

## Run another year

Run these commands from the repository root with Python 3.10 or newer and curl. The queue and catalog replay use the standard library. The deeper image and time-series checks list their scientific Python dependencies in the batch reports. The optional figures use Matplotlib. The example below starts a new 2023 search. The checked-in cohorts cover 2025 and 2024.

```sh
make verify
make verify-compass

run=runs/compass-2023
python3 compass_archive.py fetch --year 2023 --snapshot "$run/snapshot"
python3 compass.py --snapshot "$run/snapshot" --out "$run/results" \
  --exclude-ledger experiments/compass-2025/results/ledger.json \
  --exclude-ledger experiments/compass-2024/results/ledger.json
python3 agents/review.py plan --packet "$run/results/packet.json" --limit 3
python3 compass_direct.py --snapshot "$run/snapshot" --ledger "$run/results/ledger.json" \
  --sources "$run/direct" --out "$run/results/direct.json" --fetch
python3 compass_astrometry.py --snapshot "$run/snapshot" --ledger "$run/results/ledger.json" \
  --out "$run/astrometry"
python3 compass_reference.py --snapshot "$run/snapshot" --packet "$run/astrometry/packet.json" \
  --sources "$run/reference-sources" --out "$run/results/reference.json" --fetch
```

1. **Register the window and budget before reading labels.** Record the calendar year, quota and any changes to the questions. Default sampling takes up to four young, four other positive and four remaining anchors per month, then fills unused places. Use `--per-stratum 1` for a budget of three histories per month. An ID belongs to the month when it first appears. A month with no returned anomalies contributes zero; it says little about the sky or the observing schedule. [Exact sampling, eligibility and tie rules](experiments/compass-2025/protocol.md)
2. **Read the gate and the ledger before following a case.** The original seven moving-star/supernova controls test foreground contamination. Separate controls test each photometric rule. The 2025 record preserves a failed first control panel and its coverage-only amendment. The astrometric geometry tests use injected and null shifts; real false-positive rates are unknown. Previously studied pilot, October and November positions are excluded automatically. Pass `--exclude-ledger` for each previous wide run whose sampled positions should be excluded. The 5″ rule flags possible aliases; confirm names and propagated positions. [Control outcomes](experiments/compass-2025/calibration/README.md)
3. **Give the same packet to an investigator and a skeptic.** Record independent first opinions, the strongest alternatives and their predicted observational differences. The primary queue is in `results/report.md`, with the full rank and decision for each ID in `results/ledger.json`. The separate astrometric packet contains the actual paired dates and candids. When taking another backlog case, save its original rank and keep the earlier decisions. [Review roles](agents/README.md) · [Example first opinions](experiments/compass-2025/investigator.md) · [Skeptic's first opinions](experiments/compass-2025/skeptic.md)
4. **Record what would reverse the result.** Save the actual exposure or registry record, the cutoff date, how dates were paired, excluded or missing measurements and whether the observations share a telescope or reduction. Distinguish a documented source, a newly measured behavior and a confirmation of its cause. Put the open test beside the finding. [Discrepancy record](DISCREPANCIES.md)

To replay a frozen archive, use `compass_archive.py replay --snapshot PATH` and omit `--fetch` from the two photometry/reference commands. Archive capture resumes interrupted downloads and validates source requests, identities and checksums. Use a fresh snapshot directory when changing selection settings. The discovery sample is drawn from the broker's nightly anomaly selections; quiet fields, single-filter transients, sparse events, long gaps and objects missing from the broker remain poorly represented.
