# A known ATLAS variable brightened in two ZTF bands

The November search recovered ZTF18abqclig at a position the ATLAS survey had flagged for stochastic variability in 2018. The ZTF source was **0.903 magnitude brighter in red in 2023–25 than in 2021–22**. The red data suggested the date intervals. We then registered a test in green on separate exposures, correcting each exposure against nearby stars. The change measured **0.926 magnitude in green**, about 2.35 times the earlier median flux. A second ATLAS candidate lies near ZTF20aaothpc, where five matched dates leave the photometric check inconclusive. [Green test registered before retrieval](brightening-protocol.md) · [Measurements](results/brightening.json) · [ATLAS catalog paper](https://arxiv.org/abs/1804.02132)

![Annual median brightness of ZTF18abqclig in red and field-corrected green, relative to each filter's 2021 median. The 2023, 2024 and 2025 medians in both bands sit about one magnitude above the 2021 and 2022 medians. Vertical bars enclose the middle half of measured nights.](results/brightening.png)

The green comparison covers 87 dates in 2021–22 and 94 dates in 2023–25. The five nearby stars were selected by their location, mean magnitude and exposure count. Their median shifts were all smaller than 0.05 mag in absolute value. Giving each of the eight shared calendar months equal weight yields a target brightening of 0.987 mag. Seven dates had enough measured comparison stars and lacked a usable target measurement. Even assigning all six missing early dates arbitrarily bright and the missing late date arbitrarily faint leaves a median shift of at least 0.909 mag across those field dates. The plot uses the 2021 median of each band as its zero. Bars show the middle half of the nightly measurements. Seasonal gaps limit when the change can be dated, and a Gaia neighbor 2.118″ away leaves the varying component uncertain. [Field and light-curve receipts](sources/brightening/) · [Graphic source](figure.py) · [IRSA direct-photometry archive](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html)

## What the November screen returned

The archive returned 120 packets. The registered UTC interval of November 1 through November 24, 2026, contained 110 packets on 11 dates from 86 distinct ZTF IDs. We reviewed the fixed, time-spread sample of 20 IDs. Nineteen were new to both October cohorts. ZTF18aabfdwa repeated its October 2026 ID and matched the same SIMBAD source. The ID, catalog and coordinate checks found no other overlap. The seven reused controls again gave four correct moving-star identifications and three supernovae with no moving-star flag. [Selection protocol](protocol.md) · [Archive requests and hashes](snapshot/manifest.json) · [All 20 outcomes](results/audit.json)

Every one of the 110 packets carried the Fink CDS label `Fail`. All 20 sampled objects stopped at that gate. This matters for a discovery search because an unavailable catalog lookup can hide both familiar sources and genuine variability. In the separate outage audit, a new SIMBAD lookup found neighbors for 15 positions. Of the other five, three had a single VSX variable within 1″. The remaining two passed the original photometric cuts. Each had an empty 5″ cone in the Gaia DR3 variability-classifier table and both Chen et al. ZTF tables. An exploratory lookup then found one entry within 1″ of each position in the published ATLAS variable-candidate catalog. Neither ASAS-SN cone returned an entry. All 20 sampled positions have at least one catalog neighbor under these checks. [Fink archive and failure labels](https://doc.ztf.fink-broker.org/services/search/anomaly_detection/) · [Failure-audit procedure](failure-audit.md) · [Gaia and Chen matches](results/extended.json)

| Audited position | Published ATLAS candidate | Direct red-band result | Nearby Gaia DR3 sources at alert epoch |
|---|---|---|---|
| ZTF18abqclig | ATO J092.6615+20.6346, 0.19″, `STOCH` | 6 positive and 169 negative dates, rank correlation 0.998, median 0.886 mag brighter on positive dates | 3375210671016297472 at 0.165″ and 3375210671016297728 at 2.118″ |
| ZTF20aaothpc | ATO J130.6711−02.3818, 0.51″, `IRR` | 5 positive dates, zero negative dates | 3072346890316182144 at 0.472″ and 3072346890318619520 at 1.670″ |

For ZTF18abqclig, images in 2019, 2022 and 2026 show a compact residual at the central source. The telescope recorded a negative sign in all three exposures. The independently stretched displays show the location. The calibrated photometry supplies the sign and brightness. The ATLAS catalog gives both the difference-image and proximity flags a value of 1, supporting its variability detection. Fink and the direct light curve share ZTF exposures, while ATLAS observed with separate telescopes. A fainter Gaia neighbor 2.118″ away could still contribute to the ZTF signal. The ATLAS stochastic label describes the light curve and leaves the physical mechanism open. [Three image triplets](review/ZTF18abqclig/gallery.md) · [ATLAS catalog definitions](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/156/241)

At ZTF20aaothpc, the three sampled science and template images show a central source. The difference image shows a compact residual near it. On some dates the residual core appears uneven. All ten quality-selected alert dates are positive. The direct comparison recovered only five red-band dates, below the 20-date minimum, and the ATLAS difference-image flag is 0. A two-source fit to ZTF images spanning a range of seeing, with positions fixed to the two Gaia trajectories, could distinguish changes in one star from seeing-dependent blending. The same test can localize ZTF18abqclig. [Three image triplets](review/ZTF20aaothpc/gallery.md) · [ATLAS data release](https://archive.stsci.edu/hlsp/atlas-var) · [Original photometry control and thresholds](../october-census/photometry-protocol.md)

## Replay

Independent first readings identified blends, reference-image effects and previously recorded variability as competing explanations. A subsequent ALeRCE lookup returned a tentative young-stellar-object classification for ZTF18abqclig. Its classifiers disagree on subtype. The search also records which public catalog requests returned access errors. [Investigator review](investigator.md) · [Skeptic review](skeptic.md) · [Prior-identification checks](prior/README.md)

Run these commands from the repository root. They read the archived responses offline and regenerate the measurements.

```sh
python3 pilot.py replay --snapshot experiments/november-search/snapshot --out experiments/november-search/results
python3 experiments/november-search/audit.py
python3 experiments/november-search/followup.py
python3 experiments/november-search/brightening.py
```

With Matplotlib installed, run `python3 experiments/november-search/figure.py` to rebuild the graphic.
