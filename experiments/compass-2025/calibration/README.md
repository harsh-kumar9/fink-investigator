# What the three routes recovered on real sources

The wide search uses small control sets with known behavior. One source can test several routes. The effective panel contains four distinct physical objects. The moving-star test uses seven additional objects.

| Route | True positives | True negatives | False positives | False negatives | Eligible positive controls |
|---|---:|---:|---:|---:|---|
| Recent onset | 1 | 3 | 0 | 0 | SN 2022oqm |
| Repeated bright detections | 2 | 2 | 0 | 0 | ASASSN-20en, CRTS J222643.4-200536 |
| Repeated fading | 1 | 2 | 0 | 0 | CRTS J222643.4-200536 |
| Original Gaia moving-star screen | 4 | 3 | 0 | 0 | Four literature Gaia IDs recovered |

SN 2022oqm was first reported in July 2022. Through its August 7 control packet, Fink recorded positive measurements on 11 usable nights, including five with both filters. ASASSN-20en had pairs of green positive nights in five distinct years from 2020 through 2024. The CRTS RR Lyrae had repeated positive and negative green detections from 2019 through 2024. G 143-8, a moving star that leaves a negative residual at its reference position, passed none of the three morphology routes. The dwarf nova also triggered the fading route. Its fading outcome was unregistered and contributes to neither side of that route's counts. [Supernova discovery observations](https://www.wis-tns.org/index.php/astronotes/astronote/2022-142) · [Control labels and cutoff dates](expectations.json) · [Results and example dates](results.json)

**Why two panels are recorded.** V1 could evaluate three of its four registered sources. Its dwarf-nova control MGAB-V713 had 153 recorded Fink detections, all from 2026. With no observations at or before the 2024 cutoff, the V1 onset and recurrent gates stayed closed. The registered V2 amendment substituted the first externally cataloged dwarf nova in a fixed list with any timestamp in the same 2024 interval. ASASSN-20en had 95 such alerts. V2 retained the original anchor rule, expected outcomes and thresholds; the second name in the list never entered the test. The earlier result remains in the machine record. [Original test definition](protocol.md) · [Coverage amendment](replacement/protocol.md) · [Replacement control identities](replacement/expectations.json)

From the repository root, `python3 compass_calibration.py` reproduces the result offline from checksummed public histories. Its source records include each ZTF packet's contemporaneous first-detection field and exclude broker classifications. The baseline correct-identity test replays the original Gaia archive. Calibration subjects and the documented G 260-1 alert alias are excluded from candidate slots. These small, deliberately chosen historical examples establish that the rules run as intended on the named cases. The missed-source rate for other objects remains unmeasured. [Calibration code](../../../compass_calibration.py) · [First Fink receipt](manifest.json) · [Replacement receipt](replacement/manifest.json) · [Original motion controls](../../../data/controls.json)
