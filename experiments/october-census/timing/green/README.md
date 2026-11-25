# The green-band check on ZTF18abjtfpx

The 5.81178-day red-band waveform predicted part of the later green-band variation. On the 30 held-out green nights with no measurement in the selected quality-cut red series, the rank correlation was **0.561** and the squared prediction error fell **21.8%** relative to the training mean. Both fell short of the advance thresholds of 0.60 and 30%. The case retains its status as a red-band timing lead. [Frozen green test and initial investigator assessment](protocol.md) · [Red-band timing results](../results.json)

Two IRSA DR24 green series lay within 1.5″. The series with the most catalog observations was 0.389″ from the Fink anchor. We kept one measurement per date. The red model fixed the period, phase, amplitude and waveform. The only parameter estimated from the first 60 green-only dates was an additive offset of 0.881 mag. We used that offset to forecast the other 30. For the secondary analysis including dates observed in both colors, the corresponding 132 held-out dates gave a correlation of 0.681 and a 29.0% error reduction. Across all 90 green-only dates the correlation was 0.601. [IRSA light-curve interface](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) · [Results and exact model coefficients](results.json)

| Green series | Distance from anchor | Green-only dates, training / held out | Held-out rank correlation | All green-only dates, rank correlation |
|---|---:|---:|---:|---:|
| Target 538101200105064 | 0.39″ | 60 / 30 | 0.561 | 0.601 |
| Field star 538101200000658 | 28.55″ | 56 / 28 | 0.057 | −0.007 |
| Field star 538101200001178 | 40.98″ | 61 / 31 | 0.248 | 0.146 |

The field stars share the target's detector quadrant and are within 0.5 mag of its mean green brightness. The first field star has fewer than the 30 held-out dates required for a threshold decision. The larger correlation at the target lends some support to a change within its aperture. The green and red exposures share ZTF's observing window, and the crowded field leaves the variable component unresolved. [Nearby Gaia stars](../../results/extended.json) · [Field selection and archived request](field.json) · [Green measurement request and checksum](green.json)

Run `python check.py` in this directory with NumPy, SciPy and Astropy to regenerate the JSON from the two saved IRSA responses and the parent red series. The analysis checks the request URLs and response hashes before calculating predictions.
