# Does every 5.8-day interval contain the same fade?

At ZTF18abjtfpx, green-band observations distinguish the two halves of the proposed 11.62356-day cycle. The red training fit fixed the cycle, the deep-fade center at HJD 2458210.38245 and the window width of 0.92988 days. We placed a second window of identical width exactly half a cycle later. No parameter was fitted in green. Each UTC date contributed at most one measurement per source and band. [Red fit and its selection history](../dips-protocol.md) · [Source photometry](../green/README.md)

| Position and data selection | Deep window median, mag (dates) | Half-cycle median, mag (dates) | Remaining median, mag (dates) |
|---|---:|---:|---:|
| Target green, all | 18.814 (27) | 17.361 (27) | 17.261 (340) |
| Target green, dates absent from both red fields | 18.910 (6) | 17.348 (11) | 17.247 (71) |
| Comparison A green, dates absent from both red fields | 16.870 (6) | 16.869 (9) | 16.874 (67) |
| Comparison B green, dates absent from both red fields | 17.857 (6) | 17.872 (11) | 17.862 (73) |
| Target red, original field's later dates | 17.198 (7) | 16.559 (7) | 16.403 (139) |
| Target red, second field's dates absent from original | 17.489 (3) | 16.583 (5) | 16.405 (21) |

The target fades deeply once per 11.62-day cycle. Half a cycle later its median is only 0.10 mag fainter than the remaining green dates. Both red fields also have a smaller dip near that phase. Unequal eclipses offer a physical explanation for the 5.8-day peak and the alternating depths. The data support an eclipsing-binary candidate at the position. The blending and sparse phase coverage leave the responsible star, geometry and orbital period open to independent confirmation.

ATLAS published 5.796240 days from one periodogram, 5.793208 days from its long-period Fourier fit, and 2.408143 days from a separate periodogram. The short-period fit indicator was zero, and its short-fit period fields were empty. ATLAS assigned zero to the four contact/detached binary scores and labeled the case `dubious`. The archived row contains no 11.62-day period in its period fields. [Primary ATLAS record](https://vizier.cds.unistra.fr/viz-bin/VizieR-5?-source=J/AJ/156/241/table4&recno=2182150) · [Original source and blend check](../catalog/README.md)

Run `python3 compare.py` in this directory with the same packages as the timing analysis to regenerate all medians and catalog fields. [Full counts and medians](results.json)
