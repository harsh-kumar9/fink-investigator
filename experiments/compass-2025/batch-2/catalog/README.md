# Earlier identifications and observer coverage

All three positions have a named counterpart. The SIMBAD cone returned one source within 10″ of each anchor. VSX, ATLAS, ASAS-SN and Gaia DR3 variability catalogs were queried within 5″. The table gives the static SIMBAD separation and the nearest two Gaia separations after propagation to the alert epoch. [Catalog rows and query receipts](results.json) · [Gaia measurements by epoch](../images/results.json)

| Fink position | Prior identification | SIMBAD distance | First and second Gaia source at the anchor |
|---|---|---:|---:|
| ZTF25aacfjde | TCP J07222683+6220548, also ASASSN-25aj and GALEX J072227.0+622057. VSX identifies an outbursting helium binary. | 0.233″ | 0.065″, 33.427″ |
| ZTF18aagstwh | LINEAR 15938573, also ZTF J155842.59+511549.0. VSX gives an RRab fundamental pulsator at 0.62232 d. ASAS-SN gives 0.62233 d. ATLAS's Fourier-pipeline Lomb–Scargle estimate is 0.622518 d and its other estimator gives 1.244667 d. Gaia classifies it RRab. | 0.164″ | 0.141″, 7.141″ |
| ZTF18abmfrjj | V1153 Aql, a Mira. VSX gives 325 d and ASAS-SN gives 328.22 d. The ATLAS classification is irregular with discordant estimates of 0.993990 and 443.379574 d. | 0.153″ | 0.024″, 4.245″ |

Gaia's independent RRab fit gives a fundamental period of 0.62231877 ± 0.00000148 d. The first-overtone field is empty. For V1153 Aql, the Gaia long-period-variable frequency corresponds to 328.7 ± 8.9 d using first-order uncertainty propagation. [Gaia RR Lyrae query](ZTF18aagstwh-gaia.xml.json) · [Gaia long-period-variable query](ZTF18abmfrjj-gaia.xml.json)

The Mira has a Gaia neighbor inside 5″, relevant to measurements when the Mira becomes faint. The ID and the catalog periods establish earlier reports for each periodic source. [V1153 Aql in VSX](https://vsx.aavso.org/index.php?view=detail.top&oid=2292) · [LINEAR 15938573 in VSX](https://vsx.aavso.org/index.php?view=detail.top&oid=321240)

## What observers already knew about the helium binary

**January.** NMW measured the source at 12.8 mag on January 20.9416 UT. Its discovery notice records an ASAS-SN limit of g > 17 on January 20.3957. Yoshimoto's January 21 follow-up placed the source at the Gaia counterpart. The Fink team's February 28 ATel describes brightening on January 21 and subsequently fading. Neither notice singles out January 19. The VSOLJ search for January 17–21 found 785 observations, all on January 21. This leaves the timing and amplitude of any earlier excess to the ZTF measurements. [Original CBAT notice](http://www.cbat.eps.harvard.edu/unconf/followups/J07222683+6220548.html) · [Lima et al., ATel 17058](https://www.astronomerstelegram.org/?read=17058) · [VSOLJ coverage receipt](vsolj-january.json)

**February and March.** Tampo and colleagues announced the AM CVn interpretation on February 28. Tarasenkov and colleagues subsequently published the discovery and two February spectra with helium absorption. They examined archival data and reported no obvious pre-outburst rise. They measured a 46.87 ± 0.12 minute superhump during the first rebrightening and identified shorter rebrightenings at roughly days 30, 39 and 55 after discovery. Their data include ZTF through Lasair as well as independent telescopes. The middle panel of their Figure 1 shows the outburst through roughly day 60. A later INASAN report, section 3, describes two months of monitoring and illustrates the main outburst and first rebrightening. [VSNET alert 28049](http://ooruri.kusastro.kyoto-u.ac.jp/mailarchive/vsnet-alert/28049) · [Tarasenkov et al., sections 2–3](https://arxiv.org/html/2505.20842v1) · [2026 INASAN report](https://doi.org/10.51194/INASAN.2026.11.2.010)

**April.** The independent observer archive already contains a bright run on April 15. VSOLJ returned 195 observations for the month, including 180 on April 15 attributed to H. Itoh. The April 14 and 16 entries are from the same observer. The April 30 VSNET CV Circular lists those flanking dates and flags the April 16 value as uncertain. Its summary skips April 15. Itoh was also a coauthor on the discovery paper. The archive establishes prior observational coverage of the episode. Its measurement dates do not establish when every observation became public, or whether the episode has previously been analyzed in another publication. The ZTF comparison in this batch gives a separate two-band measurement on April 16. [April observing inventory](vsolj-april.json) · [VSOLJ observer codes](http://vsolj.cetus-net.org/ObserversCode-latest.txt) · [Contemporaneous April circular](http://ooruri.kusastro.kyoto-u.ac.jp/mailarchive/vsnet-cvcirc/8929)

Three AAVSO queries covering April returned zero measurements, using the VSX name, ASASSN-25aj and AUID 000-BQC-247. The same service returned 192 measurements on four dates using the VSX name for January 25–February 19. An empty search of AAVSO alone would have missed the independent April coverage. [AAVSO query](april-name.aid.json) · [VSX identity and AUID](vsx.xml.json) · [AAVSO positive-control query](active.aid.json)

## Replay

From this directory, compare the catalog summary and inspect the observing coverage.

```sh
python3 check.py | diff - results.json
python3 observer_priority.py
```

The two VSOLJ receipts record per-date observer counts, the query and a checksum of the full response. The public archive identifies its raw observations as unvalidated and asks researchers to contact VSOLJ before using them. Only the coverage inventory is stored here. To verify those two receipts against the current service, run `python3 observer_priority.py --verify-live`. [VSOLJ archive](http://www.cetus-net.org/cgi-bin/obs_search.cgi)
