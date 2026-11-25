# Fresh-month replication

The preregistered UTC interval runs from 2026-10-01 through 2026-11-01 (exclusive). Fink lists 1,655,912 processed alerts across 23 observing-night rows. The anomaly archive contains 19 UTC dates, 190 packets and 130 distinct ZTF IDs. Both coverage thresholds passed. These two services count different products. [Archived coverage](../../coverage/coverage.json) · [Fink statistics](https://doc.ztf.fink-broker.org/services/search/statistics/)

The frozen sampling rule selected 20 IDs on 16 UTC dates. 0 IDs overlap the October 2024 pilot. The original 7/7 control outcomes reproduced on the current snapshot. The sample is new, while the seven controls are reused. [Registered protocol](../../protocol.md) · [Control and object measurements](report.md)

| First-pass disposition | Objects |
|---|---:|
| broker crossmatch failed | 5 |
| catalog context | 11 |
| review | 4 |

Among 4 first-pass objects, the 5″ VSX check closed 4 as unique variable-star associations within 1″. **0 remained after the frozen rules.** [Machine-readable stage counts](summary.json) · [Bounded review packet](review_packet.json)

| First-pass case | Anchor residual | Usable UTC dates | VSX | Downstream disposition |
|---|---|---:|---|---|
| [ZTF18abjedws](https://ztf.fink-portal.org/ZTF18abjedws) | positive | 94 | [MGAB-V713](https://vsx.aavso.org/index.php?view=detail.top&oid=1499062) (UGZ, 0.15″) | closed as a new source, retained for alias review |
| [ZTF18abdjpbo](https://ztf.fink-portal.org/ZTF18abdjpbo) | positive | 292 | [ZTF18abdjpbo](https://vsx.aavso.org/index.php?view=detail.top&oid=2345383) (AM+E:, 0.91″) | closed as a new source, retained for alias review |
| [ZTF18abkwscm](https://ztf.fink-portal.org/ZTF18abkwscm) | negative | 195 | [ZTF18abkwscm](https://vsx.aavso.org/index.php?view=detail.top&oid=2344838) (UG, 0.08″) | closed as a new source, retained for alias review |
| [ZTF18abmfrmh](https://ztf.fink-portal.org/ZTF18abmfrmh) | negative | 558 | [ASASSN-20en](https://vsx.aavso.org/index.php?view=detail.top&oid=1543944) (UG, 0.08″) | closed as a new source, retained for alias review |

## Broker service failures

A separate coverage audit examined the archived SIMBAD cone for every historical broker crossmatch failure. Entries with no SIMBAD neighbor were queried against VSX. This exploratory audit leaves the frozen sample dispositions above intact.

| Object | Broker status | Independent lookup |
|---|---|---|
| [ZTF18aaaduso](https://ztf.fink-portal.org/ZTF18aaaduso) | Fail 503 | VSX · [MGAB-V258](https://vsx.aavso.org/index.php?view=detail.top&oid=838077) (UGSU, 0.31″) |
| [ZTF17aaeencn](https://ztf.fink-portal.org/ZTF17aaeencn) | Fail 504 | VSX · [MGAB-V3757](https://vsx.aavso.org/index.php?view=detail.top&oid=2215141) (UG, 0.23″) |
| [ZTF18acxgxrh](https://ztf.fink-portal.org/ZTF18acxgxrh) | Fail 503 | SIMBAD · MASTER OT J205110.36+044842.2 (CV*, 0.21″) |
| [ZTF17aaanyrs](https://ztf.fink-portal.org/ZTF17aaanyrs) | Fail | SIMBAD · EM* LkHA  290 (Y*O, 0.23″) |
| [ZTF18aabszgd](https://ztf.fink-portal.org/ZTF18aabszgd) | Fail | SIMBAD · NVSS J104939+154838 (BLL, 0.30″) |

The archive returned 11,230 history rows across controls and exploration objects. 438 contained broker metadata with null object ID, time, position, candidate ID and all photometric fields. Those rows remain in the raw archive and contribute zero observations to the analysis. No alert-selection or classification threshold changed.

Source ID overlap alone cannot exclude all astrophysical aliases. Weather and observing cadence affect coverage. Broker crossmatch failures are listed separately and may hide additional sources. These counts measure a fixed 20-object sample, so multiple unreviewed objects and observing dates remain in the month.
