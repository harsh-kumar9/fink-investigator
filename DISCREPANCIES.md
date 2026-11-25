# How we record disagreements

An alert score ranks something unusual. Deciding what changed in the sky requires a second measurement capable of contradicting the first interpretation. For every discrepancy worth follow-up, record the following fields before assigning a physical explanation.

| Field | What to preserve |
|---|---|
| Question and unit | The object or exposure, case ID, tested claim, selection rule and comparison epoch. Group repeated detections from the same source. |
| Method outputs | Each method's observed value or missing field, with units, time, sky position or search radius, original URL and frozen response when possible. Distinguish a successful empty query from service failure. |
| Competing interpretations | At least one conventional physical or instrumental explanation that could produce the observations. Record sources of shared error across products. |
| Discriminating evidence | A predicted difference between the explanations and the control that checks the test. Preserve outcomes that fail the threshold or abstain. |
| Resolution status | State exactly which question closed. A catalog alias can resolve prior identification while the physical mechanism remains open. |
| Next discriminating observation | Specify the measurement that could still reverse the conclusion, its source and the scale, sign or timing to compare. |

## October examples

* **The source was already reported.** A missing TNS field and empty stellar-catalog cones left ZTF26acbcwpz open. Caltech explicitly links its ZTF identifier to SN 2026afgd, whose TNS record includes a P60/SEDM spectrum. The identification closes the new-source reading. [Archived registry evidence](experiments/october-census/sources/registry/sn-2026afgd.json) · [TNS record](https://www.wis-tns.org/object/2026afgd)
* **The star changed position.** A fixed coordinate search could miss G 143-8 at the ZTF22abgqdcj residual. Its independently measured motion predicts the date it crossed the alert position inside the interval used for the reference image. Repeated paired residuals and their signs support the moving-star interpretation. A registered multi-epoch fit could test a superposed source. [Track and reference dates](experiments/october-census/review/ZTF22abgqdcj/reference/report.md) · [Images](experiments/october-census/review/ZTF22abgqdcj/images/gallery.md)
* **The catalog coordinate was close to the cutoff.** For ZTF18abdjpbo the fixed VSX offset was 0.91″. On 269 dates, allowing for the Gaia star's measured motion reduced the median separation to 0.134″, versus 0.896″ to the fixed VSX coordinate. Identity is supported, while the published subtype remains uncertain. [Astrometric check](experiments/replication/2026-10/astrometry/README.md)
* **The star has a position but its variability is unclassified in the queried tables.** The direct-image light curve for ZTF17aacydyv agrees with its difference-image history. The data leave open a variable at the Gaia position, an unresolved neighbor and correlated systematics. Sharper observations can test which star varies. [Independent review](experiments/october-census/review/ZTF17aacydyv/independent-review.md)

The complete October record preserves initial decisions for all 130 IDs and parallel downstream checks, including historical broker outages. The machine ledger gives source-linked observations, alternatives, dispositions and next tests for all seven cases advancing beyond the secondary catalogs and for the moving-star case. The near-cutoff alias is recorded separately above. [Month reconciliation](experiments/october-census/README.md) · [Per-object outcomes](experiments/october-census/results/outcomes.csv) · [Case ledger](experiments/october-census/results/discrepancies.json)

A parallel biodiversity census keeps the same distinction between a source record, competing interpretations, the adjudication and a decisive next observation. [Marine Ledger discrepancy cases](https://github.com/harsh-kumar9/marine-ledger/blob/main/experiments/full-page-transfer/DISAGREEMENTS.md)
