# Does the near-cutoff VSX match hold up?

ZTF18abdjpbo was the least certain positional match among the four sources admitted to VSX review. Its anchor lies 0.91″ from the listed variable, close to the registered 1″ threshold. The VSX object has the same ZTF identifier. A separate coordinate test strengthens the association while leaving its reported `AM+E:` classification open to physical verification. [VSX entry](https://vsx.aavso.org/index.php?view=detail.top&oid=2345383)

I queried Gaia DR3 via IRSA and the Pan-STARRS DR2 mean catalog via MAST in a 5″ cone around the alert. Each returned three entries. At the alert epoch, the position of Gaia DR3 4499667881166202624 propagated with its measured proper motion lies **0.186″ from the ZTF centroid**. The next Gaia source with measured motion is 4.083″ away. Pan-STARRS also records a source nearby, object 124102687712147056. [Archive requests, dates and checksums](sources.json)

The Gaia star moves south by 34.177 milliarcseconds per year. I repeated the comparison on the 269 UTC dates with usable positive ZTF residuals up to the selected alert, giving each date one median position. The median separation was **0.134″ from the epoch-adjusted Gaia position**, compared with **0.896″ from the fixed VSX position**. Propagating Gaia to J2000 gives a position 0.115″ from the listed VSX coordinate. This last comparison describes geometric consistency. The actual epoch used to determine the VSX position has not been established.

The repeated positions support the interpretation that ZTF recovered a known source. A positional match does not establish the cause of its variability. The calculation uses a small-angle approximation, neglects parallax and blends, and makes no claim about a formally measured proper motion in the ZTF images.

Run `python3 experiments/replication/2026-10/astrometry/check.py` from the repository root to reproduce the values offline.
