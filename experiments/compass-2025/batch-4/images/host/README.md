# A faint galaxy at the position of ZTF25aanlbls

The central position of the transient lies 0.133″ from the Legacy Surveys DR10 galaxy at 250.02792336°, +29.76324886°. Fourteen positive ZTF alerts passing the compass quality cuts produced nine independent nightly coordinates from April 20 to May 29, 2025. Their median and 90th-percentile radial deviations from the nightly median center were 0.126″ and 0.187″. The measured separation is about one quarter of the galaxy’s 0.52″ effective radius. These uncertainties include nightly scatter, while differences between the ZTF and Legacy Surveys astrometric frames could also affect a displacement of this size. A nuclear transient and a projected off-center event both fit the localization. The next distinct DR10 optical source is 7.04″ away. [Numerical results](legacy.json) · [Positions and PS1 counterparts](results.json)

**Pre-event identification.** Pan-STARRS DR2 has one primary stack object within 10″, object 143712500279286504 at 0.121″. Its mean-table entry has zero single-exposure detections. A separate mean entry 2.08″ away has two detections and no primary stack source there. At the matched stack position, r PSF and Kron magnitudes are 22.99 ± 0.17 and 22.59 ± 0.19. In i they are 22.93 ± 0.19 and 22.81 ± 0.28. The g stack reports a noisy PSF measurement and lacks a valid Kron magnitude. The preselected nearby Gaia control stars have PSF-minus-Kron values near zero and are substantially brighter. PS1 alone gives weak shape information. [Control selection](protocol.md) · [MAST catalog definitions](https://catalogs.mast.stsci.edu/docs/panstarrs.html) · [STScI image service](https://ps1images.stsci.edu/ps1image.html)

**Resolved light.** The Legacy Surveys grz exposures span April 2015 to May 2018. DR10 fitted the counterpart as a round exponential galaxy with effective radius 0.522 ± 0.049″ and a 34.2 improvement in χ² over the point-source model. The fluxes correspond to g = 23.36, r = 22.50 and z = 22.18 AB. DR9 fitted essentially the same profile from the same observation interval. The DR10 central pixel has three exposures in each of g, r and z and zero in i. [Tractor table definitions](https://www.legacysurvey.org/dr10/catalogs/) · [Archived DR10 table subset and full-file digest](sources/legacy-dr10.json)

For an image-based check, we scaled the nearest Gaia field star to the counterpart’s 1.8″ aperture flux and injected it into 22 blank locations in each coadd. The table uses the ratio of 0.8″ to 1.8″ aperture flux. A smaller fraction puts more of the light outside the core. In the actual images the r and z counterpart is broader than the star at the matched flux. Three g controls have a concentration as low as the target’s. The same field star and sky positions contributed to all three filters. These descriptive ranks quantify local image noise and seeing. They do not measure the probability of a chance superposition. [Comparison images](profiles.png) · [Aperture calculation and control selection](legacy.py)

| Band | Counterpart | Stellar controls, median | Controls at or below target |
| --- | ---: | ---: | ---: |
| g | 0.467 | 0.542 | 3/22 |
| r | 0.470 | 0.596 | 0/22 |
| z | 0.493 | 0.737 | 0/22 |

**Approximate redshift.** DR10.1 gives a median photometric redshift of 0.534, with a catalog 68% interval of 0.440–0.623 and a 95% interval of 0.346–0.724. The DR9 viewer gives a mean of 0.519 with a distribution width of 0.104. These releases share imaging. The catalog authors identify increasing systematic errors for galaxies fainter than z-band magnitude 21. At z-band magnitude 22.18 this host estimate is a rough scale for testing classifications. DR10.1 has no reference spectroscopic redshift for the exact object. The public DESI DR1 spectrum position query and the SDSS DR18 primary-spectrum query returned zero matches within 5″. A rest-frame timescale or luminosity inferred using this redshift requires both a physical association and a reliable galaxy distance. [DR10.1 data model](https://www.legacysurvey.org/dr10/files/#photo-z-sweeps-10-1-photo-z-sweep-brickmin-brickmax-pz-fits) · [Zhou et al.](https://arxiv.org/abs/2309.06443) · [Redshift row and byte-range receipt](sources/photo-z-row.bin.json)

The Pan-STARRS1 and Legacy Surveys teams produced the pre-event images. The U.S. Department of Energy Office of Science, Office of High Energy Physics supported the PRLS photometric-redshift catalog through grant DE-SC0007914. [Survey acknowledgments](https://www.legacysurvey.org/acknowledgment/)

## Reproduce the measurements

From the repository root, use Python with NumPy and Astropy for the two catalog queries and the image-coordinate checks. The Legacy Surveys analysis also uses SciPy. Add `--plot` with Matplotlib installed to regenerate the figure. Each command below reads only the local archive. The retrieval receipts include the data request, UTC retrieval time and checksum. A new capture can reconstruct the compact Tractor subsets from the public release FITS. The photo-z extraction reads the original header and exact row by byte range and verifies its release, brick and object identity. [Source files and receipts](sources/)

```sh
python3 experiments/compass-2025/batch-4/images/host/check.py
python3 experiments/compass-2025/batch-4/images/host/legacy.py
```
