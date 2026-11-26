# Hubble’s 2022 view of the transient site

HST program 16691 observed the future position of ZTF25absffwy on August 30, 2022. WFC3/UVIS accumulated 710 seconds in F555W and 780 seconds in F814W. The images show several faint sources on structured galaxy light near the later transient. A 0.15″ aperture placed at the registered ZTF position gives a modest positive signal in both filters. These data cannot assign the ground-based peak to a particular star. [Two-filter figure](field.png) · [Original observation footprints](sources/coverage.json)

**Registration.** The visit mosaic used the GSC242 solution with a reported fit RMS of about 0.12″ on each coordinate. Four unsaturated Gaia DR3 stars with measured proper motions and RUWE ≤ 1.4 had matches in the HST visit catalog. Propagating them to the HST date gives a median HST-minus-Gaia displacement of −0.021″ east and +0.091″ north. The largest residual after applying this translation is 0.030″ across the four stars. This measures the Hubble-to-Gaia tie. The absolute ZTF-to-Hubble alignment retains possible systematic error. In the figure, the orange circle shows the 0.211″ 90th percentile of the nightly ZTF radial scatter. It describes the observations and does not constitute a formal confidence region. [Tie stars, individual residuals and measured offsets](results.json) · [Gaia cone](sources/gaia.json)

**Nearby catalog entries.** We accepted Hubble Advanced Products sources with positive flux measurements and catalog flags 0 or 1, which allow stellar and extended profiles. The nearest accepted F555W source, #822, lies 0.394″ from the registered transient position and has 0.15″ aperture AB magnitude 26.53 ± 0.20. The nearest in F814W, #876, lies 0.731″ away at 25.82 ± 0.13. Both carry the catalog’s extended-source flag. There are additional low-significance entries closer to the position. The small cyan circles identify accepted neighboring entries in each filter. [F555W source-catalog receipt](sources/f555w-point.ecsv.json) · [F814W source-catalog receipt](sources/f814w-point.ecsv.json)

| HST filter | Net count rate at the site (e⁻/s) | Comparison apertures as bright | Flux above local median | Empirical 5σ depth (aperture AB) |
| --- | ---: | ---: | ---: | ---: |
| F555W | 0.337 | 2 / 134 | 2.58 robust σ | 26.23 |
| F814W | 0.225 | 2 / 137 | 2.58 robust σ | 26.09 |

We used the local median sky in a 0.25–0.45″ annulus and the native calibration in each drizzled science image. Controls followed the same aperture procedure within 4.5″ on the arm, with broadly similar background and centers at least 0.35″ from accepted HAP sources. Among the 126 positions usable in both bands, none reached the target’s aperture flux in both. The bands share the same astrophysical background. These exploratory ranks and the local median absolute deviation describe the observed scene. The quoted depths apply to flux inside the aperture, with no point-source aperture correction, artificial-star completeness test or extinction correction. An individual progenitor identification would need accurate transient-to-pre-event registration plus evidence distinguishing the source from its neighbors. [Photometry and controls](hubble.py) · [Numerical results](results.json)

The extraction uses byte ranges from STScI’s public Hubble Advanced Products. It verifies the SCI extension, input row range, cutout hash and world coordinate system. North is up and east is left. The two primary filter mosaics, Gaia astrometry and full HAP source catalogs retain separate receipts. [F555W parent and cutout](sources/f555w.json) · [F814W parent and cutout](sources/f814w.json) · [STScI Hubble archive](https://archive.stsci.edu/missions-and-data/hst)

From the repository root, replay the local archive with NumPy and Astropy, and SciPy for the imported field utilities. Matplotlib is required only to rebuild the figure with `--plot`.

```sh
python3 experiments/compass-2025/batch-5/images/host/hst/hubble.py
```
