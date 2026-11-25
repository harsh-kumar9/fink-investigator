# Reference-image test for ZTF18abgobyb

The calibrated foreground screen flagged Gaia DR3 4293315765165489536 (VB 10). The anomaly packet has candid 2841171223915010017 and epoch JD 2460595.6712269. [SIMBAD identity query](https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+TOP+10+b.main_id%2Cb.otype%2Cb.pmra%2Cb.pmdec+FROM+basic+AS+b+JOIN+ident+AS+i+ON+i.oidref%3Db.oid+WHERE+i.id%3D%27Gaia+DR3+4293315765165489536%27) · [Fink object](https://ztf.fink-portal.org/ZTF18abgobyb)

Projecting the alert position onto the straight Gaia track gives 2018.363, about 2018-05-13. This is a geometric estimate of when the star crossed the residual position. Its perpendicular offset is 0.131″.

The alert's 18-frame ZTF reference spans 2018-03-21 through 2018-06-24 (JD 2458199.004259–2458293.814687, field 488, readout channel 39, filter 2). The inferred crossing falls **inside the actual coadd interval**. [Fink reference metadata](https://api.ztf.fink-portal.org/api/v1/objects?objectId=ZTF18abgobyb&columns=i%3Acandid%2Ci%3Ajd%2Ci%3Afield%2Ci%3Arcid%2Ci%3Afid%2Ci%3Ajdstartref%2Ci%3Ajdendref%2Ci%3Anframesref) · [ZTF field definitions](https://zwickytransientfacility.github.io/ztf-avro-alert/schema.html)

70 of 70 quality-selected detections through this alert are negative residuals. Together the sign, positional offset and coadd dates support subtraction at the star's former position. The effective epoch depends on the weights of individual reference exposures. A superposed variable remains possible. Multi-epoch image registration and the reference input exposures can discriminate between these explanations.

The estimate uses Gaia's catalog `pmra` (already multiplied by cos δ), `pmdec`, and reference epoch. It omits parallax and perspective acceleration. This targeted follow-up used an object outside the seven calibration controls. It does not independently measure the false-positive rate of the screen.
