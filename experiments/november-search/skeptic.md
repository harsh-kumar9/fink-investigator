# Independent skeptic's first pass

These assessments use only the two entries in the [outage-audit packet](results/outage_followup_packet.json), before downstream catalog or image evidence.

## ZTF18abqclig

Inspect the source. Positive residuals on 63 of 315 quality-selected observing dates and the negative red-band anchor support a recurring change relative to the template. A variable superposed on a bright reference star is plausible. A persistent mismatch between the science image and template in a crowded field could recur on hundreds of dates as seeing changes. The fast-motion screen covers only stars passing its speed, significance and geometric thresholds. A 5″ empty VSX cone and a successful 3″ empty SIMBAD cone permit catalog aliases or a blended neighbor. The broker's failed lookup contributes no classification. Evidence IDs are `alert`, `history`, `motion`, `catalog` and `vsx-cone`.

Compare the centers and shapes of positive and negative residuals against Gaia positions on several dates, and measure the original-image brightness separately on the two groups of dates. Coincident pointlike changes accompanied by the expected direct-light variation would disfavor a subtraction artifact. The packet's `image` link supplies the anchor for that check. Keep this case open for archival follow-up.

## ZTF20aaothpc

Ten quality-selected positive dates make this a candidate brightening relative to the template. The packet supplies no dates for the first nine, so either a single episode or recurrent activity fits the counts. A previously reported transient with an absent broker name is a credible alternative, as is a recurrent source detected only above its template brightness. The narrow catalog searches do little to distinguish them. Evidence IDs are `alert`, `history`, `catalog`, `vsx-cone` and `motion`.

Retrieve the observation dates, full light curve and limiting magnitudes between detections. Compare the signed residual centroids with any source or extended host visible in the science and reference images, then query TNS by position and observing interval. A single rise and decay centered on a host would favor a transient. Well-separated episodes at a fixed stellar position would favor recurrent activity. Keep the case open for these checks. The only image evidence in this packet is a link (`image`).
