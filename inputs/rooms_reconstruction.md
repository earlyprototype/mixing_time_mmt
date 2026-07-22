# Room geometry reconstruction for the nine Lindau et al. 2010 rooms

Step 3 setup for the 2011 thesis reproduction. Tabulated data come from
thesis Table 4.3.1.1 (page 42, transcribing Lindau et al. 2010). Only the
volume is tabulated, so footprints were digitized from the floor plans of
Figure 4.3.1.1 on the same page, and height was derived from the volume.

Code: geometry in `src/mixtime/rooms.py`, digitizer in
`inputs/digitize_rooms.py` (rerunnable, writes `inputs/digitized_rooms.json`),
RIR generator in `src/mixtime/ism.py`, RT estimator in
`src/mixtime/rt_check.py`, sanity script in `analysis/check_rooms.py`.

## Digitization method

1. Figure 4.3.1.1 is embedded in `inputs/FINAL_thom_conaty.pdf` (page 42,
   zero-based index 41) as a raster image, xref 256, native size 723 x 755 px.
   That native raster was extracted with pymupdf and used directly, so no
   resampling was involved.
2. The nine subplot axis frames were located as long solid dark lines
   (grayscale < 128, run length > 120 px).
3. Axis calibration per subplot: each subplot has dotted gridlines at data
   coordinates -20, 0 and +20 m on both axes (visible at grayscale
   threshold < 210). Their pixel positions give a per-axis linear scale.
   Measured scales: x (length) 3.95 px/m in all subplots; y (width) 3.475
   to 3.675 px/m depending on subplot row. The subplot frames extend beyond
   plus and minus 20 m, so the frames themselves were NOT used for
   calibration.
4. Room outlines were segmented as solid black strokes (grayscale < 128
   with contiguous run >= 8 px, which excludes the dotted gridlines). Edge
   positions were taken as stroke-centre centroids of the edge pixel
   clusters, then converted to metres.
5. L = horizontal extent, W = vertical extent, H = V / (L * W).

Estimated digitization uncertainty: about 1 px per edge, roughly 0.3 m in
L or W, which propagates to a few percent in H and S.

## Assumptions (all of them)

- A1. The dotted gridlines in Figure 4.3.1.1 mark -20, 0, +20 m exactly.
  Supported by the tick labels "-20 0 20" on both axes.
- A2. Each room is a shoebox of height H = V / (L * W), with the tabulated
  volume exact and the digitized footprint approximate. Sanity band
  2 m < H < 25 m; all nine derived heights fall inside it, so NO
  proportional adjustment of L and W was needed for any room (the
  adjustment rule was therefore never exercised).
- A3. Room 9 (JC church) is non-rectangular (notched outline). Its L and W
  are the bounding rectangle of the outline, so the footprint area is
  overestimated and the derived H = 20.12 m underestimates the true mean
  height. Room 9 is flagged included=False, as the thesis excluded it.
- A4. RIR length default is round(rt_s * 44100) samples per room. For
  Room 1 that is 17199 samples, matching the roughly 17200-sample profile
  length in thesis Figure 5.3. This is an assumption; the thesis states the
  ISM GUI took a length input but never tabulates the per-room value.
- A5. Source at (L/2, 2.0, 1.9) m and microphone at (L/2, 1.0, 1.9) m,
  from thesis Section 5.1 text. The thesis GUI screenshot (Figure 5.1a)
  swaps the source and receiver y values; for a pure ISM the path is
  reciprocal so the RIR is identical, and the text convention is used.
  Positions are clipped to keep a 0.5 m margin from every wall; no room in
  this set triggers the clip (smallest height 2.95 m leaves 1.05 m above
  the 1.9 m source and mic height).
- A6. Absorption: pyroomacoustics is given the energy absorption
  coefficient alpha_ave on all six surfaces. The thesis fed the pressure
  reflection coefficient beta = sqrt(1 - alpha_ave) to an Allen and Berkley
  code; energy reflection beta^2 = 1 - alpha_ave, so the two are equivalent
  in energy terms.
- A7. Pure image-source model: no air absorption, no ray tracing, no
  randomized image method, speed of sound 343 m/s (pyroomacoustics
  default). The computation is fully deterministic; there is nothing
  random to seed.
- A8. max_order per room: ceil(c * t_len / min_dim) + 3, capped at 100.
  The cap exists to keep runtime and memory tractable for large rooms; with
  the default lengths the per-room requirements are 42 to 95 (see table),
  so the cap never truncates a default request. The generator also returns
  the raw untrimmed RIR so decay coverage can be verified (for Room 1 the
  raw RIR covers 3.4 times the analysis window).

## Reconstructed geometry and Sabine sanity check

Sabine RT = 0.161 V / (S alpha) with S = 2(LW + LH + WH). With V and alpha
fixed, Sabine RT depends only on S, so ratio = Sabine/tabulated measures
how wrong the reconstructed proportions are (or how non-Sabine the real
room was). S_match is the surface area that WOULD make Sabine equal the
tabulated RT. Dimensions were NOT tuned to force a match.

| Room | Name | V m3 | alpha | RT s | L m | W m | H m | S m2 | Sabine RT s | ratio | S_match m2 | order | included |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | E. Studio | 216 | 0.36 | 0.39 | 9.48 | 7.62 | 2.99 | 246.7 | 0.392 | 1.00 | 247.7 | 48 | yes |
| 2 | EN-111 | 224 | 0.26 | 0.62 | 6.95 | 7.18 | 4.49 | 226.7 | 0.612 | 0.99 | 223.7 | 51 | yes |
| 3 | EN-190 | 182 | 0.17 | 0.79 | 8.60 | 7.18 | 2.95 | 216.5 | 0.796 | 1.01 | 218.2 | 95 | yes |
| 4 | H-0104 | 3300 | 0.28 | 1.15 | 25.69 | 20.71 | 6.20 | 1639.7 | 1.157 | 1.01 | 1650.0 | 67 | yes |
| 5 | HE-101 | 5179 | 0.23 | 1.67 | 26.83 | 25.86 | 7.46 | 2174.2 | 1.667 | 1.00 | 2170.8 | 80 | yes |
| 6 | Teldex | 3647 | 0.20 | 1.83 | 28.35 | 16.14 | 7.97 | 1624.3 | 1.807 | 0.99 | 1604.3 | 82 | yes |
| 7 | UoA concert hall | 8298 | 0.33 | 1.52 | 37.72 | 16.46 | 13.37 | 2690.0 | 1.505 | 0.99 | 2663.4 | 43 | yes |
| 8 | TUB Audimax | 8500 | 0.23 | 2.08 | 29.74 | 27.89 | 10.25 | 2840.1 | 2.095 | 1.01 | 2860.6 | 73 | yes |
| 9 | JC church | 7417 | 0.23 | 2.36 | 25.82 | 14.28 | 20.12 | 2350.7 | 2.209 | 0.94 | 2200.0 | 60 | NO |

"order" is the per-room default max_order at the default RIR length.

Reading of the ratios: rooms 1 to 8 agree with Sabine to within 1 percent.
This is almost certainly because Lindau et al. derived alpha_ave from the
measured RT, V and S via Sabine's formula, so the agreement validates the
digitized S (and hence the proportions) rather than the acoustics. Room 9's
0.94 reflects the bounding rectangle overestimating the footprint of the
notched plan.

## End-to-end check, Room 1 (from analysis/check_rooms.py)

- max_order 48 (uncapped requirement also 48), generation 0.1 s.
- RIR length 17199 samples (0.390 s), raw untrimmed length 58548 samples,
  so image sources cover 3.4x the analysis window. The pyroomacoustics
  global fractional-delay offset (40 samples) is removed before trimming.
- Raw Schroeder level at the window end: -48.1 dB. Not -60 dB because the
  window is one TABULATED RT long while the ISM RIR decays with its own RT
  (see next point); the decay inside the window is smooth and fully
  developed, with no truncation artifact.
- Schroeder RT60 (backward integration, T20 line fit from -5 to -25 dB,
  extrapolated): 0.454 s versus tabulated 0.39 s, ratio 1.164.

## Warnings and caveats

- The ISM RT differs from Sabine by construction. A shoebox ISM with
  uniform absorption is not a diffuse field; axial and tangential image
  families decay more slowly than the Sabine average, so the measured
  Schroeder RT (0.454 s for Room 1) exceeds the tabulated 0.39 s by about
  16 percent even though the Sabine RT of the geometry matches to 1
  percent. The same behaviour is expected of the thesis's own Allen and
  Berkley RIRs, so this is faithful to what is being reproduced, not a bug.
- alpha_ave = 0.36 (Room 1) and 0.33 (Room 7) are large enough that the
  Sabine formula itself is questionable (Eyring would predict noticeably
  shorter RTs). Since the tabulated alphas appear Sabine-derived, they are
  used as-is.
- Room 2's digitized footprint (6.95 x 7.18 m) gives H = 4.49 m, notably
  taller than rooms 1 and 3. Plausible for a seminar room but it is the
  least certain height in the small-room group.
- Room 9 (excluded): bounding-rectangle footprint, H underestimated, and
  a shoebox model of it would misrepresent the notched plan. Kept only for
  completeness with included=False.
- No positions were clipped in any of the nine rooms.
