# Thesis parameters: Conaty 2011, "The Fractal Dimension of Room Impulse Responses as a Predictor for the Perceptual Mixing Time"

Extracted from `inputs/FINAL_thom_conaty.pdf` (101 pages) and its text dump
`inputs/thesis_extracted_text.txt`. Page references below give the printed
page number as it appears on the page, with the PDF text-dump marker in
parentheses (marker `PAGE N` corresponds to printed page N minus 11, e.g.
`PAGE 55` is printed p.44). Every number carries a provenance tag:

- `exact`: transcribed from a thesis table or the thesis text.
- `digitized`: measured programmatically from a figure by
  `inputs/digitize_figures.py`; audit overlays in `inputs/digitization_debug/`.
- `inferred`: not stated in the thesis, derived here and flagged as such.
- `missing`: could not be recovered.

## 1. Room table (Lindau et al. 2010 rooms)

Volumes, average absorption and RT are exact from Table 4.3.1.1, printed p.31
(PAGE 42); verified against the native-resolution embedded table image
(`digitization_debug/table4311_xref254.png`). Room names are exact from the
panel titles of the floor plan figure, Figure 4.3.1.1, printed p.31
(`digitization_debug/floorplans_fig4311_xref256.png`). Footprint length and
width are digitized from those floor plans (about plus minus 0.5 m); heights
are inferred as V/(L*W) since the thesis never tabulates room dimensions.

| Room | Name | V (m3) | alpha_ave | RT (s) | L x W (m, digitized) | H = V/LW (m, inferred) | In 2011 analysis |
|---|---|---|---|---|---|---|---|
| 1 | E. Studio | 216 | 0.36 | 0.39 | 9.6 x 7.6 | 2.96 | yes |
| 2 | EN-111 | 224 | 0.26 | 0.62 | 7.1 x 7.2 | 4.38 | yes |
| 3 | EN-190 | 182 | 0.17 | 0.79 | 8.6 x 7.2 | 2.94 | yes |
| 4 | H-0104 | 3300 | 0.28 | 1.15 | 26.0 x 20.7 | 6.13 | yes |
| 5 | HE-101 | 5179 | 0.23 | 1.67 | 27.2 x 25.9 | 7.35 | yes |
| 6 | Teldex Studio | 3647 | 0.20 | 1.83 | 28.3 x 16.1 | 8.00 | yes |
| 7 | UoA concert hall | 8298 | 0.33 | 1.52 | 38.2 x 17.3 | 12.56 | yes |
| 8 | TUB Audimax | 8500 | 0.23 | 2.08 | 30.1 x 29.3 | 9.64 | yes |
| 9 | JC church | 7417 | 0.23 | 2.36 | 25.8 x 15.0 (bbox) | n/a | NO (excluded) |

An independent digitization of the same floor plans exists in
`inputs/digitized_rooms.json` (script `inputs/digitize_rooms.py`, from the
geometry work stream), calibrated on the dotted gridlines rather than the
axis labels. The two agree within 0.5 m on most rooms; the widths of rooms
7, 8 and 9 differ by 0.7 to 1.4 m. Treat the footprints as good to about
plus minus 1 m; for ISM geometry prefer the gridline calibrated values.

Room 9 was excluded because it "deviated from the typical shoebox shape ...
beyond the modelling capability of the ISM" (printed p.38, PAGE 49). Its floor
plan is non rectangular (cross shaped); the L x W above is only its bounding
box. Table 4.3.1.1 also gives column averages: small V 207 m3, medium V
4042 m3, large V 8072 m3; row averages alpha (RT): 0.32 (1 s), 0.24 (1.45 s),
0.2 (1.66 s).

Lindau measurement context (printed p.32, PAGE 43, exact): BRIRs measured
with the FABIAN head and torso simulator, 3 way dodecahedron source on a
stage at one end of each room, horizontal head orientations within plus and
minus 80 degrees in 1 degree steps; 24 subjects, 20 trials per room, only 10
subjects retained for analysis (printed p.32, PAGE 43).

## 2. Perceptual ground truth: per room tmp50%

The thesis NEVER tabulates the tmp50% values it regressed against. They were
recovered here by programmatic digitization of two independent figures.

### 2.1 Primary: Figure 6.2 scatter panels (printed p.46, PAGE 57)

Four panels (Criteria I to IV) plot the same eight tmp50% values (y, in
samples at 44100 per s) against the known detected sample numbers (x, Table
6.1b). Markers were detected by exact colour (the marker glyphs are flat
pure blue or purple, distinct from text antialiasing), axes calibrated from
tick marks (calibration residual under 0.25 px, about 1 sample), and markers
assigned to rooms by the known x values (Hungarian assignment, near
coincident x pairs disambiguated across panels by y consensus). The four
panels agree to within 0.5 samples per room, so the per panel spread is not
a useful uncertainty measure; realistic uncertainty is about plus minus 5
samples (0.11 ms) from marker centroid and calibration bias.

Reconciled values: penalized least squares over the 8 unknowns against the
12 reported regression statistics (slope, intercept, R2 for four criteria,
figure annotation versions), after dropping the two statistics that are
unreproducible from any consistent data set (see Section 6). Reconciliation
moved every value by at most 2.3 samples, confirming the digitization.

| Room | tmp50% digitized (samples) | (ms) | reconciled (samples) | (ms) | dig minus rec (samples) | round ms hypothesis (inferred) |
|---|---|---|---|---|---|---|
| 1 | 1656.0 | 37.55 | 1654.8 | 37.52 | 1.2 | 37.5 |
| 2 | 1325.0 | 30.04 | 1324.0 | 30.02 | 1.0 | 30.0 |
| 3 | 993.8 | 22.54 | 992.3 | 22.50 | 1.6 | 22.5 |
| 4 | 2643.7 | 59.95 | 2645.7 | 59.99 | -2.0 | 60.0 |
| 5 | 2328.7 | 52.80 | 2326.4 | 52.75 | 2.3 | 52.8 |
| 6 | 3086.9 | 70.00 | 3087.4 | 70.01 | -0.5 | 70.0 |
| 7 | 4309.6 | 97.72 | 4310.8 | 97.75 | -1.2 | 97.8 |
| 8 | 2868.0 | 65.03 | 2867.3 | 65.02 | 0.7 | 65.0 |
| 9 | missing | missing | missing | missing | n/a | n/a |

Round ms hypothesis (inferred, NOT thesis text): every digitized marker sits
within 3.5 samples (under 1 px) of a round ms value times 44.1, namely
37.5, 30.0, 22.5, 60.0, 52.8, 70.0, 97.8, 65.0 ms, and these values
reproduce the reported regressions as well as the raw digitization does.
The endpoints 37.5 and 97.8 ms are exactly the figures quoted in the thesis
abstract (printed p.3, PAGE 14). This is very likely the exact data vector
the thesis used (values read from Lindau, converted to samples at 44.1 per
ms).

### 2.2 Secondary: Figure 4.3.2 error bar plot (printed p.33, PAGE 44)

Independent digitization of the Lindau figure reproduced in the thesis: mean
tmp50% (small widened dot on each bar; in all 9 cases the dot is present) and
95 percent CI (serif caps). Uncertainty about plus minus 1 ms (2 px). This
also provides Room 9, which is absent from Figure 6.2.

| Room | mean (ms) | 95 percent CI (ms) | CI midpoint (ms) |
|---|---|---|---|
| 1 | 36.3 | [26.5, 46.4] | 36.4 |
| 2 | 29.4 | [22.7, 36.6] | 29.6 |
| 3 | 24.7 | [15.5, 33.4] | 24.5 |
| 4 | 54.5 | [42.8, 66.8] | 54.8 |
| 5 | 46.4 | [36.8, 56.3] | 46.5 |
| 6 | 63.9 | [48.6, 79.1] | 63.9 |
| 7 | 91.2 | [57.8, 124.1] | 90.9 |
| 8 | 63.4 | [40.4, 86.9] | 63.6 |
| 9 | 73.7 | [51.6, 95.6] | 73.6 |

IMPORTANT: the Figure 6.2 y values (the ones the 2011 regressions used) do
NOT equal an accurate reading of Figure 4.3.2. Rooms 4 to 8 are 4 to 7 ms
higher in Figure 6.2 than the means digitized from Figure 4.3.2 (e.g. Room 7:
97.7 vs 91.2 ms; Room 5: 52.8 vs 46.4 ms). Rooms 1 to 3 agree within about
1 to 2 ms. Whatever source the thesis used for its y vector, it was not a
precise read of this error bar figure. For reproducing the 2011 regressions,
use the Section 2.1 values.

### 2.3 Sanity anchors

- Figure 5.3 (printed p.42, PAGE 53) draws the Room 1 tmp50% as a vertical
  line; digitized at 1647.2 samples = 37.35 ms. Matches Room 1 above
  (1656 samples, 37.55 ms) within about 9 samples (2 px of that figure).
- Abstract (printed p.3, PAGE 14) quotes tmp50% "ranging from 37.5ms to
  97.8ms". The measured maximum (Room 7, 97.7 ms) matches, but the measured
  MINIMUM in Figure 6.2 is Room 3 at about 22.5 ms (about 994 samples), not
  37.5 ms. 37.5 ms is Room 1. The abstract's "range" statement does not
  describe the actual min across rooms in the thesis's own data; do not use
  it as a constraint.

## 3. Method constants (2011 pipeline)

All exact, from thesis text and the MATLAB code in Appendix 9.1 (printed
pp.53 to 56, PAGES 64 to 67) unless noted.

- Sample rate: fs = 44100 Hz (printed p.39, PAGE 50).
- FD profile: HMW.m moving window, window parameter y = 50 samples
  (about 1 ms; printed p.39). NOTE: the MATLAB segment `k(i:(i+y))` is 51
  samples inclusive. The window advances 1 sample per step; the RIR is zero
  padded with y samples at the end; loop breaks when the window would pass
  the end.
- Higuchi: hfd.m (Salai Selvam implementation), called by HMW.m with NO kmax
  argument, so the default kmax = 5 applies. FD is the slope of a first
  order polyfit of ln(L(k)) against ln(1/k), k = 1..5.
- Smoothing: `sl = smooth(l, 200)` before analysis (printed p.39 says
  "moving window of 200 samples"). MATLAB smooth is a centred moving
  average and reduces even spans to the next lower odd value, so the
  effective span is 199.
- Tail statistics: last 10 percent of the profile,
  `x((length(x) - round(length(x)/10)) : length(x))`, i.e. the last 10
  percent plus one sample; mean and standard deviation of that slice
  (FDthres.m).
- Criteria thresholds (printed p.41, PAGE 52; FDthres.m):
  I = tail mean; II = mean minus 1 SD; III = mean minus 2 SD;
  IV = mean minus 3 SD.
- Detection: first sample i, scanning i = 1000 upward in steps of 1, where
  profile(i) >= threshold ("Start at sample 1000 to avoid FD error at start
  of profile"). Detection is the FIRST crossing at or after sample 1000.
- ISM: MATLAB re-code of the Fortran in Allen and Berkley 1979 (sroom.m,
  Appendix 9.1.4). Inputs: room dimensions, source and receiver Cartesian
  positions, per wall reflection coefficients beta, fs, response length.
- Source and receiver (thesis text, printed p.39, PAGE 50): receiver at
  x = L/2, y = 1 m, z = 1.9 m; source at x = L/2, y = 2 m, z = 1.9 m
  (x is the length axis). See discrepancy note in Section 6.
- Reflection coefficients: alpha_ave "converted to beta and inputted as
  each wall's coefficient" (printed p.39). The conversion formula is NOT
  stated. beta = sqrt(1 - alpha) is the standard Allen and Berkley
  convention, but this is an assumption to be tested, not thesis text. The
  GUI example (Figure 5.1a) shows all betas = 0.9 for a 10x10x10 m example
  room with fs 44100 and response length 0.256 s.
- RIR length: UNCERTAIN, see Section 6.

## 4. 2011 results to reproduce

### 4.1 Table 6.1a: tail statistics of the FD profiles (printed p.44, PAGE 55; exact)

| Room | tail mean (FD) | tail SD (FD) | Thr I | Thr II | Thr III | Thr IV |
|---|---|---|---|---|---|---|
| 1 | 2.0364 | 0.0252 | 2.0364 | 2.0112 | 1.9860 | 1.9608 |
| 2 | 2.0450 | 0.0265 | 2.0450 | 2.0185 | 1.9921 | 1.9656 |
| 3 | 2.0250 | 0.0321 | 2.0250 | 1.9927 | 1.9605 | 1.9284 |
| 4 | 1.9833 | 0.0394 | 1.9833 | 1.9439 | 1.9046 | 1.8652 |
| 5 | 1.9782 | 0.0404 | 1.9782 | 1.9377 | 1.8973 | 1.8569 |
| 6 | 1.9866 | 0.0432 | 1.9866 | 1.9435 | 1.9003 | 1.8571 |
| 7 | 1.9682 | 0.0394 | 1.9682 | 1.9288 | 1.8894 | 1.8500 |
| 8 | 1.9644 | 0.0506 | 1.9644 | 1.9138 | 1.8631 | 1.8125 |

(Thesis table prints Thr III for Room 1 as "1.986" and Thr I equals the
mean by definition.)

### 4.2 Table 6.1b: detected sample numbers (printed p.44, PAGE 55; exact)

Rooms 1 to 8 in order; samples at 44100 per s (ms in the thesis table are
sample/44.1):

- Criterion I: 6277, 4573, 6056, 13724, 7996, 12366, 14850, 16194
- Criterion II: 5772, 3197, 4335, 7027, 7971, 9486, 12308, 12571
- Criterion III: 3136, 3170, 2494, 7013, 7906, 8060, 11708, 8465
- Criterion IV: 2053, 3049, 1567, 7004, 7661, 2052, 2192, 5209

(The running text on printed p.43, PAGE 54, says "4571 (103.70ms)" for Room
2 Criterion I; 4573/44.1 = 103.70 ms, so the table value 4573 is the
consistent one and 4571 is a typo.)

### 4.3 Regressions, tmp50% (samples) against detected sample number

Table 6.2 (printed p.46, PAGE 57; exact transcription):

| Criterion | equation | R2 |
|---|---|---|
| I | tmp50% = 0.2090 * I + 331 | 72.29 percent |
| II | tmp50% = 0.2760 * II + 293 | 78.82 percent |
| III | tmp50% = 0.3197 * III + 325 | 93.49 percent |
| IV | tmp50% = 0.0517 * IV + 2202 | 00.14 percent |

Figure 6.2 panel annotations (exact transcription; more digits, produced by
the MATLAB curve fitting toolbox):

| Criterion | equation | R2 |
|---|---|---|
| I | y = 0.2019x + 330.9 | 72.29 percent |
| II | y = 0.276x + 293.2 | 78.82 percent |
| III | y = 0.3197x + 324.9 | 93.49 percent |
| IV | y = 0.05172x + 2202 | 00.14 percent |

Our regression check, fitting the DIGITIZED y against the known x:

| Criterion | fitted equation | fitted R2 | verdict |
|---|---|---|---|
| I | y = 0.2017x + 333.5 | 72.26 percent | matches figure (0.2019), so Table 6.2's 0.2090 is a digit transposition typo |
| II | y = 0.2758x + 240.9 | 78.87 percent | slope and R2 match; reported intercept 293.2 is NOT reproducible (see Section 6) |
| III | y = 0.3195x + 326.8 | 93.52 percent | matches |
| IV | y = 0.0517x + 2202.4 | 1.37 percent | slope and intercept match; reported R2 00.14 percent is NOT reproducible (actual about 1.37 percent) |

### 4.4 Lindau comparison context (exact)

- Best Lindau model predictor: tmp50% = 20.08 * V/S + 12 (ms), R2 = 81.5
  percent (mean free path; printed p.33, PAGE 44).
- Other Lindau R2 for tmp50%: reflection density sqrt(V) model 78.6 percent,
  V 77.4 percent, S 73.3 percent, RT 53.4 percent, alpha_ave 0.8 percent
  (printed p.34, PAGE 45; S = 73.3 percent from Table 6.3, printed p.47).
- Empirical predictors (Lindau study): Abel and Huang I 74.7 percent,
  Hidaka 1 kHz 57.3 percent, Hidaka 500 Hz 56.5 percent, Abel and Huang II
  50.7 percent, Stewart and Sandler 37.6 percent (printed p.34).
- Thesis ranking table 6.3 (printed p.47): Conaty III 93.49, Conaty II
  78.82, Abel and Huang I 74.7, Conaty I 72.29, Hidaka 1k 57.3, Hidaka 500
  56.5, Abel and Huang II 50.7, Stewart and Sandler 37.6 (percent).

## 5. Digitization deliverables

- `inputs/digitize_figures.py`: deterministic script (identical output
  across runs) that extracts the embedded figure images at native
  resolution from the PDF, digitizes Figures 6.2, 4.3.2, 4.3.1.1 and 5.3,
  runs the regression checks and reconciliation, and writes
  `inputs/ground_truth.json`.
- `inputs/digitization_debug/`: native resolution figure extractions,
  per panel overlay images with every detected marker circled and labelled
  by room (`overlay_fig62_c1..4.png`), the error bar overlay with detected
  caps and means (`overlay_fig432.png`), and the raw numeric results
  (`digitization_results.json`).

## 6. Known discrepancies, uncertainties and gaps

1. Criterion I slope: Table 6.2 says 0.2090, Figure 6.2 annotation says
   0.2019. The digitized fit gives 0.2017, so the FIGURE is correct and the
   table transposed digits.
2. Criterion II intercept: both table (293) and figure (293.2) report an
   intercept that is inconsistent with the reported slope 0.276, the
   reported R2 78.82 percent, and the marker positions in the figure
   itself. The self consistent intercept is about 239 to 241. All three
   other panels imply the same mean y (about 2401 samples); Criterion II's
   reported intercept implies 2455. Treated as a thesis typo; the
   reconciliation dropped this constraint.
3. Criterion IV R2: reported as 00.14 percent, but the marker data with the
   reported slope and intercept give about 1.37 percent. Unreproducible;
   dropped in reconciliation. (Possibly a misreading of 0.0137.)
4. Room 2 Criterion I detected sample: text says 4571, table says 4573; the
   quoted ms value supports 4573.
5. Source and receiver positions: thesis text (printed p.39) places the
   receiver at y = 1 and source at y = 2; the GUI screenshot (Figure 5.1a,
   printed p.40) shows source (5, 1, 1.9) and receiver (5, 2, 1.9), i.e.
   swapped. Because the image source RIR is reciprocal under exchange of
   source and receiver, the swap does not change the RIR; either reading
   gives the same result.
6. RIR length per room: UNCERTAIN and never stated. Evidence: GUI default
   response length 0.256 s (also in the GUI code, sroomgui.m); Figure 5.1b
   plots the Room 1 RIR to 0.2 s; Figure 5.3's Room 1 FD profile extends to
   about 17200 samples, about 0.39 s = RT * fs for Room 1. The tail
   statistics depend on profile length, so a reproduction should treat the
   length as a free parameter (RT * fs per room is the best supported
   guess for the analysis stage).
7. alpha to beta conversion formula not stated (Section 3); beta =
   sqrt(1 - alpha_ave), identical on all six walls, is the standard
   assumption to test.
8. Room dimensions are not tabulated; footprints digitized from small floor
   plans (plus minus about 0.5 m) and heights inferred from V. Errors here
   propagate into any ISM reproduction.
9. Figure 6.2 y values disagree with Figure 4.3.2 means by up to about 7 ms
   for rooms 4 to 8 (Section 2.2). The thesis's regression target is the
   Figure 6.2 data.
10. The abstract's tmp50% range "37.5ms to 97.8ms" matches Rooms 1 and 7,
    not the actual minimum (Room 3, about 22.5 ms).
11. Room 9 (JC church): no tmp50% appears in Figure 6.2 and no criteria
    were computed for it (excluded, non shoebox). Its only recovered ground
    truth is the Figure 4.3.2 digitization (mean 73.7 ms, CI [51.6, 95.6]).
12. tmp95% values (Lindau) are shown only in cascade plots (Figures 4.3.2.1
    and 4.3.2.2, printed pp.35 to 36) at small scale; not digitized here:
    missing.
13. The GUI code default fs in sroomgui.m is 8000; the screenshot and the
    thesis text confirm 44100 was used for the actual RIRs.
