# Pre-registered analysis plan

Committed before any RIR was generated or any HFD profile computed. The room
list, ground truth digitization, and pipeline code are being prepared in
parallel, but no detection results or regressions exist at the time this file
is written. This document fixes the analysis so that nothing downstream can be
quietly tuned to the answer.

## Primary test (confirmatory)

Reproduce Criterion III of Conaty 2011 exactly as defined:

1. For each of the eight included Lindau rooms (rooms 1 to 8, room 9 excluded
   as in the thesis), generate an image-source RIR with pyroomacoustics at
   fs = 44100 Hz, geometry and absorption from the reconstruction recorded in
   `inputs/rooms_reconstruction.md`, source (L/2, 2.0, 1.9), mic (L/2, 1.0, 1.9),
   RIR length = round(RT * fs) samples (assumption documented in
   `inputs/thesis_params.md`).
2. Compute the moving-window Higuchi FD profile, window 50 (51-sample segments,
   MATLAB convention), kmax 5, smoothed with a centered 199-point moving
   average (MATLAB smooth(l, 200) convention).
3. Tail statistics over the last 10 percent of the profile (MATLAB convention,
   round(len/10) + 1 samples): mean m and standard deviation s (ddof = 1).
4. Detection: first 1-based sample index >= 1000 where the smoothed profile
   >= m - 2s (Criterion III, k = 2).
5. Ordinary least squares regression of tmp50 percent (ground truth, in samples
   at 44100/s, from `inputs/ground_truth.json`, reconciled values) on the
   detected sample number, across the eight rooms. Report R squared, slope,
   intercept, and compare to the 2011 values (R squared 93.49 percent, slope
   0.3197, intercept 325).

The primary outcome is this single R squared. It is declared before results
are seen. No other criterion, window, kmax, tail fraction, or measure may be
substituted for it afterwards.

## Reported alongside the primary (also confirmatory, no selection)

- Full threshold sweep: k from 0 to 4 in steps of 0.25 (17 values), same
  pipeline, R squared for each. The whole curve is reported so that the 2011
  choice of the best of four criteria is visible in context.
- Criteria I, II, IV (k = 0, 1, 3) called out from the sweep for direct
  comparison with the 2011 table.

## Honest small-n statistics

n = 8 rooms (the thesis's own count). Declared in advance:

- Leave-one-out cross-validation of the Criterion III regression: refit on 7
  rooms, predict the held-out room, report predicted vs actual, LOOCV R
  squared (1 - PRESS/SStot) and LOOCV RMSE in ms.
- Bootstrap 95 percent confidence intervals (percentile, 10000 resamples of
  rooms with replacement, seed 20260722, resamples with fewer than 3 unique
  rooms discarded) on R squared and slope of the primary regression.
- Wide intervals are expected with n = 8 and will be stated as such.

## Sensitivity analysis (confirmatory of fragility, not of the effect)

Recompute the primary R squared (Criterion III only) over a small grid around
the arbitrary constants, one factor at a time from the primary configuration:

- tail fraction: 0.05, 0.10, 0.20
- window: 25, 50, 100
- kmax: 3, 5, 8
- RIR length: 0.75, 1.00, 1.25 times round(RT * fs) (the length is our
  assumption, so its influence must be shown)
- detection start index: 500, 1000, 2000

Report the table. Stability is reassurance; large swings are themselves the
finding.

## Secondary, exploratory only

Clearly labeled exploratory, reported separately, never merged into the
primary claim: (a) sample entropy profile with the same windowing, (b) Katz
fractal dimension profile. Same detection rule, k = 2 only. No sweep over
their parameters, no selection of a better-looking variant.

## Ground truth

tmp50 percent per room is digitized from thesis figures (the thesis never
tabulates it) and reconciled against the four reported 2011 regressions. The
digitized-vs-reconciled spread is carried as ground truth uncertainty and the
primary regression is additionally run with the raw digitized values as a
robustness check. If the two disagree materially, both results are reported.

## Verdict rules

- "Reproduces" requires the primary R squared within roughly 10 points of
  93.49 percent with a same-sign slope of similar magnitude.
- "Reproduces but fragile" if the point estimate is high but LOOCV collapses,
  bootstrap intervals span near-zero, or the sensitivity table swings by tens
  of points under small constant changes.
- "Does not reproduce" otherwise. All three verdicts are acceptable outcomes.
