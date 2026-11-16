# Pre-registered follow-up: window-scaled HFD mixing-time detection on measured RIRs

Written and committed 2026-07-26, BEFORE: (a) the author's recovery of any
2011 archive material, (b) any analysis of any measured RIR, and (c) the
ISM window-sweep groundwork below (which runs only after this file is
committed). Prior knowledge declared in full: everything already in this
repository, including the single sensitivity point that window = 100
raised the band-limited track's Criterion III R squared from 47.2 to 74.7
percent, which is the observation that motivates this study. Nothing else
about window scaling has been computed at the time of writing.

## Motivation (one paragraph)

The 2011 method's feature, moving-window Higuchi FD with a 51-sample
window, measures spike sparsity at a scale shorter than the impulse
spread of any band-limited system, so it saturates instantly on measured
or band-limited simulated RIRs (REPORT.md, docs/ISM_COMPARISON.md). The
follow-up hypothesis, due to the author, is that the failure is one of
SCALE, not of concept: HFD windows longer than the system's impulse
spread should again see echo-density texture rather than kernel texture.

## Hypotheses

- H1 (feasibility gate): on measured RIRs, the smoothed HFD profile
  computed with window length W well above the impulse spread exhibits a
  measurable rise: (tail mean minus profile minimum) > 3 times the tail
  SD, in at least three quarters of tested rooms at some W in the grid.
  If H1 fails, the method is dead on measured data and H2 is not tested.
- H2 (prediction): the floor-free detection defined below, at the
  pre-fixed window W* and k = 2, predicts perceptual mixing time tmp50
  across rooms with LOOCV R squared > 0.5 and Spearman rho > 0.7.
  H2 requires perceptual ground truth (the Lindau BRIR set with its
  tmp50 values, or a new listening test); it cannot be tested on
  unlabeled corpora.

## Fixed analysis pipeline (no free constants at test time)

1. Profile: Higuchi FD, kmax = 5, moving window of W samples (segment
   W + 1), same code as src/mixtime/profile.py.
2. Smoothing: centered moving average with span = 4 * W (MATLAB
   smooth semantics, nan tolerant), scaling the original 200-for-50
   ratio.
3. Detection (floor-free): locate the profile minimum over finite
   samples AT OR AFTER THE ONSET, where the onset is defined (Amendment
   1) as the first sample whose absolute RIR amplitude reaches 5 percent
   of the RIR's peak absolute amplitude; threshold = tail mean minus 2
   times tail SD over the last 10 percent; detected sample = first
   crossing at or after the argmin. k = 2 is retained from the 2011
   Criterion III, deliberately unchanged.
4. Window grid, all reported, none selected post hoc: W in {50, 75,
   100, 150, 200, 300, 500, 750, 1000, 1500, 2000}.
5. W* (the single confirmatory window for H2) is fixed by PROCEDURE,
   not by hand, from the ISM groundwork that runs after this commit:
   W* = the smallest W whose Criterion III R squared lies within 10
   points of the maximum over the grid on BOTH ISM tracks
   simultaneously. Once the groundwork script first runs, W* is frozen
   and recorded in its output; any later change would be a logged
   amendment.
6. Primary metric for H2: LOOCV R squared (1 minus PRESS over SStot).
   Secondary: Spearman rho, OLS R squared, RMSE in ms. The full W grid
   is reported alongside as context, exploratory.

## Data, in order of preference

1. The Lindau et al. measured BRIR set with its per-room tmp50 values
   (full H2 test). Left ear, frontal head orientation, energy onset
   aligned; details to be amended (logged) once the data's actual
   format is known.
2. Failing that, any measured RIR corpus WITH perceptual mixing-time
   labels (none known at time of writing).
3. Failing both, open measured corpora without labels (openAIR, Aachen
   AIR, MIT IR Survey): H1 only, explicitly no correlational claims.

## Honesty rules

- n will be small in every scenario; bootstrap CIs reported, wide
  intervals stated as such.
- Every window in the grid is reported in every output. No best-of
  quoting: the confirmatory number is W* alone, chosen by the frozen
  procedure above.
- If H1 passes and H2 fails, that is the result. If archive material
  recovered later suggests different constants, testing them is a new,
  amended analysis, not a revision of this one.
- ISM-based results (including the groundwork) are hypothesis
  generation only and can never confirm H2.

## Amendment 1 (2026-07-27, logged before any H2 data exists)

Triggered by external code review noting two gaps; recorded here per the
honesty rules. No measured labeled data has been seen.

1. Onset definition. The original text said the floor-free minimum is
   taken "after the direct sound" without operationalizing it, and the
   first groundwork run took the argmin over the whole profile (some
   detections landed at sample 1). Fixed definition: the onset is the
   first sample whose absolute RIR amplitude reaches 5 percent of the
   RIR's peak absolute amplitude, and the argmin is restricted to
   profile samples at or after that index. The groundwork sweep is
   re-run under this rule and W* re-frozen from the re-run; both the
   original (unguarded) and amended outputs remain in git history.
2. Metric mismatch, clarified not changed: W* is deliberately frozen on
   OLS R squared, as a plateau-location point estimate over the grid,
   while H2's confirmatory metric remains LOOCV R squared. Freezing the
   window on one metric and confirming on another avoids selecting the
   window that happens to maximize the confirmatory statistic, which
   would reintroduce the in-sample selection this project exists to
   avoid.

## Groundwork permitted before measured data (runs after this commit)

analysis/exploratory/window_sweep.py: the W grid above on the eight
EXISTING ISM RIRs (both tracks, both the original fixed-start rule and
the floor-free rule), reporting R squared, Spearman, LOOCV and
floor-pinned counts per cell. Purpose: freeze W* by the stated
procedure and characterize the plateau shape. Labeled
hypothesis-generating in all outputs.
