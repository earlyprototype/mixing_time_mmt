# Corrections log: places where this project overrides the printed thesis

Authorized by the author (session instruction, 2026-07-25): corrections
that contradict the printed 2011 text may be actioned, provided every one
is logged. This file is that log. Rules: the original printed value is
always preserved here, the corrected value is the operative one used in
tables and analysis, and each entry states the evidence. New corrections
must be appended here in the same commit that applies them.

## C1. Table 6.2, Criterion I regression slope

- Printed: tmp50% = 0.2090 * I + 331
- Corrected (operative): tmp50% = 0.2019 * I + 330.9
- Evidence: the Figure 6.2 panel annotation reads 0.2019x + 330.9, and an
  independent least squares fit of the digitized markers against the known
  Table 6.1b x values gives 0.2017x + 333.5 with R squared 72.26 percent,
  matching the printed R squared (72.29). The 0.2090 in Table 6.2 is a
  digit transposition of the figure's 0.2019. Applied in
  inputs/thesis_params.md Section 5 and analysis/analyze.py (shown as the
  corrected comparator alongside the printed value).

## C2. Table 6.2 and Figure 6.2, Criterion II regression intercept

- Printed: 293 (table), 293.2 (figure annotation)
- Corrected (operative): about 241 (self-consistent value 240.9)
- Evidence: no data vector consistent with the printed slope (0.2760) and
  R squared (78.82 percent) can produce an intercept near 293 with the
  known x values; the digitized markers give 240.9. Slope and R squared
  are reproducible, the intercept is not. The reconciliation in
  inputs/digitize_figures.py drops this constraint (recorded in
  inputs/ground_truth.json under reconciliation_dropped_constraints).

## C3. Table 6.2, Criterion IV explained variance

- Printed: R squared = 00.14 percent
- Corrected (operative): about 1.37 percent
- Evidence: the printed slope (0.0517) and intercept (2202) reproduce from
  the digitized markers, and that same data gives R squared 1.37 percent.
  0.14 percent is not obtainable from any data consistent with the other
  eleven reported statistics. Materially irrelevant (both values mean "no
  correlation") but corrected for the record.

## C4. Table 6.1a, Room 3 threshold values

- Printed: Thr II 1.9927, Thr III 1.9605, Thr IV 1.9284
- Corrected (operative): 1.9929, 1.9608, 1.9287
- Evidence: the printed mean (2.0250) and SD (0.0321) imply mean minus
  1, 2, 3 SD of 1.9929, 1.9608, 1.9287. The printed thresholds carry a
  uniform 0.0003 deficit consistent with an unrounded mean near 2.0247
  that was rounded up for display. The corrected values are what the
  stated definition produces from the stated inputs. Applied in
  inputs/thesis_params.md Section 4.1. No computation consumes these
  values (they are comparison-only), so no results change.

## C5. Running text, Room 2 Criterion I detected sample

- Printed (text, p.43): 4571 (103.70 ms)
- Corrected (operative): 4573 (Table 6.1b), since 4573 / 44100 = 103.70 ms
  is the self-consistent pair and 4571 is a text typo.
- Applied throughout: all analysis uses 4573.

## C6. Abstract, claimed tmp50% range

- Printed: "mean perceptual mixing times (tmp50%) ranging from 37.5 ms to
  97.8 ms"
- Corrected (operative): the thesis's own data (Figure 6.2 markers) range
  from about 22.5 ms (Room 3) to 97.8 ms (Room 7). The 37.5 ms minimum is
  Room 1, not the minimum of the set.
- Evidence: digitized markers, inputs/ground_truth.json. The ground truth
  used for regression already reflects the corrected values.

## Non-corrections (looked wrong, left alone, for the record)

- GUI screenshot swaps source and receiver y relative to the text. By
  acoustic reciprocity the RIR is identical, so there is nothing to
  correct; both are recorded in thesis_params.md.
- Figure 5.1b shows a Room 1 RIR plotted to 0.2 s while the profile in
  Figure 5.3 spans about 0.39 s. Treated as a plotting window, not a data
  error; the RT times fs length assumption is documented and sensitivity
  tested instead.
