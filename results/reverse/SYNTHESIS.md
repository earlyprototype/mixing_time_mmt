# Reverse-engineering synthesis: can the 2011 inputs be recovered by deduction?

EXPLORATORY DIAGNOSTIC. Eight per-room searches (about 6000 pipeline
evaluations each plus local refinement) over footprint aspect, height,
RIR length, beta convention and axis interpretation, scored against the
thesis's printed fingerprints: Table 6.1a tail mean and SD (4 decimals)
and Table 6.1b crossings (exact integers). Per-room detail in
room{N}.json and room{N}.log.

## Result: not invertible. Zero exact matches in eight rooms.

| Room | score | tail match | crossings match | best npts / RT*fs |
|---|---|---|---|---|
| 1 | 0.963 | 4 decimals | no (misses 100 to 1700) | 0.70 |
| 2 | 0.979 | 3 to 4 dec | no | 0.55 |
| 3 | 0.484 | 4 decimals | closest (misses 60 to 620) | 0.85 |
| 4 | 0.597 | 4 decimals | no | 0.85 |
| 5 | 0.751 | 3 to 4 dec | no (all four high) | 0.80 |
| 6 | 0.543 | 4 decimals | no | 0.80 |
| 7 | 1.057 | 4 decimals | no (compressed spread) | 0.80 |
| 8 | 0.877 | 4 decimals | no | 0.55 |

## What the failure pattern proves

1. **The plateau chain is right.** In every room the search reproduces
   the printed tail mean and SD to the fourth decimal. The Higuchi
   port, windowing, smoothing and tail statistics behave like the 2011
   chain at the plateau. This is strong evidence the near-misses are
   not a bug in our port of the printed appendix code.
2. **The crossing region depends on something the document does not
   contain.** No configuration in a broad, refined search lands the
   integer crossings, and the misses are structured, not random: our
   best fits often compress the spread between Criteria I to IV where
   the thesis values are widely spread (Rooms 5, 7, 8), meaning the
   2011 profiles rose through the mid-region with a different shape
   than any configuration of our forward model produces.
3. **The optimizer signals model mismatch, not parameter mismatch.**
   Where crossings pull hardest, the search contorts geometry to
   implausible shapes (Room 4 best fit is 47 x 29 x 2.4 m for a lecture
   hall) without reaching the targets. When the best fit inside a
   parameter space is both implausible and still wrong, the discrepancy
   lives outside the parameter space.
4. **Two secondary deductions did land.** beta = sqrt(1 - alpha) beats
   beta = 1 - alpha in all eight rooms (the assumption is now
   evidence-backed), and best-fit RIR lengths are consistently SHORTER
   than RT times fs (0.55 to 0.85 of it), suggesting the 2011 response
   lengths were set by hand per room rather than derived from RT,
   consistent with the GUI's free-text "Response length (sec)" field.
   The per-room values are not consistent enough to identify.

## Candidate locations of the un-printed difference

Not decidable from the document; listed for the record: the GUI
wrapper's handling of the RIR before HMW.m (normalization, trimming,
padding), MATLAB smooth() NaN propagation over the pre-arrival windows
(our port averages around NaNs), exact per-room response lengths typed
into the GUI, and exact room dimensions entered (possibly rounded
values not derivable from the floor plans).

## Implication for the verdict

Unchanged, and sharpened: a direct search aimed at the thesis's own
numbers, with the thesis's own algorithm, cannot find ANY input
configuration that reproduces them. The 2011 detections are therefore
not recoverable from the thesis alone, which upgrades "we could not
rebuild the inputs" from a limitation of this project to a demonstrated
property of the document. The tail-statistic agreement simultaneously
shows the method's plateau behavior IS reproducible; only the
mixing-time-bearing transition region is not.
