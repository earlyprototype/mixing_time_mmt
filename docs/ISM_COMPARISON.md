# Why two image-source models produce such different RIRs, and why it matters here

Companion note to REPORT.md, written for the thesis author. The measured
numbers in Section 5 come from the exploratory diff analysis in
results/exploratory/ (clearly labeled exploratory, n = 8).

## 1. What the two implementations agree on

Everything geometric. Both start from Allen and Berkley's 1979 insight:
a rectangular room's reflections are exactly equivalent to a lattice of
mirror-image sources in an infinite grid of mirrored rooms. Both compute
the same image positions, the same arrival delays (distance over c), and
the same per-arrival attenuations (wall reflection losses times spherical
spreading). Given the same room, the two produce essentially the same
LIST of arrivals: when each echo lands and how strong it is. The Spearman
correlation between the two tracks' detected crossings (0.93) reflects
this shared geometry.

## 2. Where they part ways: rendering the list into a sampled signal

An arrival time is a real number, say 5.2837 ms = 233.01 samples at
44.1 kHz. A sampled signal has no sample 233.01. Each implementation must
decide what to write into the sample grid, and this is the entire
difference.

The Allen-Berkley Fortran code (and the thesis's MATLAB port of it)
ROUNDS: the arrival goes into sample 233 as a single-sample spike of the
arrival's amplitude. The result is a sparse train of isolated deltas with
true digital silence between them, then a butter(3, 0.01) highpass. This
was the standard and computationally sane choice in 1979, and still the
natural one in a 2011 MATLAB port.

pyroomacoustics INTERPOLATES: the arrival is written as a windowed sinc
kernel (81 taps, about 1.8 ms wide) centred at the fractional position
233.01. This is the textbook-correct band-limited representation: a
sampled system cannot actually contain an instantaneous click, and a
real measurement chain (loudspeaker, air, microphone, anti-alias filter)
would spread every arrival in exactly this way. Thousands of overlapping
kernels sum into a dense, noise-like waveform from the first
milliseconds.

So: identical physics, identical echo list, two legitimate answers to
"what does this look like at 44.1 kHz". The AB rendering is faithful to
the MATHEMATICAL idealization (a train of Dirac impulses). The
pyroomacoustics rendering is faithful to what a MEASUREMENT would
record.

## 3. Why the 2011 choice of AB was reasonable

- It is the published reference implementation; the thesis ported the
  code printed in the original JASA paper, which is the most defensible
  reproduction of "the image source method" circa 2011.
- Fractional-delay ISM renderers were not packaged and accessible then;
  pyroomacoustics dates from 2017.
- For most uses of an ISM RIR (reverberation time, energy decay, echo
  density at coarse scales, convolution reverb), the rendering choice is
  a second-order detail. Both tracks here produce nearly identical
  Schroeder decay curves and RT60s. Nothing about the 2011 project
  flagged rendering as a live risk.

The catch was not foreseeable without a modern renderer to compare
against: the chosen FEATURE happens to live at exactly the scale where
the two renderings differ.

## 4. Why the feature reacts so violently to the rendering

The moving-window Higuchi FD is computed over 51-sample windows, about
1.16 ms. That is SHORTER than one pyroomacoustics sinc kernel (1.8 ms).

- On the AB rendering, an early window contains a few isolated spikes
  separated by silence: a curve that fills space poorly, FD near 1. As
  echo density builds, windows fill up and FD climbs toward 2. The
  profile's rise IS the physics the thesis wanted: growing echo density.
- On the pyroomacoustics rendering, every window from the first
  milliseconds already contains overlapping smooth oscillating kernels:
  a dense curve, FD near 2 immediately. The profile starts near its
  plateau, and the detector, searching from sample 1000 for the profile
  to reach the plateau region, fires almost at once. Detections pin at
  the search floor for the small rooms.

In one sentence: the HFD profile is a millisecond-scale
spike-sparsity meter, the AB rendering preserves spike sparsity, the
band-limited rendering erases it, and the mixing-time information rides
on the sparsity.

CORRECTION (added after the measured-BRIR feasibility check,
results/exploratory/measured_feasibility.md): the inference that
measured BRIRs would behave like the pyroomacoustics track because both
are band-limited turned out to be WRONG. Real rooms have true
pre-arrival silence and sparse, strong early reflections, so measured
HFD profiles exhibit a large rise (12 to 194 times the tail SD on the
IoSR rooms at every window tested), behaving like the AB track in the
early region rather than the dense-from-the-first-sample synthesis.
pyroomacoustics' early density is a property of its synthetic onset,
not of band-limiting per se. Whether the measured-data rise lands at
perceptually meaningful times remains untested (H2 in
analysis/PREREGISTRATION_FOLLOWUP.md).

The original (superseded) claim continued: The baseline detections do rank the rooms plausibly
(Spearman 0.86 against tmp50), but Section 6 shows that ranking is
mostly an artifact of how fast each room's profile saturates past the
fixed search start, not a readable mixing-time signal. A precise
wording of the pathology: detections do not pin exactly at sample 1000
(only Room 2 does), they cluster just above it, rooms 1 to 4 at 1000 to
1091, rooms 5 to 8 at 1542 to 3830.

## 5. Measured differences (from the exploratory diff analysis)

Full tables in results/exploratory/profile_diff.md, plots in
results/figures/exploratory/ (rir_zoom_room1.png is the clearest single
picture). Headlines:

- Same echo set, different rendering: per-room peak counts match almost
  one-for-one between tracks. Every AB RIR has exactly 129 leading zero
  samples and up to 43 percent near-zero samples between spikes; the
  pyroomacoustics RIRs contain not a single exactly-zero or near-zero
  sample anywhere.
- Profile shape: the pyroomacoustics profile recovers from its minimum
  to 50 percent of its rise within about 40 samples in six of eight
  rooms (about 1 ms, i.e. instantly at this scale). The AB profile takes
  500 to 5200 samples, and that slow rise is the quantity the detector
  actually measures. Tail plateaus agree closely between tracks.
- Cross-track detections: ab = 1.58 * pra + 1718, R2 0.83, Spearman
  0.93, so the two renderings rank the rooms almost identically while
  disagreeing on absolute values.

## 6. Can the pyroomacoustics output be rescued? Tested, and no

The floor-free detector study (results/exploratory/detectors.md, all
variants reported, none endorsed, n = 8, post-hoc) is blunt:

- Min-referenced start, rise-fraction detection (q = 0.5, 0.7, 0.9), and
  peak-train re-sparsification ALL collapse on the pyroomacoustics
  track: R2 0.002 to 0.009, Spearman near zero, negative LOOCV.
- This retroactively explains the seemingly promising start=2000 result
  (R2 80.7 percent) in the sensitivity table: with three rooms clamped
  exactly at the start value, the floor constant itself was carrying the
  fit. Remove the floor and there is nothing underneath on this track.
  The baseline Spearman of 0.86 was likewise largely floor-ordering.
- On the AB track the rise-fraction q = 0.9 variant is the only one that
  edges past the original rule (R2 0.652, LOOCV 0.385 vs 0.571 / 0.233),
  a modest, exploratory-only observation.

Conclusion for the rendering question: the mixing-time information the
2011 method reads is carried by the spike-sparsity texture of the
signal. Once that texture is gone, as in the pyroomacoustics synthesis,
no affine rescale, re-referencing, or re-sparsification of the same HFD
feature recovers it. The later measured-BRIR check (see the correction
in Section 4) showed real rooms DO retain early sparsity, so the
follow-up question is not whether the feature deflects on measured data
(it does) but whether its crossing times track perception, which is
exactly what the pre-registered H2 in
analysis/PREREGISTRATION_FOLLOWUP.md will decide when labeled data
exists.
