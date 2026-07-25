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

This also explains the practical implication in REPORT.md: measured
BRIRs are band-limited like the pyroomacoustics track, so the method as
specified in 2011 (fixed constants, fixed start) would underperform on
real measurements for the same reason it underperforms on
pyroomacoustics. The signal is not gone (rank order survives, Spearman
0.86 against tmp50), but the fixed-constant detector cannot read it.

## 5. Measured differences (from the exploratory diff analysis)

See results/exploratory/profile_diff.md and detectors.md for the full
tables and results/figures/exploratory/ for the plots. Headline numbers
are summarized in REPORT.md's exploratory section.

## 6. What would make the method rendering-robust (hypotheses only)

The exploratory detector study (results/exploratory/detectors.md) tests
floor-free variants: referencing the search to the profile's own minimum
instead of a fixed sample, detecting a fraction of the profile's rise
instead of an absolute threshold, and re-sparsifying a band-limited RIR
into its peak train before computing HFD. All of it is post-hoc at
n = 8: it can propose a pre-registered follow-up (ideally on measured
BRIRs), it cannot confirm anything.
