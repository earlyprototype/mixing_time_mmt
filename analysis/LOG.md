# Running log of subagent results

Audit trail of what each parallel subagent was asked to do and what it
returned. Filled in as results arrive.

## Round 1 (parallel, launched together)

### Agent A: thesis parameter extraction and ground truth digitization
Status: running. Deliverables: inputs/thesis_params.md,
inputs/ground_truth.json, inputs/digitize_figures.py.

### Agent B: HFD pipeline and unit tests
Status: DONE. Delivered src/mixtime/{higuchi,profile,detect}.py matching the
thesis MATLAB semantics line by line, plus tests (22 passed, re-verified in
the main thread). Profile length convention resolved from HMW.m: profile has
len(rir) - 2 samples, value s is HFD of the 51-sample zero-padded segment at
s. Validation vs known FD: white noise 2.0026 (expect 2), sine 1.0003
(expect 1), fBm H=0.2/0.5/0.8 gave 1.8008/1.4975/1.1780 (expect
1.8/1.5/1.2, kmax=10), Weierstrass D=1.5 gave 1.4894. Fast profile matches
reference to 2.7e-15; 100k-sample profile in ~0.02 s.

### Agent C: room reconstruction and ISM generator
Status: running. Deliverables: src/mixtime/{rooms,ism,rt_check}.py,
inputs/rooms_reconstruction.md, Room 1 end-to-end RT check.

### Agent C result (room reconstruction and ISM generator)
Status: DONE. Rooms digitized from Figure 4.3.1.1 floor plans; Sabine RT
from reconstructed geometry matches tabulated RT within 1 percent for rooms
1 to 8 (untuned), validating the proportions. Room 1 end to end: RT60 0.454 s
vs 0.39 tabulated (expected non-diffuse ISM excess).

## Round 2 (per-room fan-out, analysis/run_room.py, both RIR tracks)

Main thread first ran Room 1 and found the pyroomacoustics profile saturates
by sample ~1070 (dense sinc kernels), while a faithful Allen-Berkley port
lands near the thesis crossings. Pre-registration Amendment 1 added the AB
track before any regression. See results/rooms/room1.json.

### Batch agent rooms 2, 3, 4: DONE
Room 2 ab crossings 4405/4223/1829/1680 (thesis 4573/3197/3170/3049),
room 3 ab 3498/2900/2878/1815 (thesis 6056/4335/2494/1567),
room 4 ab 10650/8338/5237/5222 (thesis 13724/7027/7013/7004).
pra track floors at the search start (1000) for most k, confirming the
Room 1 saturation pattern. Schroeder RT60 overshoots tabulated RT by up to
about 1.3x for the larger rooms, both tracks.

### Batch agent rooms 5, 6, 9: rooms 5 delivered by agent; the container was
recycled mid-run (3 day gap), killing rooms 6 and 9; both rerun to
completion in the main thread with the identical deterministic script.
### Batch agent rooms 7, 8: room 7 delivered by agent; room 8 killed by the
same container recycle and rerun in the main thread.

## Round 3 (main thread): analysis/analyze.py over rooms 1 to 8, both
tracks. Outputs results/analysis.json, results/summary.md, figures.
Headline: Criterion III R2 47.2 percent (pra) / 57.1 percent (ab) vs 93.49
claimed; LOOCV 3.5 / 23.3 percent; bootstrap CIs span [4, 91] percent.
Best-of-sweep criterion flips from III to I. See REPORT.md.

## Round 4 (author-requested exploratory, two parallel agents)

### Agent D: AB vs pra diff characterization: DONE
Same echo set, different rendering: per-room peak counts match nearly
one-for-one; AB RIRs have 129 leading zeros and up to 43 percent
near-zero samples, pra RIRs have none anywhere. pra profiles reach 50
percent rise within about 40 samples of their minimum in 6 of 8 rooms;
AB takes 500 to 5200. Cross-track: ab = 1.58 pra + 1718, R2 0.83,
Spearman 0.93.

### Agent E: floor-free detector variants: DONE
Min-referenced start, rise-fraction (q 0.5/0.7/0.9) and peak-train
re-sparsification all collapse on pra (R2 0.002 to 0.009, negative
LOOCV, Spearman near or below zero), proving the pra baseline and the
start=2000 sensitivity result were floor artifacts. AB q=0.9 variant:
R2 0.652, LOOCV 0.385, exploratory only. All variants reported, none
endorsed. Full tables in results/exploratory/detectors.md.

### Agent A (ground truth digitization): DONE
tmp50 recovered from Figure 6.2 at native resolution, four panels agreeing
within 0.5 samples; independent check against Figure 4.3.2; reconciled
against the four reported regressions. Reconciled values (samples): 1654.8,
1324.0, 992.3, 2645.7, 2326.4, 3087.4, 4310.8, 2867.3 for rooms 1 to 8.
All markers sit within 3.5 samples of round-ms values (37.5, 30.0, 22.5,
60.0, 52.8, 70.0, 97.8, 65.0), almost certainly the exact 2011 vector.
Two thesis typos identified (Criterion I slope in Table 6.2 is 0.2090 vs
figure 0.2019 which reproduces; Criterion II intercept 293.2 not
reproducible, self-consistent value about 241). Criterion IV R2 is
actually about 1.37 percent, not the reported 0.14. Room 9 tmp50 and all
tmp95 marked missing. Notable: the thesis abstract's own minimum (37.5 ms)
is wrong for its data, the measured minimum is Room 3 at about 22.5 ms.

## Round 5 (reverse engineering, author-requested)

Search harness analysis/reverse/search_room.py; four agents (rooms 1, 3,
5, 7) plus a main-thread background rerun for the stalled rooms 2, 4, 6,
8. All eight rooms: no exact fingerprint match; tail statistics match to
the fourth decimal everywhere; integer crossings never reproduced; beta
= sqrt(1 - alpha) confirmed in all rooms; best-fit lengths 0.55 to 0.85
of RT*fs. Synthesis in results/reverse/SYNTHESIS.md, folded into
REPORT.md.

## Round 6 (follow-up groundwork, after PR #1 merged)

Branch restarted from merged main. Committed BEFORE any groundwork ran:
analysis/PREREGISTRATION_FOLLOWUP.md (window-scale follow-up, frozen W*
procedure, floor-free detection, LOOCV primary) and
docs/ARCHIVE_CHECKLIST.md (ranked list for the author's archive hunt).

Window sweep on ISM data (hypothesis-generating): W* procedure returned
None. No window rescues the pyroomacoustics track under floor-free
detection (LOOCV deeply negative everywhere); the AB track degrades sharply above W=50 and stays near zero; the earlier W=100 74.7 percent point was a
fixed-start floor artifact. results/exploratory/window_sweep.md.

Measured-BRIR feasibility (IoSR RealRoomBRIRs, Rooms A to D, 48 kHz,
fetched via raw.githubusercontent, not committed, provenance recorded):
H1 PASSES at every window in all four rooms, rises 12 to 194 times the
tail SD. Real measured rooms have true pre-arrival silence and sparse
strong early reflections, so their profiles rise like the AB track, not
like the dense-from-t0 pyroomacoustics synthesis. The instrument
deflects on real data; whether it points at perceptual mixing time is
exactly what H2 (labeled data) must decide.
results/exploratory/measured_feasibility.md.

## Round 7 (kernel track, author-suggested)

analysis/exploratory/kernel_track.py: third rendering, AB spike train
convolved with a causal measurement-like bandpass pulse (100 Hz to 16
kHz, 90 percent of energy within 0.09 ms), replacing both the bare
spikes and the ideal sinc. Result: the profile RISE survives in all
eight rooms (12.8 to 24.3 times the tail SD), confirming the texture
reasoning, but the TIMING information at window 50 does not: crossings
compress to 852 to 2803 samples, the floor-free detections are nearly
constant (852 for five rooms), and the honest floor-free regression is
useless (R2 0.39 with LOOCV -0.94 and negative Spearman). The
fixed-start numbers (R2 0.68) are partially floor-inflated again. Even
a 0.09 ms kernel is enough to accelerate window filling and destroy
most of the room-dependence the 51-sample feature reads from bare
spikes. results/exploratory/kernel_track.md.

## Round 8 (CodeRabbit review of PR 2, amended reruns)

Review round on the phase 2 PR: 8 actionable comments, all addressed.
Substantive: the pre-registration's "after the direct sound" onset was
under-specified and unimplemented (some floor-free detections landed at
sample 1), fixed by Amendment 1 (onset = first sample at 5 percent of
peak absolute amplitude, argmin restricted to at or after onset) and
all three groundwork scripts re-run under the amended rule. Verdicts
unchanged: W* still None, H1 still passes 4 of 4 rooms at every window,
and the kernel track's floor-free result is starker (rooms 1 to 3
detect nothing, rooms 4 to 8 detect the identical sample 852, zero
between-room information). Also fixed: kernel order labeling (4th order
prototype, 8th order realized), full convolution tail retained, atomic
checksummed downloads, figures directory creation, frontal azimuth
tolerance assert, order-independent ground truth lookup, per-cell n
columns in the sweep table, W* metric choice justified in the
amendment, markdown and gitignore nits.
