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

### Agent C: room reconstruction and ISM generator
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

### Batch agent rooms 5, 6, 9: running
### Batch agent rooms 7, 8: running
### Agent A (ground truth digitization): running
