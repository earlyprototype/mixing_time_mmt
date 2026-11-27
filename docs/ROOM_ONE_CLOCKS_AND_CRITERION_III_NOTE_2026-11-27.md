# Room 1's two clocks and the missing share of the Criterion III fit

*Written on 27 November 2026 for the General, this note sits beside the repository's replication report. What did Pathfinder Scarab's “two clocks” measure in the studio known as Room 1? How much of the thesis's agreement with listeners survives when the room responses are rebuilt? Those questions determine what to look for in the author's archive.*

> **Provenance.** On 2026-11-27 I read Scarab's 2026-11-26 field report, Conaty's thesis tables and printed simulation code, the two reconstruction routines and the saved room results at repository revision `fec6cee`. I computed the first arrival twice in memory for Room 1, measured the first visible peak on the published plot, and refitted straight lines to the stored room values. The listener observations and full-room simulations came from earlier work. I checked recalled details against the sources below. [established] marks a source or check named here, [inferred] a conclusion drawn from it, and [speculation] a possible explanation still to test.

## 1. The answers in brief

**The same one-metre path takes 129 samples under the reconstruction's air-speed setting and 44 under the printed interface's metre conversion.** [established] At 44,100 samples per second, the 85-sample difference puts the first sound about 1.93 milliseconds earlier. [established] The thesis plot favours the later arrival, but the original Room 1 settings are still unknown. [inferred]

**The published detections recover the thesis's 93.49% fit to listening times across eight rooms.** [established] With the same listening times, the rebuilt detections explain 47.22% of the variation for pyroomacoustics and 57.11% for Allen-Berkley, shortfalls of **46.27 and 36.38 percentage points**, respectively. [established] A timing difference in one room cannot tell us how much of either shortfall it caused. [inferred]

## 2. Where the 85 samples come from

Room 1 is E. Studio. The thesis puts the sound source and microphone at the same height and the same position along the room's length, respectively two metres and one metre from the same side wall. [established] The first sound travels the one metre between them before any wall reflection arrives. [established]

Both “clocks” count 44,100 samples, or equally spaced sound values, per second from emission. They differ in the time allowed for sound to travel one metre. [established] The Allen-Berkley reconstruction converts distance to samples by multiplying metres by 44,100 and dividing by 343, its assumed speed of sound in metres per second. [established] The printed MATLAB interface multiplies distance by the sample rate and a unit factor, then divides by 1,000. With the factor set to one and distances entered in metres, that amounts to 1,000 metres per second. [established]

| Conversion for the one-metre path | Delay before rounding, in samples | First nonzero sample, counted from zero at emission |
| --- | ---: | ---: |
| Allen-Berkley default | 44,100 / 343 = 128.57 | 129 (2.93 milliseconds) [established] |
| Printed interface, assuming metres and factor one | 44,100 / 1,000 = 44.10 | 44 (1.00 millisecond) [established] |

I changed only the speed argument in the same Allen-Berkley routine and inspected the first nonzero value in each short response. [established] Subtracting the rounded delays gives 129 minus 44, or 85 samples. The printed MATLAB code adds one to store each delay because its arrays start at one, so the matching storage positions would be 130 and 45. [established]

Figure 5.1a, the example setup screen, shows “metres”. [established] In Figure 5.1b, the Room 1 response beneath it, the first tall peak sits about 20 pixels along a 1,274-pixel axis spanning 200 milliseconds. That puts it about 3.1 milliseconds after zero, closer to 2.93 than to 1.00 milliseconds. [established] The drawing supports the later arrival to the precision of the plot. Which interface settings produced the results table is still open. [inferred]

## 3. What the eight-room fit measures

Conaty's 2011 Trinity College Dublin thesis fitted eight of Lindau's nine rooms. The ninth, a church, fell outside the rectangular room model. [established] Each room supplied a listening time, the time after which listeners had difficulty distinguishing a change to the direction of the late echoes, and a detection time calculated from a simulated echo trace. [established]

Criterion III is the thesis's third detection rule. It scans a smoothed Higuchi fractal dimension profile, a running estimate of how rough the echo trace is, starting at sample 1,000 (about 22.7 milliseconds from the start). [established] It records the first sample at or above the average roughness of the final tenth of the profile minus twice that portion's standard deviation, a measure of its spread. [established] In Room 1 this happens well after the first sound under all three sets of detections. [established]

For each row below I fitted a straight line from detection time to listening time, allowing the line to cross the vertical axis above or below zero. The score, R², measures the share of variation in the eight listening times explained by that line. Predicting their mean for every room scores 0%, and matching every room exactly scores 100%. [established]

| Source of detection times | Room 1 detection (sample number, counting from one) | R² across the same eight rooms | Shortfall from the printed 93.49% |
| --- | ---: | ---: | ---: |
| Thesis Table 6.1b, the published crossing times | 3,136 | 93.49% | 0.00 percentage points [established] |
| pyroomacoustics, which spreads each modelled arrival across nearby samples | 1,070 | 47.22% | 46.27 percentage points [established] |
| Allen-Berkley, which places each arrival on a whole sample before filtering | 4,032 | 57.11% | 36.38 percentage points [established] |

Table 6.2, the thesis's summary of the fitted lines, gives 93.49% for the third rule. The listening times used in that table survive as plotted points in Figure 6.2. [established] The repository read their heights from the figure, then adjusted them slightly to agree with the printed fits. Using the unadjusted heights with the published detections gives 93.52%, against 93.49% after adjustment. [established] The close agreement checks the transcription. [inferred] When the repository withheld each room in turn and predicted its listening time from the other seven, the scores fell to 3.52% for pyroomacoustics and 23.28% for Allen-Berkley, against the same eight-room mean. [established]

The printed line adds 0.3197 samples of listening time for each extra sample of detection time, starting from 325 samples at a detection time of zero. [established] Its slope resembles 0.343, the ratio of the two assumed speeds, but a slope alone cannot identify a speed setting. Adding the same delay to every detection shifts where the fitted line crosses the vertical axis. Multiplying every detection by the same positive factor changes the slope. In either case each room keeps the same fitted listening time and the same share of variation is explained. [established]

Changing the simulated speed also moves the echoes relative to the fixed span of samples used to measure roughness. That can change which echo first crosses the threshold in each room. [inferred] This could account for part of the shortfall. [speculation]

## 4. What remains, and what needs the operator's decision

**What happened.** The Room 1 check reproduced Scarab's first arrivals. The eight-room check recovered the printed score only when it used the thesis's detections. [established]

**What it means.** The unit conversion gives us a specific setting to investigate. The aggregate scores alone cannot identify what the author ran. [inferred]

**What remains.** The following costs describe the proposed work in archive searches and model runs.

1. Search the author's archive once for the saved Room 1 waveform, the settings used to generate it and the original plotted data. The appendix already contains `sroomguicbs.m`, the interface routine that applies the unit conversion. [established] A saved response with a known time origin and the copy of the routine used that day would let us check both the first arrival and the later detection against the published record. [inferred]
2. If authorised, run one additional Allen-Berkley reconstruction for each of the eight fitted rooms using the printed metre conversion. Hold the geometry, response length and detection rule fixed. Compare first arrivals, crossings and the eight-room fit against the existing reconstruction and the published values. [inferred]

**What needs the operator's decision.** Close the replication with the discrepancy recorded, or authorise the archive search and the eight-room comparison. I recommend the archive search first. Even if a new simulation matches the printed score, only a contemporary record can tie that setting to the original run. [inferred]

## Sources

- Pathfinder Scarab, *DrJones Recon Brief for earlyprototype/mixing_time_mmt*, field report of 26 November 2026, “Test the thesis GUI’s delay conversion as a separate reconstruction”.
- Thom Conaty, *The Fractal Dimension of Room Impulse Responses as a Predictor for the Perceptual Mixing Time*, Trinity College Dublin, 2011. Printed pp. 32, 38 to 46, 53, 57 to 58 and 70 cover the listening task, room setup, detection rule, results, integer delays and interface conversion. [Thesis](../inputs/FINAL_thom_conaty.pdf) · [Room 1 figures](../inputs/figures/page51.png) · [Regression figure and table](../inputs/figures/page57.png)
- Implementation and stored inputs. [Allen-Berkley simulator](../src/mixtime/ab_ism.py) · [pyroomacoustics simulator](../src/mixtime/ism.py) · [Recovered listening times and tabulated crossings](../inputs/ground_truth.json) · [Room 1 detections](../results/rooms/room1.json) · [Eight-room detection vectors](../results/analysis.json) · [Fit calculation](../analysis/analyze.py)
- Findings and archive leads. [Replication report](../REPORT.md) · [All eight rooms' fitted results](../results/summary.md) · [Thesis review](DeepReviewReport.md) · [Archive checklist](ARCHIVE_CHECKLIST.md)
