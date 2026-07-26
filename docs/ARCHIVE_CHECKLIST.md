# Archive hunt checklist: what to look for and what each item resolves

For the author's search of the 2011 project materials. Ranked by value.
Context: the reverse-engineering study (results/reverse/SYNTHESIS.md)
showed the printed thesis reproduces the pipeline's plateau behavior to
four decimals but cannot regenerate the Table 6.1b detections, so the
missing information is in the artifacts below, not in the document.

## 1. Saved RIRs: the jackpot

Files like room1.mat, rir*.mat, *.wav, or workspace .mat dumps
containing the eight ISM outputs.

Resolves EVERYTHING at once: push them through the rebuilt pipeline
(one command) and either Table 6.1b regenerates, closing every open
question, or the exact processing difference becomes directly
observable by comparing the stored waveform against our AB port for
the same geometry. Even one room's RIR is decisive.

## 2. The GUI callback file: sroomguicbs.m

The thesis appendix prints sroomgui.m (layout) but NOT sroomguicbs.m,
the callbacks it references, which is where the RIR is actually
generated, possibly normalized, trimmed, padded, or resampled before
saving. This file is the prime suspect for the un-printed processing
step the reverse engineering detected. Also valuable: any edited copies
of sroom.m / lthimage.m / HMW.m / hfd.m that differ from the printed
versions.

## 3. Original MATLAB figure files (.fig)

Figures 5.1b, 5.2, 5.3, 6.2 as .fig files. MATLAB .fig files embed the
plotted data arrays: Figure 5.3's .fig contains the entire Room 1 FD
profile, and Figure 6.2's contain the exact tmp50 and detection values.
Extraction is trivial (a .fig is a .mat). The Room 1 profile alone
would let us binary-search the processing difference.

## 4. Anything recording the typed inputs

Lab notebook, spreadsheet, script, or screenshots recording: the room
dimensions entered per room (L, W, H), the Response length (sec) typed
per room, and source and receiver coordinates. The reverse engineering
suggests lengths were hand-set per room (best fits 0.55 to 0.85 of
RT times fs) and could not pin them down.

## 5. Environment details

MATLAB version (smooth() and butter() behavior), operating system, and
the exact Curve Fitting Toolbox version used for the regressions.
Lowest priority, only matters if items 1 to 3 surface and a residual
discrepancy remains.

## What NOT to worry about

- The tmp50 ground truth: already recovered to round-millisecond
  precision from Figure 6.2; nothing in the archive can improve it.
- The alpha to beta conversion: confirmed as sqrt(1 - alpha) by the
  reverse engineering, all eight rooms.
- The Higuchi / windowing / threshold code semantics: verified to four
  decimals against Table 6.1a.

## Handling note

Per analysis/PREREGISTRATION_FOLLOWUP.md, the follow-up analysis
constants are frozen before any archive material is examined. When
material is found, commit it (or its extracted numbers) to the repo
BEFORE running anything against it, so the record shows what was known
when.
