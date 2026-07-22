# Re-validation of Conaty 2011: mixing time from Higuchi fractal dimension profiles

This repository re-implements and honestly re-validates the method of Conaty
(2011, MPhil thesis, Trinity College Dublin): detecting the perceptual mixing
time of a room from the moving-window Higuchi fractal dimension profile of an
image-source-model room impulse response, regressed against the perceptual
mixing times (tmp50 percent) measured by Lindau et al. (2010) for nine rooms.
The thesis reported R squared = 93.49 percent for its Criterion III.

## Layout

- `inputs/` thesis PDF, extracted text and figures, extracted parameters,
  digitized ground truth with provenance
- `src/mixtime/` pipeline: Higuchi FD, moving-window profile, threshold
  detection, room reconstruction, ISM RIR generation
- `tests/` unit tests, including HFD validation against signals of known
  fractal dimension
- `analysis/` pre-registration (committed before results), analysis scripts
- `data/rirs/` generated RIRs and per-room profiles
- `results/` regression outputs, sweep curves, figures
- `REPORT.md` the verdict

## Reproducibility

Python 3, package versions pinned in `requirements.txt`. All random seeds are
fixed in the scripts that use them. Run `pytest` for the unit tests, then the
scripts in `analysis/` in numbered order.
