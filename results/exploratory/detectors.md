# Floor-free detector variants

**EXPLORATORY, NOT CONFIRMATORY. Post-hoc detector variants at n = 8, devised after seeing the confirmatory results. No variant is selected or endorsed; all are reported identically.**

Target: tmp50_samples_reconciled (samples at 44100 Hz), rooms 1-8.
OLS is tmp50 on detected sample. LOOCV R2 is 1 - PRESS/SStot over
leave-one-room-out refits. Floor-pinned counts rooms whose
detection lies within 5 samples of that detector's own search start
(fixed sample 1000 for D0 and D3, the profile argmin for D1 and D2).
Rooms with no crossing are null and excluded from the fit (reported n).
D3 is pra-only by design. All detections are 1-based sample numbers.

| detector | track | n | floor-pinned | slope | intercept | R2 | Spearman rho | LOOCV R2 |
|---|---|---|---|---|---|---|---|---|
| D0 | pra | 8 | 1 | 0.5611 | 1228.3 | 0.472 | 0.857 | 0.035 |
| D0 | ab | 8 | 0 | 0.3546 | 620.1 | 0.571 | 0.714 | 0.233 |
| D1 | pra | 8 | 0 | 0.0436 | 2354.5 | 0.003 | 0.167 | -0.453 |
| D1 | ab | 8 | 0 | 0.3546 | 620.1 | 0.571 | 0.714 | 0.233 |
| D2q50 | pra | 8 | 0 | 0.0336 | 2373.4 | 0.002 | -0.528 | -0.454 |
| D2q50 | ab | 8 | 0 | 0.4325 | 1150.7 | 0.500 | 0.667 | 0.115 |
| D2q70 | pra | 8 | 0 | 0.0314 | 2374.4 | 0.002 | -0.503 | -0.454 |
| D2q70 | ab | 8 | 0 | 0.2386 | 1568.0 | 0.296 | 0.714 | -0.143 |
| D2q90 | pra | 8 | 0 | 0.0452 | 2352.8 | 0.003 | 0.228 | -0.453 |
| D2q90 | ab | 8 | 0 | 0.3595 | 603.7 | 0.652 | 0.810 | 0.385 |
| D3 | pra | 8 | 0 | 0.1470 | 2137.3 | 0.009 | 0.000 | -1.581 |

Full per-room detail (detections, search starts, thresholds, tail
statistics) is in results/exploratory/detectors.json. Figure:
results/figures/exploratory/detectors_scatter.png.

*EXPLORATORY, NOT CONFIRMATORY. Post-hoc detector variants at n = 8, devised after seeing the confirmatory results. No variant is selected or endorsed; all are reported identically.*
