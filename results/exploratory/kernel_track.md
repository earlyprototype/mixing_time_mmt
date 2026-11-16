# AB + measurement-like kernel track (EXPLORATORY)

Kernel: causal Butterworth bandpass pulse, 100 Hz to 16 kHz, 90 percent of energy within 0.09 ms.

| room | rise/tailSD | near-zero frac (first 2000) | Crit III (fixed start) | k=2 (floor-free) |
|---|---|---|---|---|
| 1 | 4.3 | 0.13 | 1011 | None |
| 2 | 5.0 | 0.18 | 1000 | None |
| 3 | 5.5 | 0.12 | 1010 | None |
| 4 | 7.2 | 0.57 | 1000 | 852 |
| 5 | 8.8 | 0.62 | 1289 | 852 |
| 6 | 7.5 | 0.66 | 1418 | 852 |
| 7 | 8.7 | 0.75 | 2803 | 852 |
| 8 | 10.0 | 0.75 | 2002 | 852 |

| regression | n | R2 | LOOCV R2 | Spearman |
|---|---|---|---|---|
| fixedstart_k2 | 8 | 0.721 | 0.587 | 0.755 |
| floorfree_k2 | 5 | 0.000 | -0.562 | nan |
| fixedstart_I | 8 | 0.571 | 0.289 | 0.548 |
| fixedstart_II | 8 | 0.706 | 0.565 | 0.755 |
| fixedstart_IV | 8 | 0.724 | 0.592 | 0.850 |
