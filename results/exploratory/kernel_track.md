# AB + measurement-like kernel track (EXPLORATORY)

Kernel: causal Butterworth bandpass pulse, 100 Hz to 16 kHz, 90 percent of energy within 0.09 ms.

| room | rise/tailSD | near-zero frac (first 2000) | Crit III (fixed start) | k=2 (floor-free) |
|---|---|---|---|---|
| 1 | 15.8 | 0.13 | 1269 | 1269 |
| 2 | 24.3 | 0.18 | 1000 | 948 |
| 3 | 15.4 | 0.12 | 1069 | 1069 |
| 4 | 14.4 | 0.57 | 1000 | 852 |
| 5 | 17.6 | 0.62 | 1289 | 852 |
| 6 | 12.8 | 0.66 | 1418 | 852 |
| 7 | 18.4 | 0.75 | 2803 | 852 |
| 8 | 18.8 | 0.75 | 2002 | 852 |

| regression | n | R2 | LOOCV R2 | Spearman |
|---|---|---|---|---|
| fixedstart_k2 | 8 | 0.678 | 0.517 | 0.755 |
| floorfree_k2 | 8 | 0.391 | -0.944 | -0.791 |
| fixedstart_I | 8 | 0.018 | -1.865 | -0.024 |
| fixedstart_II | 8 | 0.571 | 0.287 | 0.548 |
| fixedstart_IV | 8 | 0.682 | 0.524 | 0.755 |
