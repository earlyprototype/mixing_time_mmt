"""Window-scale sweep on the existing ISM RIRs. HYPOTHESIS-GENERATING.

Runs the grid fixed in analysis/PREREGISTRATION_FOLLOWUP.md on the eight
already-generated ISM RIRs (both tracks), with both the original
fixed-start detection and the floor-free min-referenced detection, and
freezes W* by the pre-registered procedure: the smallest W whose
Criterion III (k = 2, floor-free) R squared lies within 10 points of the
grid maximum on BOTH tracks simultaneously.

Outputs: results/exploratory/window_sweep.json, window_sweep.md,
results/figures/exploratory/window_sweep.png. Deterministic.
"""

import json
import pathlib
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analysis"))

from mixtime.profile import hfd_profile
from mixtime.detect import tail_stats, detect_crossing
from run_room import moving_average_nan

ROOMS = list(range(1, 9))
WINDOWS = [50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000]
K = 2.0
TAIL = 0.1


def ols_stats(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    sstot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ((y - yhat) ** 2).sum() / sstot if sstot > 0 else np.nan
    n = x.size
    preds = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        c, _, _, _ = np.linalg.lstsq(A[m], y[m], rcond=None)
        preds[i] = A[i] @ c
    loocv = 1 - ((y - preds) ** 2).sum() / sstot if sstot > 0 else np.nan
    rho = spearmanr(x, y)[0] if n >= 3 else np.nan
    return float(coef[0]), float(coef[1]), float(r2), float(loocv), float(rho)


def detect_floorfree(smoothed, k=K, tail_frac=TAIL):
    """Min-referenced start: search from the argmin over finite samples."""
    finite = np.isfinite(smoothed)
    if finite.sum() < 100:
        return None
    idx = np.arange(smoothed.size)[finite]
    amin = idx[np.argmin(smoothed[finite])]
    return detect_crossing(smoothed, k, tail_frac=tail_frac,
                           start=int(amin) + 1)


def main():
    gt = json.loads((ROOT / "inputs/ground_truth.json").read_text())
    y = np.array([r["tmp50_samples_reconciled"] for r in gt["rooms"]
                  if r["room"] in ROOMS], float)

    out = {"note": "HYPOTHESIS-GENERATING (ISM data, not measured). "
                   "See analysis/PREREGISTRATION_FOLLOWUP.md.",
           "windows": WINDOWS, "k": K, "smooth_span_rule": "4*W",
           "cells": []}

    detections = {}
    for track, suffix in (("pra", ""), ("ab", "_ab")):
        rirs = {n: np.load(ROOT / f"data/rirs/room{n}_rir{suffix}.npy")
                for n in ROOMS}
        for W in WINDOWS:
            t0 = time.time()
            det_ff, det_fs = [], []
            rises = []
            for n in ROOMS:
                prof = hfd_profile(rirs[n], window=W, kmax=5)
                sm = moving_average_nan(prof, span=4 * W)
                det_ff.append(detect_floorfree(sm))
                det_fs.append(detect_crossing(sm, K, tail_frac=TAIL,
                                              start=1000))
                m, s = tail_stats(sm, TAIL)
                finite = sm[np.isfinite(sm)]
                rises.append(float((m - finite.min()) / s) if s > 0
                             else np.nan)
            cell = {"track": track, "W": W,
                    "elapsed_s": round(time.time() - t0, 1),
                    "detect_floorfree": det_ff,
                    "detect_fixedstart": det_fs,
                    "rise_over_tailsd": rises}
            for label, det in (("floorfree", det_ff), ("fixedstart",
                                                       det_fs)):
                ok = [i for i, v in enumerate(det) if v is not None]
                if len(ok) >= 4:
                    xs = [det[i] for i in ok]
                    sl, ic, r2, loocv, rho = ols_stats(xs, y[ok])
                    cell[label] = {"n": len(ok), "slope": sl,
                                   "intercept": ic, "r2": r2,
                                   "loocv_r2": loocv, "spearman": rho}
                else:
                    cell[label] = {"n": len(ok)}
            out["cells"].append(cell)
            detections[(track, W)] = cell
            print(f"{track} W={W}: ff R2="
                  f"{cell['floorfree'].get('r2', float('nan')):.3f} "
                  f"loocv={cell['floorfree'].get('loocv_r2', float('nan')):.3f}")

    # Freeze W* per the pre-registered procedure (floor-free rule).
    r2_by_track = {}
    for track in ("pra", "ab"):
        r2_by_track[track] = {W: detections[(track, W)]["floorfree"].get(
            "r2", np.nan) for W in WINDOWS}
    wstar = None
    for W in WINDOWS:
        ok = True
        for track in ("pra", "ab"):
            vals = np.array([v for v in r2_by_track[track].values()
                             if np.isfinite(v)])
            v = r2_by_track[track][W]
            if not np.isfinite(v) or vals.size == 0 or \
                    v < vals.max() - 0.10:
                ok = False
        if ok:
            wstar = W
            break
    out["W_star"] = wstar
    out["W_star_rule"] = ("smallest W with floor-free R2 within 0.10 of "
                          "the grid max on both tracks; frozen at first "
                          "run per PREREGISTRATION_FOLLOWUP.md")

    (ROOT / "results/exploratory").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "results/exploratory/window_sweep.json", "w") as f:
        json.dump(out, f, indent=1)

    # Figure: R2 and LOOCV vs W, both tracks, floor-free rule.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for ax, metric, title in ((axes[0], "r2", "OLS R2"),
                              (axes[1], "loocv_r2", "LOOCV R2")):
        for track, color in (("pra", "tab:blue"), ("ab", "tab:orange")):
            vals = [detections[(track, W)]["floorfree"].get(metric,
                                                            np.nan)
                    for W in WINDOWS]
            ax.plot(WINDOWS, vals, "o-", color=color, label=track, ms=4)
        if wstar:
            ax.axvline(wstar, color="gray", ls="--", lw=0.8,
                       label=f"W* = {wstar}")
        ax.set_xscale("log")
        ax.set_xlabel("window W (samples)")
        ax.set_ylabel(title)
        ax.set_ylim(-1, 1)
        ax.axhline(0, color="black", lw=0.5)
        ax.legend(fontsize=8)
    fig.suptitle("Window-scale sweep, floor-free detection, ISM data "
                 "(hypothesis-generating)")
    fig.tight_layout()
    (ROOT / "results/figures/exploratory").mkdir(parents=True,
                                                 exist_ok=True)
    fig.savefig(ROOT / "results/figures/exploratory/window_sweep.png",
                dpi=130)

    # Markdown summary.
    L = ["# Window-scale sweep (HYPOTHESIS-GENERATING, ISM data)\n\n",
         "Grid and rules per analysis/PREREGISTRATION_FOLLOWUP.md. ",
         f"Frozen W* = {wstar}.\n\n",
         "| track | W | floor-free R2 | LOOCV R2 | Spearman | "
         "fixed-start R2 |\n|---|---|---|---|---|---|\n"]
    for track in ("pra", "ab"):
        for W in WINDOWS:
            c = detections[(track, W)]
            ff = c["floorfree"]
            fs = c["fixedstart"]
            L.append(f"| {track} | {W} | "
                     f"{ff.get('r2', float('nan')):.3f} | "
                     f"{ff.get('loocv_r2', float('nan')):.3f} | "
                     f"{ff.get('spearman', float('nan')):.3f} | "
                     f"{fs.get('r2', float('nan')):.3f} |\n")
    with open(ROOT / "results/exploratory/window_sweep.md", "w") as f:
        f.write("".join(L))
    print(f"W* frozen at {wstar}")


if __name__ == "__main__":
    main()
