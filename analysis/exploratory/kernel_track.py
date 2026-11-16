"""Third ISM rendering track: AB spike train + measurement-like kernel.

EXPLORATORY. Tests the hypothesis from the discussion with the author:
a real measurement chain applies a SHORT, CAUSAL, band-passed pulse to
each arrival (a thump), not an ideal symmetric sinc (a ring), so an ISM
that convolves the geometrically-exact Allen-Berkley spike train with
such a kernel should reproduce the early sparsity seen in measured
BRIRs, and with it the profile rise the 2011 feature needs.

Kernel: impulse response of a Butterworth bandpass (4th order lowpass
prototype, 8th order realized bandpass transfer function), 100 Hz to
16 kHz at fs = 44100, truncated at 3 ms. Causal and minimum-phase-like
by construction (IIR impulse response), front-loaded energy. It replaces
the thesis's butter(3, 0.01) highpass (its 100 Hz corner subsumes it).

For rooms 1 to 8: convolve the UNFILTERED sparse AB train with the
kernel, run the standard pipeline (window 50, kmax 5, smooth 200, tail
10 percent), detect Criteria I to IV with both the original fixed-start
rule and the floor-free rule, regress k = 2 against tmp50 (reconciled),
and record texture and rise statistics for comparison with the ab, pra
and measured tracks. Deterministic.
"""

import json
import pathlib
import sys

import numpy as np
from scipy.signal import butter, lfilter
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analysis"))

from mixtime.rooms import get_room
from mixtime.ab_ism import ab_rir
from mixtime.profile import hfd_profile
from mixtime.detect import tail_stats, detect_crossing
from run_room import moving_average_nan

FS = 44100
ROOMS = list(range(1, 9))
KERNEL_BAND = (100.0, 16000.0)
KERNEL_LEN = int(round(0.003 * FS))  # 3 ms truncation


def make_kernel():
    b, a = butter(4, [KERNEL_BAND[0], KERNEL_BAND[1]], btype="bandpass",
                  fs=FS)
    imp = np.zeros(KERNEL_LEN)
    imp[0] = 1.0
    k = lfilter(b, a, imp)
    e = np.cumsum(k ** 2)
    e /= e[-1]
    width_90 = int(np.searchsorted(e, 0.90)) + 1  # samples holding 90% energy
    return k, width_90


def ols_stats(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    sstot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ((y - A @ coef) ** 2).sum() / sstot
    n = x.size
    preds = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        c, _, _, _ = np.linalg.lstsq(A[m], y[m], rcond=None)
        preds[i] = A[i] @ c
    loocv = 1 - ((y - preds) ** 2).sum() / sstot
    return (float(coef[0]), float(coef[1]), float(r2), float(loocv),
            float(spearmanr(x, y)[0]))


def onset_index(rir, frac=0.05):
    """Amendment 1 onset: first sample at 5 percent of peak |amplitude|."""
    a = np.abs(np.asarray(rir, float))
    return int(np.argmax(a >= frac * a.max()))


def detect_floorfree(sm, rir, k=2.0):
    onset = onset_index(rir)
    finite = np.isfinite(sm)
    finite[:onset] = False
    if finite.sum() < 100:
        return None
    idx = np.arange(sm.size)[finite]
    amin = idx[np.argmin(sm[finite])]
    return detect_crossing(sm, k, tail_frac=0.1, start=int(amin) + 1)


def main():
    kernel, width_90 = make_kernel()
    gt = json.loads((ROOT / "inputs/ground_truth.json").read_text())
    by_room = {r["room"]: r["tmp50_samples_reconciled"]
               for r in gt["rooms"]}
    y = np.array([by_room[n] for n in ROOMS], float)

    out = {"note": "EXPLORATORY third rendering track: AB spikes + "
                   "causal measurement-like kernel.",
           "kernel": {"type": "butterworth bandpass impulse response, 4th order prototype, 8th order realized",
                      "band_hz": KERNEL_BAND,
                      "trunc_ms": 1000.0 * KERNEL_LEN / FS,
                      "width90_samples": width_90,
                      "width90_ms": 1000.0 * width_90 / FS},
           "rooms": {}}

    crossings = {"I": [], "II": [], "III": [], "IV": []}
    crossings_ff = []
    room1_rirs = {}
    for n in ROOMS:
        room = get_room(n)
        g = ab_rir(room, fs=FS)
        sparse = g["sparse_rir"]
        rir = np.convolve(sparse, kernel)  # full, keeps causal tails
        if n == 1:
            room1_rirs = {"sparse": sparse, "kernel_track": rir}
        prof = hfd_profile(rir, window=50, kmax=5)
        sm = moving_average_nan(prof, span=200)
        m, s = tail_stats(sm, 0.1)
        finite = sm[np.isfinite(sm)]
        seg = rir[:2000]
        entry = {
            "tail_mean": m, "tail_sd": s,
            "rise_over_tailsd": float((m - finite.min()) / s),
            "nearzero_frac_first2000":
                float(np.mean(np.abs(seg) < 1e-9 * np.abs(rir).max())),
            "crossings_fixedstart": {}, "crossing_k2_floorfree": None,
        }
        for crit, k in (("I", 0), ("II", 1), ("III", 2), ("IV", 3)):
            c = detect_crossing(sm, k, tail_frac=0.1, start=1000)
            entry["crossings_fixedstart"][crit] = c
            crossings[crit].append(c)
        ff = detect_floorfree(sm, rir)
        entry["crossing_k2_floorfree"] = ff
        crossings_ff.append(ff)
        out["rooms"][n] = entry

    # Regressions (rooms with detections; expect all 8).
    out["regressions"] = {}
    for label, xs in (("fixedstart_k2", crossings["III"]),
                      ("floorfree_k2", crossings_ff)):
        ok = [i for i, v in enumerate(xs) if v is not None]
        if len(ok) >= 4:
            sl, ic, r2, loocv, rho = ols_stats([xs[i] for i in ok], y[ok])
            out["regressions"][label] = {
                "n": len(ok), "slope": sl, "intercept": ic, "r2": r2,
                "loocv_r2": loocv, "spearman": rho}
    for crit in ("I", "II", "IV"):
        xs = crossings[crit]
        ok = [i for i, v in enumerate(xs) if v is not None]
        if len(ok) >= 4:
            sl, ic, r2, loocv, rho = ols_stats([xs[i] for i in ok], y[ok])
            out["regressions"][f"fixedstart_{crit}"] = {
                "n": len(ok), "r2": r2, "loocv_r2": loocv,
                "spearman": rho}

    (ROOT / "results/exploratory").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "results/exploratory/kernel_track.json", "w") as f:
        json.dump(out, f, indent=1)

    # Figure: Room 1 zoom, sparse vs kernel track vs pra, first 30 ms.
    pra1 = np.load(ROOT / "data/rirs/room1_rir.npy")
    nz = int(0.030 * FS)
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for ax, (label, sig) in zip(axes, [
            ("AB sparse spike train", room1_rirs["sparse"]),
            ("AB + measurement-like kernel (this track)",
             room1_rirs["kernel_track"]),
            ("pyroomacoustics (ideal sinc)", pra1)]):
        ax.plot(np.arange(nz) / FS * 1000, sig[:nz], lw=0.5)
        ax.set_title(label, fontsize=9)
        ax.set_ylabel("amp")
    axes[-1].set_xlabel("time (ms)")
    fig.suptitle("Room 1, first 30 ms: three renderings of the same "
                 "echo list")
    fig.tight_layout()
    fig.savefig(ROOT / "results/figures/exploratory/kernel_track_zoom.png",
                dpi=130)

    L = ["# AB + measurement-like kernel track (EXPLORATORY)\n\n",
         f"Kernel: causal Butterworth bandpass pulse, "
         f"{KERNEL_BAND[0]:.0f} Hz to {KERNEL_BAND[1]/1000:.0f} kHz, "
         f"90 percent of energy within {out['kernel']['width90_ms']:.2f} "
         "ms.\n\n",
         "| room | rise/tailSD | near-zero frac (first 2000) | "
         "Crit III (fixed start) | k=2 (floor-free) |\n|---|---|---|---|---|\n"]
    for n in ROOMS:
        e = out["rooms"][n]
        L.append(f"| {n} | {e['rise_over_tailsd']:.1f} | "
                 f"{e['nearzero_frac_first2000']:.2f} | "
                 f"{e['crossings_fixedstart']['III']} | "
                 f"{e['crossing_k2_floorfree']} |\n")
    L.append("\n| regression | n | R2 | LOOCV R2 | Spearman |\n"
             "|---|---|---|---|---|\n")
    for label, r in out["regressions"].items():
        L.append(f"| {label} | {r['n']} | {r['r2']:.3f} | "
                 f"{r['loocv_r2']:.3f} | {r['spearman']:.3f} |\n")
    with open(ROOT / "results/exploratory/kernel_track.md", "w") as f:
        f.write("".join(L))
    print(json.dumps(out["regressions"], indent=1))
    print("rises:", [round(out['rooms'][n]['rise_over_tailsd'], 1)
                     for n in ROOMS])


if __name__ == "__main__":
    main()
