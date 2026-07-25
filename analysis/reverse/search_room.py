"""Reverse-engineer the 2011 per-room inputs from the printed outputs.

Usage: python3 analysis/reverse/search_room.py <room_number> [budget_evals]

For one room, search the unrecorded 2011 inputs (footprint aspect, height,
RIR length, beta convention, axis interpretation) for the configuration
whose full pipeline output reproduces the thesis's printed fingerprints:
Table 6.1a tail mean and SD (4 decimals) and Table 6.1b detected samples
for Criteria I to IV (exact integers). Six observables against roughly
three continuous unknowns, so an exact hit identifies the inputs.

Honest interpretation rules, also encoded in the output JSON:
- exact_match true requires all four crossings equal and both tail stats
  within 5e-5 (print rounding) of the printed values.
- Anything else is a best fit with residuals; it constrains but does not
  identify the 2011 inputs, and a systematic near-miss pattern points to
  a forward-model difference (something the thesis GUI did that the
  printed appendix code does not show).

Deterministic: fixed grids, deterministic Nelder-Mead starts, no RNG.
Writes results/reverse/room<N>.json and a progress log to
results/reverse/room<N>.log.
"""

import json
import math
import pathlib
import sys
import time
import types

import numpy as np
from scipy.optimize import minimize

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analysis"))

from mixtime.rooms import get_room
from mixtime.ab_ism import ab_rir
from mixtime.profile import hfd_profile
from mixtime.detect import tail_stats, detect_crossing
from run_room import moving_average_nan

FS = 44100

# Printed fingerprints. Table 6.1a tail mean, SD; Table 6.1b crossings.
TARGETS = {
    1: {"mean": 2.0364, "sd": 0.0252, "cross": [6277, 5772, 3136, 2053]},
    2: {"mean": 2.0450, "sd": 0.0265, "cross": [4573, 3197, 3170, 3049]},
    3: {"mean": 2.0250, "sd": 0.0321, "cross": [6056, 4335, 2494, 1567]},
    4: {"mean": 1.9833, "sd": 0.0394, "cross": [13724, 7027, 7013, 7004]},
    5: {"mean": 1.9782, "sd": 0.0404, "cross": [7996, 7971, 7906, 7661]},
    6: {"mean": 1.9866, "sd": 0.0432, "cross": [12366, 9486, 8060, 2052]},
    7: {"mean": 1.9682, "sd": 0.0394, "cross": [14850, 12308, 11708, 2192]},
    8: {"mean": 1.9644, "sd": 0.0506, "cross": [16194, 12571, 8465, 5209]},
}


def make_room(base, L, W, H, beta_formula):
    """Lightweight room stand-in accepted by ab_rir."""
    r = types.SimpleNamespace()
    r.L, r.W, r.H = float(L), float(W), float(H)
    r.volume_m3 = base.volume_m3
    r.rt_s = base.rt_s
    # ab_rir reads alpha_ave and applies sqrt(1 - alpha). To emulate the
    # alternative convention beta = 1 - alpha, feed the alpha whose sqrt
    # gives that beta.
    if beta_formula == "sqrt":
        r.alpha_ave = base.alpha_ave
    else:  # beta = 1 - alpha
        beta = 1.0 - base.alpha_ave
        r.alpha_ave = 1.0 - beta * beta
    return r


def forward(base, L, W, H, npts, beta_formula, swap_axis):
    """Full pipeline for one candidate configuration."""
    if swap_axis:
        L, W = W, L
    room = make_room(base, L, W, H, beta_formula)
    g = ab_rir(room, fs=FS, length_samples=int(npts))
    prof = hfd_profile(g["rir"], window=50, kmax=5)
    sm = moving_average_nan(prof, span=200)
    mean, sd = tail_stats(sm, 0.1)
    cross = [detect_crossing(sm, k, tail_frac=0.1, start=1000)
             for k in (0, 1, 2, 3)]
    return mean, sd, cross


def score(target, mean, sd, cross):
    """Lower is better. Relative log distance on crossings plus scaled
    absolute distance on tail stats. A missing crossing is heavily
    penalized."""
    s = 0.0
    for got, want in zip(cross, target["cross"]):
        if got is None:
            s += 10.0
        else:
            s += abs(math.log(got / want))
    s += abs(mean - target["mean"]) / 0.001
    s += abs(sd - target["sd"]) / 0.001
    return s


def is_exact(target, mean, sd, cross):
    return (all(g == w for g, w in zip(cross, target["cross"]))
            and abs(mean - target["mean"]) <= 5e-5
            and abs(sd - target["sd"]) <= 5e-5)


def candidate_npts(base):
    """RIR length candidates: fractions of RT*fs plus the two lengths
    the thesis document itself shows (0.2 s figure, 0.256 s GUI)."""
    rt = int(round(base.rt_s * FS))
    vals = {8820, 11290, rt}
    for f in np.arange(0.50, 1.51, 0.05):
        vals.add(int(round(rt * f)))
    return sorted(v for v in vals if v >= 5000)


def main(room_number, budget=6000):
    t0 = time.time()
    base = get_room(room_number)
    target = TARGETS[room_number]
    outdir = ROOT / "results/reverse"
    outdir.mkdir(parents=True, exist_ok=True)
    logf = open(outdir / f"room{room_number}.log", "w")

    def log(msg):
        logf.write(msg + "\n")
        logf.flush()

    ratio_dig = base.L / base.W
    ratios = ratio_dig * np.exp(np.linspace(math.log(0.7), math.log(1.3), 9))
    heights = [h for h in (2.4, 2.7, 3.0, 3.3, 3.6, 4.0, 4.5, 5.0, 6.0,
                           7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 17.0, 20.0)
               if 2.2 <= h <= 22.0 and base.volume_m3 / h >= 16.0]
    npts_list = candidate_npts(base)

    log(f"room {room_number}: grid {len(ratios)} ratios x {len(heights)} "
        f"heights x {len(npts_list)} lengths x 2 beta x 2 swap")

    evals = 0
    results = []
    best_exact = None
    # Coarse pass. Loop order puts npts outermost so early exact matches
    # on plausible lengths are found soon.
    for npts in npts_list:
        for beta_formula in ("sqrt", "one_minus"):
            for swap in (False, True):
                for r in ratios:
                    for H in heights:
                        A = base.volume_m3 / H
                        W = math.sqrt(A / r)
                        L = r * W
                        if min(L, W) < 4.0:
                            continue
                        mean, sd, cross = forward(base, L, W, H, npts,
                                                  beta_formula, swap)
                        evals += 1
                        sc = score(target, mean, sd, cross)
                        results.append((sc, L, W, H, npts, beta_formula,
                                        swap, mean, sd, cross))
                        if is_exact(target, mean, sd, cross):
                            best_exact = results[-1]
                            log(f"EXACT at eval {evals}: {results[-1]}")
                        if evals >= budget:
                            break
                    if evals >= budget or best_exact:
                        break
                if evals >= budget or best_exact:
                    break
            if evals >= budget or best_exact:
                break
        if evals >= budget or best_exact:
            break
    results.sort(key=lambda t: t[0])
    log(f"coarse done: {evals} evals in {time.time()-t0:.0f} s, "
        f"best score {results[0][0]:.4f}")

    # Local refinement of the top candidates over (log ratio, H).
    refined = []
    for cand in results[:12]:
        sc0, L0, W0, H0, npts, bf, swap, *_ = cand

        def obj(x):
            lr, H = x
            if not (2.0 <= H <= 24.0):
                return 1e6
            r = math.exp(lr)
            A = base.volume_m3 / H
            W = math.sqrt(A / r)
            L = r * W
            if min(L, W) < 3.5:
                return 1e6
            mean, sd, cross = forward(base, L, W, H, npts, bf, swap)
            return score(target, mean, sd, cross)

        x0 = [math.log(L0 / W0), H0]
        res = minimize(obj, x0, method="Nelder-Mead",
                       options={"maxfev": 120, "xatol": 1e-3,
                                "fatol": 1e-4})
        lr, H = res.x
        r = math.exp(lr)
        A = base.volume_m3 / H
        W = math.sqrt(A / r)
        L = r * W
        mean, sd, cross = forward(base, L, W, H, npts, bf, swap)
        refined.append((score(target, mean, sd, cross), L, W, H, npts, bf,
                        swap, mean, sd, cross))
        if is_exact(target, mean, sd, cross):
            best_exact = refined[-1]
            log(f"EXACT after refine: {refined[-1]}")
            break
    refined.sort(key=lambda t: t[0])
    best = best_exact or (refined[0] if refined[0][0] < results[0][0]
                          else results[0])

    def pack(t):
        sc, L, W, H, npts, bf, swap, mean, sd, cross = t
        return {"score": sc, "L": L, "W": W, "H": H, "npts": int(npts),
                "beta_formula": bf, "swap_axis": bool(swap),
                "tail_mean": mean, "tail_sd": sd,
                "crossings": cross,
                "target_mean": target["mean"], "target_sd": target["sd"],
                "target_crossings": target["cross"],
                "exact_match": is_exact(target, mean, sd, cross)}

    out = {
        "room": room_number,
        "evals": evals,
        "elapsed_s": round(time.time() - t0, 1),
        "digitized_geometry": {"L": base.L, "W": base.W, "H": base.H},
        "best": pack(best),
        "top10_coarse": [pack(t) for t in results[:10]],
        "top5_refined": [pack(t) for t in refined[:5]],
        "note": ("exact_match false means the printed fingerprints were "
                 "NOT reproduced bit-for-bit; treat the best fit as a "
                 "constraint, not an identification."),
    }
    with open(outdir / f"room{room_number}.json", "w") as f:
        json.dump(out, f, indent=1)
    log(f"done in {time.time()-t0:.0f} s; best score {best[0]:.4f}, "
        f"exact={out['best']['exact_match']}")
    print(json.dumps(out["best"], indent=1))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python3 analysis/reverse/search_room.py <room> "
                 "[budget_evals]")
    main(int(sys.argv[1]),
         int(sys.argv[2]) if len(sys.argv) > 2 else 6000)
