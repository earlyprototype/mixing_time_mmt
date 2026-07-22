"""Process one Lindau room end to end: ISM RIR, HFD profile, detections.

Usage: python3 analysis/run_room.py <room_number>

Writes:
    data/rirs/room<N>_rir.npy            primary RIR (length round(RT*fs))
    data/profiles/room<N>_profile.npy    unsmoothed HFD profile (primary)
    data/profiles/room<N>_smoothed.npy   smoothed HFD profile (primary)
    results/rooms/room<N>.json           all numbers (detections, sweeps,
                                         sensitivity variants, exploratory)

Everything is deterministic (pure ISM, no randomness). The configuration
grid follows analysis/PREREGISTRATION.md exactly. This script computes and
records; it draws no conclusions.
"""

import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from mixtime.rooms import get_room
from mixtime.ism import generate_rir, default_length_samples
from mixtime.ab_ism import ab_rir
from mixtime.rt_check import schroeder_rt60
from mixtime.profile import hfd_profile, smooth_profile
from mixtime.detect import tail_stats, detect_crossing

FS = 44100
K_SWEEP = [round(0.25 * i, 2) for i in range(17)]  # 0.0 .. 4.0
PRIMARY = dict(window=50, kmax=5, span=200, tail_frac=0.10, start=1000,
               length_factor=1.00)
# One factor at a time relative to PRIMARY, per the pre-registration.
SENSITIVITY = (
    [("window", w) for w in (25, 100)]
    + [("kmax", km) for km in (3, 8)]
    + [("tail_frac", t) for t in (0.05, 0.20)]
    + [("length_factor", f) for f in (0.75, 1.25)]
    + [("start", s) for s in (500, 2000)]
)


def moving_average_nan(x, span=200):
    """Centered moving average with MATLAB smooth() edge shrinking,
    tolerant of nans (nans are ignored inside each window)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    half_span = (span - 1) // 2 if span % 2 == 0 else span // 2
    good = np.isfinite(x)
    xz = np.where(good, x, 0.0)
    cs = np.concatenate([[0.0], np.cumsum(xz)])
    cn = np.concatenate([[0.0], np.cumsum(good.astype(float))])
    idx = np.arange(n)
    half = np.minimum(np.minimum(idx, n - 1 - idx), half_span)
    lo = idx - half
    hi = idx + half + 1
    tot = cs[hi] - cs[lo]
    cnt = cn[hi] - cn[lo]
    out = np.full(n, np.nan)
    nz = cnt > 0
    out[nz] = tot[nz] / cnt[nz]
    return out


def sampen_profile(rir, window=50, m=2, r_factor=0.2):
    """Sample entropy of each 51-sample moving window (zero padded at the
    end, same segmentation as the HFD profile). Exploratory measure.

    SampEn(m, r) = -log(A / B), r = r_factor * std of the window. Windows
    with B == 0 or A == 0 get nan (no matching templates, entropy
    undefined). Vectorized over windows in chunks.
    """
    x = np.asarray(rir, dtype=float)
    n = x.size
    seg_len = window + 1
    padded = np.concatenate([x, np.zeros(window)])
    n_out = n - 2  # same length convention as hfd_profile
    from numpy.lib.stride_tricks import sliding_window_view
    segs = sliding_window_view(padded, seg_len)[:n_out]
    out = np.full(n_out, np.nan)
    chunk = 256
    for c0 in range(0, n_out, chunk):
        S = segs[c0:c0 + chunk]  # (B, seg_len)
        B_ = S.shape[0]
        sd = S.std(axis=1, ddof=0)
        r = r_factor * sd
        # Template vectors of length m and m+1.
        n_tm = seg_len - m       # number of length-m templates
        Tm = sliding_window_view(S, m, axis=1)        # (B, n_tm+1, m)
        Tm1 = sliding_window_view(S, m + 1, axis=1)   # (B, n_tm, m+1)
        Tm = Tm[:, :n_tm, :]  # use first n_tm so counts B and A align
        # Chebyshev distances between all template pairs, per window.
        dm = np.abs(Tm[:, :, None, :] - Tm[:, None, :, :]).max(axis=3)
        dm1 = np.abs(Tm1[:, :, None, :] - Tm1[:, None, :, :]).max(axis=3)
        iu = np.triu_indices(n_tm, k=1)
        rb = r[:, None]
        Bcount = (dm[:, iu[0], iu[1]] <= rb).sum(axis=1)
        n_tm1 = Tm1.shape[1]
        iu1 = np.triu_indices(n_tm1, k=1)
        Acount = (dm1[:, iu1[0], iu1[1]] <= rb).sum(axis=1)
        valid = (Bcount > 0) & (Acount > 0) & (sd > 0)
        vals = np.full(B_, np.nan)
        vals[valid] = -np.log(Acount[valid] / Bcount[valid])
        out[c0:c0 + chunk] = vals
    return out


def katz_profile(rir, window=50):
    """Katz fractal dimension of each 51-sample moving window (zero padded
    at the end). Exploratory measure.

    Katz FD = log10(nsteps) / (log10(nsteps) + log10(d / Ltot)) with Ltot
    the summed point-to-point path length and d the maximum distance from
    the first point. Windows with zero path length get nan.
    """
    x = np.asarray(rir, dtype=float)
    n = x.size
    seg_len = window + 1
    padded = np.concatenate([x, np.zeros(window)])
    n_out = n - 2
    from numpy.lib.stride_tricks import sliding_window_view
    segs = sliding_window_view(padded, seg_len)[:n_out]
    steps = np.abs(np.diff(segs, axis=1))
    Ltot = steps.sum(axis=1)
    d = np.abs(segs - segs[:, [0]]).max(axis=1)
    nsteps = seg_len - 1
    out = np.full(n_out, np.nan)
    ok = (Ltot > 0) & (d > 0)
    logn = np.log10(nsteps)
    out[ok] = logn / (logn + np.log10(d[ok] / Ltot[ok]))
    return out


def crossings_for_profile(smoothed, tail_frac, start, ks):
    stats = tail_stats(smoothed, tail_frac)
    res = {}
    for k in ks:
        res[str(k)] = detect_crossing(smoothed, k, tail_frac=tail_frac,
                                      start=start)
    return {"tail_mean": stats[0], "tail_sd": stats[1], "crossings": res}


def process(room_number):
    t_start = time.time()
    room = get_room(room_number)
    base_len = default_length_samples(room, FS)
    out = {
        "room": room_number,
        "name": room.name,
        "volume_m3": room.volume_m3,
        "alpha_ave": room.alpha_ave,
        "rt_s": room.rt_s,
        "included_in_2011": room.included,
        "dims_m": [room.L, room.W, room.H],
        "fs": FS,
        "base_length_samples": base_len,
        "primary": None,
        "ab": None,
        "k_sweep_ks": K_SWEEP,
        "sensitivity": [],
        "sensitivity_ab": [],
        "exploratory": {},
    }

    rir_cache = {}

    def rir_for(track, length_factor):
        """track 'pra' is the pre-registered pyroomacoustics path, track
        'ab' the thesis-faithful Allen-Berkley port (integer delays,
        beta = sqrt(1 - alpha), butter(3, 0.01) highpass)."""
        key = (track, length_factor)
        if key not in rir_cache:
            L = int(round(base_len * length_factor))
            if track == "pra":
                rir_cache[key] = generate_rir(room, fs=FS, length_samples=L)
            else:
                rir_cache[key] = ab_rir(room, fs=FS, length_samples=L)
        return rir_cache[key]

    (ROOT / "data/rirs").mkdir(parents=True, exist_ok=True)
    (ROOT / "data/profiles").mkdir(parents=True, exist_ok=True)
    (ROOT / "results/rooms").mkdir(parents=True, exist_ok=True)

    smoothed_by_track = {}
    for track, result_key, sens_key in (("pra", "primary", "sensitivity"),
                                        ("ab", "ab", "sensitivity_ab")):
        g = rir_for(track, 1.00)
        rir = g["rir"]
        prof = hfd_profile(rir, window=PRIMARY["window"],
                           kmax=PRIMARY["kmax"])
        smoothed = moving_average_nan(prof, span=PRIMARY["span"])
        smoothed_by_track[track] = smoothed
        res = crossings_for_profile(smoothed, PRIMARY["tail_frac"],
                                    PRIMARY["start"], K_SWEEP)
        res["length_samples"] = g["length_samples"]
        res["profile_len"] = int(smoothed.size)
        res["schroeder_rt60"] = schroeder_rt60(
            g["raw_rir"] if track == "pra" else rir, FS)
        if track == "pra":
            res["max_order"] = g["max_order"]
            res["order_capped"] = bool(g["order_capped"])
            res["raw_length"] = g["raw_length"]
        else:
            res["n_images"] = g["n_images"]
            res["beta"] = g["beta"]
        out[result_key] = res

        suffix = "" if track == "pra" else "_ab"
        np.save(ROOT / f"data/rirs/room{room_number}_rir{suffix}.npy", rir)
        np.save(ROOT / f"data/profiles/room{room_number}_profile{suffix}.npy",
                prof)
        np.save(ROOT / f"data/profiles/room{room_number}_smoothed{suffix}.npy",
                smoothed)

        # Sensitivity variants, one factor at a time, k=2 only.
        for name, value in SENSITIVITY:
            cfg = dict(PRIMARY)
            cfg[name] = value
            gv = rir_for(track, cfg["length_factor"])
            if name in ("window", "kmax", "length_factor"):
                p = hfd_profile(gv["rir"], window=cfg["window"],
                                kmax=cfg["kmax"])
            s = moving_average_nan(p, span=cfg["span"]) \
                if name in ("window", "kmax", "length_factor") else smoothed
            r = crossings_for_profile(s, cfg["tail_frac"], cfg["start"],
                                      [2.0])
            out[sens_key].append({
                "varied": name, "value": value,
                "crossing_k2": r["crossings"]["2.0"],
                "tail_mean": r["tail_mean"], "tail_sd": r["tail_sd"],
            })
    rir = rir_for("pra", 1.00)["rir"]
    smoothed = smoothed_by_track["pra"]

    # Exploratory measures, primary config, k=2 only.
    se = sampen_profile(rir, window=PRIMARY["window"])
    se_s = moving_average_nan(se, span=PRIMARY["span"])
    ka = katz_profile(rir, window=PRIMARY["window"])
    ka_s = moving_average_nan(ka, span=PRIMARY["span"])
    for label, arr in (("sampen", se_s), ("katz", ka_s)):
        finite = np.isfinite(arr)
        if finite.sum() > 2000:
            r = crossings_for_profile(np.where(finite, arr, np.nanmedian(arr)),
                                      PRIMARY["tail_frac"], PRIMARY["start"],
                                      [2.0])
            out["exploratory"][label] = {
                "crossing_k2": r["crossings"]["2.0"],
                "tail_mean": r["tail_mean"], "tail_sd": r["tail_sd"],
                "nan_fraction": float(1 - finite.mean()),
            }
        else:
            out["exploratory"][label] = {"crossing_k2": None,
                                         "nan_fraction": float(1 - finite.mean()),
                                         "note": "profile mostly undefined"}
    np.save(ROOT / f"data/profiles/room{room_number}_sampen.npy", se_s)
    np.save(ROOT / f"data/profiles/room{room_number}_katz.npy", ka_s)

    out["elapsed_s"] = round(time.time() - t_start, 1)
    path = ROOT / f"results/rooms/room{room_number}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    return out, path


if __name__ == "__main__":
    number = int(sys.argv[1])
    out, path = process(number)
    print(f"room {number} done in {out['elapsed_s']} s -> {path}")
    for key in ("primary", "ab"):
        p = out[key]
        print(f"  [{key}] length {p['length_samples']}, "
              f"RT60 {p['schroeder_rt60']:.3f} vs tab {out['rt_s']}, "
              f"tail mean {p['tail_mean']:.4f} sd {p['tail_sd']:.4f}")
        for c in ("0.0", "1.0", "2.0", "3.0"):
            print(f"    k={c}: sample {p['crossings'][c]}")
