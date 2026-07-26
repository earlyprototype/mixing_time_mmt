"""EXPLORATORY, NOT CONFIRMATORY. Floor-free detector variants, post hoc, n = 8.

The 2011 rule searches the smoothed Higuchi-FD profile from a FIXED start
(1-based sample 1000) for the first crossing of tail_mean - 2 * tail_sd.
On the band-limited pyroomacoustics track ("pra") the profile saturates
early and the detections pin at the search-start floor. This script
implements and evaluates floor-free variants. It selects nothing and
highlights nothing; every variant is reported identically.

Detectors (all return MATLAB-style 1-based sample numbers):
    D0     baseline for reference, the original fixed-start rule, read
           straight from results/rooms/room{N}.json (crossings key "2.0").
    D1     min-referenced start: same threshold rule, but the search
           starts at the argmin of the smoothed profile (finite samples).
    D2     rise fraction: with m = finite profile minimum and
           M = tail mean, first sample at or after the argmin where
           (p - m) / (M - m) >= q, for q in {0.5, 0.7, 0.9}.
    D3     peak-train re-sparsification, pra only: local maxima of
           abs(rir), a same-length sparse signal holding only the peak
           amplitudes, then the standard pipeline (hfd_profile window 50
           kmax 5, moving_average_nan span 200, fixed-start-1000 rule).

Evaluation, identical for every detector and track, no selection:
OLS of tmp50 (reconciled, samples) on detected samples, R2, Spearman rho,
LOOCV R2 (1 - PRESS/SStot), and the count of rooms whose detection lies
within 5 samples of its own search start (floor-pinned). Rooms where a
detector never crosses are recorded as null and excluded from that fit.

Deterministic: pure numpy/scipy on files already on disk, no randomness.

Writes (and nothing else):
    results/exploratory/detectors.json
    results/exploratory/detectors.md
    results/figures/exploratory/detectors_scatter.png
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analysis"))

from mixtime.profile import hfd_profile  # noqa: E402
from mixtime.detect import tail_stats  # noqa: E402
from run_room import moving_average_nan  # noqa: E402

EXPLORATORY_LABEL = (
    "EXPLORATORY, NOT CONFIRMATORY. Post-hoc detector variants at n = 8, "
    "devised after seeing the confirmatory results. No variant is selected "
    "or endorsed; all are reported identically."
)

ROOMS = list(range(1, 9))
TAIL_FRAC = 0.1
K = 2.0
FIXED_START = 1000  # 1-based, thesis convention
PIN_TOL = 5  # samples: detection within this of its search start = pinned
Q_LIST = [0.5, 0.7, 0.9]


def load_ground_truth():
    gt = json.loads((ROOT / "inputs/ground_truth.json").read_text())
    out = {}
    for r in gt["rooms"]:
        if r["room"] in ROOMS:
            out[r["room"]] = r["tmp50_samples_reconciled"]
    return out


def load_profile(room, track):
    suffix = "" if track == "pra" else "_ab"
    return np.load(ROOT / f"data/profiles/room{room}_smoothed{suffix}.npy")


def safe_tail_stats(profile):
    """Tail statistics with a nan guard.

    Uses detect.tail_stats (the exact MATLAB tail slice) when the tail is
    fully finite. If the tail contains nan (D3's sparse signal produces
    all-zero windows whose profile value is nan), the same slice is used
    but non-finite samples are excluded from the mean and sd (ddof=1).
    """
    mean, sd = tail_stats(profile, TAIL_FRAC)
    if np.isfinite(mean) and np.isfinite(sd):
        return mean, sd
    x = np.asarray(profile, dtype=float)
    n = x.size
    m = int(np.floor(n * TAIL_FRAC + 0.5))
    tail = x[n - m - 1:]
    tail = tail[np.isfinite(tail)]
    if tail.size < 2:
        raise ValueError("tail has fewer than 2 finite samples")
    return float(tail.mean()), float(tail.std(ddof=1))


def finite_argmin(profile):
    """0-based index of the minimum over finite samples, or None."""
    finite = np.isfinite(profile)
    if not finite.any():
        return None
    masked = np.where(finite, profile, np.inf)
    return int(np.argmin(masked))


def first_at_or_after(condition, i0):
    """1-based sample of the first True in condition at 0-based i0 on."""
    hits = np.nonzero(condition[i0:])[0]
    if hits.size == 0:
        return None
    return int(i0 + hits[0] + 1)


def d0_from_results(room, track):
    data = json.loads((ROOT / f"results/rooms/room{room}.json").read_text())
    key = "primary" if track == "pra" else "ab"
    det = data[key]["crossings"]["2.0"]
    return {
        "detected": None if det is None else int(det),
        "search_start": FIXED_START,
        "tail_mean": data[key]["tail_mean"],
        "tail_sd": data[key]["tail_sd"],
        "source": f"results/rooms/room{room}.json:{key}.crossings.2.0",
    }


def d1_min_start(profile):
    mean, sd = safe_tail_stats(profile)
    threshold = mean - K * sd
    i_min = finite_argmin(profile)
    if i_min is None:
        return {"detected": None, "search_start": None,
                "tail_mean": mean, "tail_sd": sd, "threshold": threshold}
    with np.errstate(invalid="ignore"):
        det = first_at_or_after(profile >= threshold, i_min)
    return {"detected": det, "search_start": i_min + 1,
            "argmin_0based": i_min, "profile_min": float(profile[i_min]),
            "tail_mean": mean, "tail_sd": sd, "threshold": threshold}


def d2_rise_fraction(profile, q):
    mean, _sd = safe_tail_stats(profile)
    i_min = finite_argmin(profile)
    if i_min is None:
        return {"detected": None, "search_start": None, "q": q,
                "tail_mean": mean}
    m = float(profile[i_min])
    denom = mean - m
    if denom <= 0:
        return {"detected": None, "search_start": i_min + 1, "q": q,
                "profile_min": m, "tail_mean": mean,
                "note": "tail mean not above profile minimum"}
    with np.errstate(invalid="ignore"):
        det = first_at_or_after((profile - m) / denom >= q, i_min)
    return {"detected": det, "search_start": i_min + 1, "q": q,
            "argmin_0based": i_min, "profile_min": m, "tail_mean": mean,
            "rise_target": m + q * denom}


def d3_peak_train(room):
    rir = np.load(ROOT / f"data/rirs/room{room}_rir.npy")
    absr = np.abs(np.asarray(rir, dtype=float))
    peaks, _ = find_peaks(absr, distance=1)
    sparse = np.zeros_like(absr)
    sparse[peaks] = absr[peaks]
    prof = hfd_profile(sparse, window=50, kmax=5)
    smoothed = moving_average_nan(prof, span=200)
    mean, sd = safe_tail_stats(smoothed)
    threshold = mean - K * sd
    with np.errstate(invalid="ignore"):
        det = first_at_or_after(smoothed >= threshold, FIXED_START - 1)
    return {"detected": det, "search_start": FIXED_START,
            "n_peaks": int(peaks.size), "tail_mean": mean, "tail_sd": sd,
            "threshold": threshold}


def evaluate(per_room, tmp50):
    """OLS, Spearman, LOOCV R2, floor-pinned count. No selection."""
    rooms_used, x, y, pinned = [], [], [], 0
    for room in ROOMS:
        rec = per_room[room]
        det = rec["detected"]
        if det is None:
            continue
        rooms_used.append(room)
        x.append(float(det))
        y.append(float(tmp50[room]))
        start = rec.get("search_start")
        if start is not None and abs(det - start) <= PIN_TOL:
            pinned += 1
    n = len(x)
    ev = {"n": n, "rooms_used": rooms_used, "floor_pinned_count": pinned,
          "floor_pinned_tolerance_samples": PIN_TOL,
          "slope": None, "intercept": None, "r2": None,
          "spearman_rho": None, "loocv_r2": None}
    if n < 3:
        ev["note"] = "fewer than 3 rooms with a detection, fit not attempted"
        return ev
    x = np.array(x)
    y = np.array(y)
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    sstot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / sstot if sstot > 0 else None
    rho = spearmanr(x, y).statistic
    press = 0.0
    for i in range(n):
        keep = np.ones(n, bool)
        keep[i] = False
        s_i, b_i = np.polyfit(x[keep], y[keep], 1)
        press += (y[i] - (s_i * x[i] + b_i)) ** 2
    loocv = 1.0 - press / sstot if sstot > 0 else None
    ev.update({"slope": float(slope), "intercept": float(intercept),
               "r2": r2, "spearman_rho": float(rho),
               "loocv_r2": None if loocv is None else float(loocv)})
    return ev


def main():
    tmp50 = load_ground_truth()
    profiles = {t: {r: load_profile(r, t) for r in ROOMS}
                for t in ("pra", "ab")}

    detectors = {
        "D0": {"description": "Baseline, original fixed-start rule "
                              "(first sample >= 1000 with profile >= "
                              "tail_mean - 2 tail_sd), read from "
                              "results/rooms/room{N}.json.",
               "tracks": {}},
        "D1": {"description": "Min-referenced start: same threshold rule, "
                              "search starts at the argmin of the smoothed "
                              "profile over finite samples.",
               "tracks": {}},
    }
    for q in Q_LIST:
        name = f"D2q{int(round(q * 100))}"
        detectors[name] = {
            "description": f"Rise fraction q = {q}: with m = finite profile "
                           "minimum and M = tail mean, first sample at or "
                           "after the argmin with (p - m)/(M - m) >= q.",
            "tracks": {}}
    detectors["D3"] = {
        "description": "Peak-train re-sparsification, pra only: local "
                       "maxima of abs(rir) (find_peaks, no height "
                       "threshold, distance 1) kept at their positions, "
                       "zeros elsewhere, then hfd_profile (window 50, "
                       "kmax 5), moving_average_nan span 200, and the "
                       "original fixed-start-1000 threshold rule.",
        "tracks": {}}

    for track in ("pra", "ab"):
        d0 = {r: d0_from_results(r, track) for r in ROOMS}
        d1 = {r: d1_min_start(profiles[track][r]) for r in ROOMS}
        detectors["D0"]["tracks"][track] = {
            "per_room": d0, "evaluation": evaluate(d0, tmp50)}
        detectors["D1"]["tracks"][track] = {
            "per_room": d1, "evaluation": evaluate(d1, tmp50)}
        for q in Q_LIST:
            name = f"D2q{int(round(q * 100))}"
            d2 = {r: d2_rise_fraction(profiles[track][r], q) for r in ROOMS}
            detectors[name]["tracks"][track] = {
                "per_room": d2, "evaluation": evaluate(d2, tmp50)}
    d3 = {r: d3_peak_train(r) for r in ROOMS}
    detectors["D3"]["tracks"]["pra"] = {
        "per_room": d3, "evaluation": evaluate(d3, tmp50)}

    order = ["D0", "D1", "D2q50", "D2q70", "D2q90", "D3"]

    out = {
        "label": EXPLORATORY_LABEL,
        "config": {
            "rooms": ROOMS, "tail_frac": TAIL_FRAC, "k": K,
            "fixed_start_1based": FIXED_START,
            "floor_pinned_tolerance_samples": PIN_TOL,
            "rise_fractions": Q_LIST,
            "target": "tmp50_samples_reconciled from inputs/ground_truth.json",
            "sample_convention": "1-based (MATLAB/thesis convention)",
            "fs_hz": 44100,
            "tail_stats_nan_guard": "detect.tail_stats on the exact MATLAB "
                                    "tail slice; if the tail contains nan "
                                    "(D3 sparse signal, all-zero windows), "
                                    "non-finite samples are excluded from "
                                    "the same slice (mean, sd ddof=1)",
        },
        "detector_order": order,
        "detectors": detectors,
    }
    (ROOT / "results/exploratory").mkdir(parents=True, exist_ok=True)
    (ROOT / "results/figures/exploratory").mkdir(parents=True, exist_ok=True)
    (ROOT / "results/exploratory/detectors.json").write_text(
        json.dumps(out, indent=1) + "\n")

    write_markdown(out, order)
    make_figure(detectors, tmp50, order)


def fmt(v, nd=3):
    return "null" if v is None else f"{v:.{nd}f}"


def write_markdown(out, order):
    lines = [
        "# Floor-free detector variants",
        "",
        f"**{EXPLORATORY_LABEL}**",
        "",
        "Target: tmp50_samples_reconciled (samples at 44100 Hz), rooms 1-8.",
        "OLS is tmp50 on detected sample. LOOCV R2 is 1 - PRESS/SStot over",
        "leave-one-room-out refits. Floor-pinned counts rooms whose",
        f"detection lies within {PIN_TOL} samples of that detector's own "
        "search start",
        "(fixed sample 1000 for D0 and D3, the profile argmin for D1 and "
        "D2).",
        "Rooms with no crossing are null and excluded from the fit "
        "(reported n).",
        "D3 is pra-only by design. All detections are 1-based sample "
        "numbers.",
        "",
        "| detector | track | n | floor-pinned | slope | intercept | R2 | "
        "Spearman rho | LOOCV R2 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for name in order:
        det = out["detectors"][name]
        for track in ("pra", "ab"):
            if track not in det["tracks"]:
                continue
            ev = det["tracks"][track]["evaluation"]
            lines.append(
                f"| {name} | {track} | {ev['n']} | "
                f"{ev['floor_pinned_count']} | {fmt(ev['slope'], 4)} | "
                f"{fmt(ev['intercept'], 1)} | {fmt(ev['r2'])} | "
                f"{fmt(ev['spearman_rho'])} | {fmt(ev['loocv_r2'])} |")
    lines += [
        "",
        "Full per-room detail (detections, search starts, thresholds, tail",
        "statistics) is in results/exploratory/detectors.json. Figure:",
        "results/figures/exploratory/detectors_scatter.png.",
        "",
        f"*{EXPLORATORY_LABEL}*",
        "",
    ]
    (ROOT / "results/exploratory/detectors.md").write_text("\n".join(lines))


def make_figure(detectors, tmp50, order):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"pra": "#2a78d6", "ab": "#eb6834"}
    surface = "#fcfcfb"
    ink = "#0b0b0b"
    ink2 = "#52514e"

    fig, axes = plt.subplots(2, 3, figsize=(13, 8.2))
    fig.patch.set_facecolor(surface)
    for ax, name in zip(axes.ravel(), order):
        det = detectors[name]
        ax.set_facecolor(surface)
        ax.grid(True, color="#e6e5e1", linewidth=0.8, zorder=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(ink2)
        ax.tick_params(colors=ink2, labelsize=8)
        for track in ("pra", "ab"):
            if track not in det["tracks"]:
                continue
            per_room = det["tracks"][track]["per_room"]
            for room in ROOMS:
                d = per_room[room]["detected"]
                if d is None:
                    continue
                ax.scatter(d, tmp50[room], s=48, color=colors[track],
                           edgecolors=surface, linewidths=1.2, zorder=3,
                           label=track if room == ROOMS[0] else None)
                off = (5, 4) if track == "pra" else (5, -10)
                ax.annotate(str(room), (d, tmp50[room]),
                            textcoords="offset points", xytext=off,
                            fontsize=7, color=ink2, zorder=4)
        tracks = [t for t in ("pra", "ab") if t in det["tracks"]]
        sub = " (pra only)" if tracks == ["pra"] else ""
        ax.set_title(f"{name}{sub}", fontsize=10, color=ink, loc="left")
        ax.set_xlabel("detected sample (1-based)", fontsize=8, color=ink2)
        ax.set_ylabel("tmp50 (samples)", fontsize=8, color=ink2)
        if len(tracks) > 1:
            leg = ax.legend(fontsize=8, frameon=False, loc="lower right")
            for txt in leg.get_texts():
                txt.set_color(ink2)
    fig.suptitle("Floor-free detector variants, detected sample vs "
                 "reconciled tmp50 (rooms 1-8)\n"
                 "EXPLORATORY, NOT CONFIRMATORY: post-hoc variants at "
                 "n = 8, no variant selected",
                 fontsize=11, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(ROOT / "results/figures/exploratory/detectors_scatter.png",
                dpi=150, facecolor=surface)
    plt.close(fig)


if __name__ == "__main__":
    main()
