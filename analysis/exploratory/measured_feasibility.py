"""H1 feasibility gate on measured BRIRs. Per PREREGISTRATION_FOLLOWUP.md.

Data: University of Surrey IoSR RealRoomBRIRs (github.com/IoSR-Surrey/
RealRoomBRIRs, Rooms A to D, 48 kHz SOFA), measured binaural RIRs of
four real rooms. No perceptual mixing-time labels exist for these rooms,
so ONLY hypothesis H1 (does the smoothed HFD profile exhibit a
measurable rise on measured data) is tested here. H2 is untestable on
this corpus by design and no correlational statistic is computed.

H1 criterion (fixed in the pre-registration before this data was seen):
rise = (tail mean minus profile minimum) / tail SD > 3, in at least
three quarters of tested rooms at some window W in the grid.

Frontal head orientation, left ear, per the pre-registration's stated
preference. Files are fetched by fetch() if absent (they are not
committed, about 150 MB total; provenance URLs recorded).
"""

import hashlib
import json
import pathlib
import subprocess
import sys

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "analysis"))

from mixtime.profile import hfd_profile
from mixtime.detect import tail_stats
from run_room import moving_average_nan

DATA = ROOT / "data/measured"
BASE = ("https://raw.githubusercontent.com/IoSR-Surrey/RealRoomBRIRs/"
        "master/")
FILES = {r: f"UniS_Room_{r}_BRIR_48k.sofa" for r in "ABCD"}
WINDOWS = [50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000]
RISE_THRESHOLD = 3.0


def fetch():
    """Atomic download: curl -f fails on HTTP errors, temp file renamed
    only on success, so an interrupted or 404 response can never poison
    the cache."""
    DATA.mkdir(parents=True, exist_ok=True)
    for name in FILES.values():
        p = DATA / name
        if p.exists():
            continue
        tmp = p.with_suffix(p.suffix + ".part")
        subprocess.run(["curl", "-fsSL", "--max-time", "300", "-o",
                        str(tmp), BASE + name], check=True)
        tmp.replace(p)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frontal_left(path, tolerance_deg=5.0):
    """Left-ear IR at the frontal source azimuth (0 degrees). Asserts
    the nearest available azimuth is within tolerance_deg so an off-axis
    IR can never silently stand in for frontal."""
    with h5py.File(path, "r") as f:
        az = f["SourcePosition"][:][:, 0] % 360.0
        dist = np.minimum(az, 360.0 - az)
        idx = int(np.argmin(dist))
        if dist[idx] > tolerance_deg:
            raise ValueError(
                f"no frontal measurement within {tolerance_deg} deg in "
                f"{path} (nearest {dist[idx]:.1f} deg)")
        ir = f["Data.IR"][idx, 0, :].astype(np.float64)
        fs = float(f["Data.SamplingRate"][:][0])
    return ir, fs, float(az[idx])


def main():
    fetch()
    out = {"note": "H1 feasibility only; no perceptual labels exist for "
                   "this corpus, H2 not testable here.",
           "corpus": "IoSR RealRoomBRIRs Rooms A to D, 48 kHz, frontal "
                     "azimuth, left ear",
           "provenance": BASE,
           "sha256": {name: sha256(DATA / name)
                      for name in FILES.values()},
           "windows": WINDOWS,
           "rise_threshold": RISE_THRESHOLD, "rooms": {}}
    profiles_for_plot = {}
    for room, name in FILES.items():
        ir, fs, az = frontal_left(DATA / name)
        entry = {"file": name, "fs": fs, "azimuth_used": az,
                 "n_samples": int(ir.size), "windows": {}}
        for W in WINDOWS:
            prof = hfd_profile(ir, window=W, kmax=5)
            sm = moving_average_nan(prof, span=4 * W)
            m, s = tail_stats(sm, 0.1)
            finite = sm[np.isfinite(sm)]
            rise = float((m - finite.min()) / s) if s > 0 else float("nan")
            entry["windows"][W] = {"tail_mean": m, "tail_sd": s,
                                   "profile_min": float(finite.min()),
                                   "rise_over_tailsd": rise,
                                   "passes": bool(rise > RISE_THRESHOLD)}
            if W in (50, 500):
                profiles_for_plot[(room, W)] = sm
        out["rooms"][room] = entry

    # H1 verdict per the pre-registered criterion.
    per_window_pass = {W: sum(1 for r in out["rooms"].values()
                              if r["windows"][W]["passes"])
                       for W in WINDOWS}
    n_rooms = len(FILES)
    h1 = any(c >= int(np.ceil(0.75 * n_rooms))
             for c in per_window_pass.values())
    out["per_window_pass_count"] = per_window_pass
    out["h1_pass"] = bool(h1)

    (ROOT / "results/exploratory").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "results/exploratory/measured_feasibility.json",
              "w") as f:
        json.dump(out, f, indent=1)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, room in zip(axes.ravel(), "ABCD"):
        for W, color in ((50, "tab:blue"), (500, "tab:orange")):
            ax.plot(profiles_for_plot[(room, W)], lw=0.6, color=color,
                    label=f"W={W}")
        ax.set_title(f"Room {room} (measured BRIR)", fontsize=10)
        ax.set_xlabel("sample")
        ax.set_ylabel("HFD")
        ax.legend(fontsize=8)
    fig.suptitle("Smoothed HFD profiles of measured BRIRs "
                 "(H1 feasibility, exploratory)")
    fig.tight_layout()
    (ROOT / "results/figures/exploratory").mkdir(parents=True,
                                                 exist_ok=True)
    fig.savefig(ROOT / "results/figures/exploratory/"
                       "measured_feasibility.png", dpi=130)

    L = ["# H1 feasibility on measured BRIRs (IoSR Rooms A to D)\n\n",
         f"H1 criterion: rise > {RISE_THRESHOLD} x tail SD in >= 3 of 4 "
         "rooms at some window. ",
         f"**H1 {'PASSES' if h1 else 'FAILS'}.**\n\n",
         "| W | rooms passing | rises (A, B, C, D) |\n|---|---|---|\n"]
    for W in WINDOWS:
        rises = ", ".join(f"{out['rooms'][r]['windows'][W]['rise_over_tailsd']:.1f}"
                          for r in "ABCD")
        L.append(f"| {W} | {per_window_pass[W]}/4 | {rises} |\n")
    with open(ROOT / "results/exploratory/measured_feasibility.md",
              "w") as f:
        f.write("".join(L))
    print(f"H1 pass: {h1}")
    for W in WINDOWS:
        print(W, per_window_pass[W])


if __name__ == "__main__":
    main()
