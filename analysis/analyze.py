"""Cross-room validation analysis, per analysis/PREREGISTRATION.md.

Inputs: results/rooms/room{1..8}.json (both RIR tracks) and
inputs/ground_truth.json (digitized and reconciled tmp50 percent).

Outputs: results/analysis.json, results/summary.md, and figures under
results/figures/.

Everything deterministic: bootstrap seed 20260722, 10000 resamples.
"""

import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

FS = 44100
SEED = 20260722
N_BOOT = 10000
INCLUDED = [1, 2, 3, 4, 5, 6, 7, 8]
THESIS = {
    "samples": {  # Table 6.1b, rooms 1..8
        "I": [6277, 4573, 6056, 13724, 7996, 12366, 14850, 16194],
        "II": [5772, 3197, 4335, 7027, 7971, 9486, 12308, 12571],
        "III": [3136, 3170, 2494, 7013, 7906, 8060, 11708, 8465],
        "IV": [2053, 3049, 1567, 7004, 7661, 2052, 2192, 5209],
    },
    "regressions": {  # Table 6.2: slope, intercept, R2 percent
        "I": (0.2090, 331, 72.29),
        "II": (0.2760, 293, 78.82),
        "III": (0.3197, 325, 93.49),
        "IV": (0.0517, 2202, 0.14),
    },
}
CRIT_K = {"I": "0.0", "II": "1.0", "III": "2.0", "IV": "3.0"}


def ols(x, y):
    """OLS y on x. Returns slope, intercept, r2."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(coef[0]), float(coef[1]), r2


def loocv(x, y):
    """Leave-one-out CV of the linear fit. Returns predictions, r2, rmse."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    preds = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        s, b, _ = ols(x[m], y[m])
        preds[i] = s * x[i] + b
    press = float(((y - preds) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - press / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(press / n))
    return preds, r2, rmse


def bootstrap_ci(x, y, n_boot=N_BOOT, seed=SEED):
    """Percentile bootstrap CIs for R2 and slope, resampling rooms."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    r2s, slopes = [], []
    draws = 0
    while len(r2s) < n_boot and draws < n_boot * 20:
        idx = rng.integers(0, n, n)
        draws += 1
        if np.unique(idx).size < 3:
            continue
        s, _, r2 = ols(x[idx], y[idx])
        if np.isfinite(r2):
            r2s.append(r2)
            slopes.append(s)
    r2s = np.array(r2s)
    slopes = np.array(slopes)
    pct = lambda a: [float(np.percentile(a, 2.5)),
                     float(np.percentile(a, 97.5))]
    return {"r2_ci": pct(r2s), "slope_ci": pct(slopes),
            "n_effective": int(r2s.size)}


def main():
    rooms = {}
    for n in range(1, 10):
        p = ROOT / f"results/rooms/room{n}.json"
        if p.exists():
            rooms[n] = json.loads(p.read_text())
    gt = json.loads((ROOT / "inputs/ground_truth.json").read_text())
    gt_rooms = {int(r["room"]): r for r in gt["rooms"]}

    missing = [n for n in INCLUDED if n not in rooms]
    if missing:
        raise SystemExit(f"missing room results: {missing}")

    y_rec = np.array([gt_rooms[n]["tmp50_samples_reconciled"]
                      for n in INCLUDED], float)
    y_dig = np.array([gt_rooms[n]["tmp50_samples_digitized"]
                      for n in INCLUDED], float)

    out = {"included_rooms": INCLUDED, "seed": SEED, "n_boot": N_BOOT,
           "ground_truth_samples_reconciled": y_rec.tolist(),
           "ground_truth_samples_digitized": y_dig.tolist(),
           "tracks": {}}

    # Thesis baseline, recomputed rather than quoted: regress the thesis's
    # OWN detected samples (its Table 6.1b) against the recovered ground
    # truth. This substantiates the "2011 arithmetic reproduces" claim
    # inside this script instead of relying on the digitization stage.
    out["thesis_baseline"] = {}
    for crit in ("I", "II", "III", "IV"):
        x = np.array(THESIS["samples"][crit], float)
        s, b, r2 = ols(x, y_rec)
        sd_, bd_, r2d_ = ols(x, y_dig)
        ts, tb, tr2 = THESIS["regressions"][crit]
        out["thesis_baseline"][crit] = {
            "slope": s, "intercept": b, "r2": r2,
            "r2_digitized_gt": r2d_,
            "published": {"slope": ts, "intercept": tb, "r2_percent": tr2},
        }

    k_sweep = [str(k) for k in rooms[1]["k_sweep_ks"]]

    for track_key, res_key, sens_key in (("pra", "primary", "sensitivity"),
                                         ("ab", "ab", "sensitivity_ab")):
        T = {"crossings": {}, "k_sweep": [], "criteria": {},
             "sensitivity": [], "exploratory": {}}

        def xs_for(kstr):
            vals = [rooms[n][res_key]["crossings"].get(kstr)
                    for n in INCLUDED]
            return vals

        # Full k sweep on reconciled ground truth.
        for kstr in k_sweep:
            vals = xs_for(kstr)
            ok = [i for i, v in enumerate(vals) if v is not None]
            entry = {"k": float(kstr), "n_detected": len(ok)}
            if len(ok) >= 3:
                s, b, r2 = ols([vals[i] for i in ok], y_rec[ok])
                entry.update(slope=s, intercept=b, r2=r2)
            T["k_sweep"].append(entry)

        # The four named criteria, compared with 2011.
        for crit, kstr in CRIT_K.items():
            vals = xs_for(kstr)
            T["crossings"][crit] = vals
            ok = [i for i, v in enumerate(vals) if v is not None]
            entry = {"detected": vals, "thesis": THESIS["samples"][crit],
                     "n_detected": len(ok)}
            if len(ok) >= 3:
                x = np.array([vals[i] for i in ok], float)
                s, b, r2 = ols(x, y_rec[ok])
                entry.update(slope=s, intercept=b, r2=r2)
                ts, tb, tr2 = THESIS["regressions"][crit]
                entry["thesis_regression"] = {"slope": ts, "intercept": tb,
                                              "r2_percent": tr2}
                # Correlation of our detected samples with the thesis's.
                tvals = np.array([THESIS["samples"][crit][i] for i in ok],
                                 float)
                if len(ok) >= 3:
                    entry["corr_with_thesis_samples"] = float(
                        np.corrcoef(x, tvals)[0, 1])
            T["criteria"][crit] = entry

        # Primary confirmatory numbers: Criterion III.
        vals = T["crossings"]["III"]
        ok = [i for i, v in enumerate(vals) if v is not None]
        prim = {"n_detected": len(ok)}
        if len(ok) >= 4:
            x = np.array([vals[i] for i in ok], float)
            y = y_rec[ok]
            s, b, r2 = ols(x, y)
            preds, r2cv, rmse = loocv(x, y)
            boot = bootstrap_ci(x, y)
            prim.update(slope=s, intercept=b, r2=r2,
                        loocv_r2=r2cv, loocv_rmse_samples=rmse,
                        loocv_rmse_ms=rmse / FS * 1000.0,
                        loocv_predictions=preds.tolist(),
                        bootstrap=boot)
            sd, bd, r2d = ols(x, y_dig[ok])
            prim["digitized_gt_r2"] = r2d
            prim["digitized_gt_slope"] = sd
        T["primary_criterion_iii"] = prim

        # Sensitivity table: R2 at k=2 under each one-factor variant.
        sens_defs = rooms[1][sens_key]
        for j, sv in enumerate(sens_defs):
            vals = [rooms[n][sens_key][j]["crossing_k2"] for n in INCLUDED]
            ok = [i for i, v in enumerate(vals) if v is not None]
            entry = {"varied": sv["varied"], "value": sv["value"],
                     "n_detected": len(ok), "detected": vals}
            if len(ok) >= 3:
                s, b, r2 = ols([vals[i] for i in ok], y_rec[ok])
                entry.update(slope=s, r2=r2)
            T["sensitivity"].append(entry)

        # Exploratory measures (pra RIR only in run_room, recorded there).
        if track_key == "pra":
            for meas in ("sampen", "katz"):
                vals = [rooms[n]["exploratory"].get(meas, {}).get(
                        "crossing_k2") for n in INCLUDED]
                ok = [i for i, v in enumerate(vals) if v is not None]
                entry = {"detected": vals, "n_detected": len(ok)}
                if len(ok) >= 3:
                    s, b, r2 = ols([vals[i] for i in ok], y_rec[ok])
                    entry.update(slope=s, intercept=b, r2=r2)
                T["exploratory"][meas] = entry

        out["tracks"][track_key] = T

    (ROOT / "results").mkdir(exist_ok=True)
    with open(ROOT / "results/analysis.json", "w") as f:
        json.dump(out, f, indent=1)

    make_figures(rooms, gt_rooms, out)
    write_summary(out)
    print("wrote results/analysis.json, results/summary.md, figures")
    return out


def make_figures(rooms, gt_rooms, out):
    figdir = ROOT / "results/figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # 1. Per-room smoothed profiles with Criterion III threshold, AB track.
    for suffix, label in (("", "pyroomacoustics"), ("_ab", "Allen-Berkley")):
        fig, axes = plt.subplots(4, 2, figsize=(13, 14), sharey=True)
        for ax, n in zip(axes.ravel(), INCLUDED):
            sm = np.load(ROOT / f"data/profiles/room{n}_smoothed{suffix}.npy")
            res = rooms[n]["ab" if suffix else "primary"]
            ax.plot(sm, lw=0.6, color="tab:blue")
            thr = res["tail_mean"] - 2 * res["tail_sd"]
            ax.axhline(thr, color="tab:red", lw=0.8,
                       label="Criterion III threshold")
            c = res["crossings"].get("2.0")
            if c is not None:
                ax.axvline(c, color="tab:red", ls="--", lw=0.8)
            gt = gt_rooms[n]["tmp50_samples_reconciled"]
            ax.axvline(gt, color="black", ls=":", lw=1.0,
                       label="tmp50 (Lindau, reconciled)")
            ax.set_title(f"Room {n} ({rooms[n]['name']})", fontsize=9)
            ax.set_xlabel("sample")
            ax.set_ylim(0.8, 2.3)
            if n == 1:
                ax.legend(fontsize=7)
        fig.suptitle(f"Smoothed HFD profiles, {label} track")
        fig.tight_layout()
        fig.savefig(figdir / f"profiles{suffix or '_pra'}.png", dpi=130)
        plt.close(fig)

    # 2. Criterion III regression scatter, both tracks.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, tk in zip(axes, ("pra", "ab")):
        T = out["tracks"][tk]
        vals = T["crossings"]["III"]
        y = np.array(out["ground_truth_samples_reconciled"])
        ok = [i for i, v in enumerate(vals) if v is not None]
        x = np.array([vals[i] for i in ok], float)
        yy = y[ok]
        ax.scatter(x, yy, marker="x", color="tab:blue")
        for i, n in enumerate(INCLUDED):
            if i in ok:
                ax.annotate(str(n), (vals[i], y[i]), fontsize=8,
                            xytext=(3, 3), textcoords="offset points")
        p = T["primary_criterion_iii"]
        if "slope" in p:
            xs = np.linspace(x.min(), x.max(), 10)
            ax.plot(xs, p["slope"] * xs + p["intercept"], "r-", lw=1)
            ax.set_title(f"{tk}: R2 = {p['r2']*100:.1f}%  "
                         f"(2011: 93.49%)", fontsize=10)
        ax.set_xlabel("Criterion III detected sample")
        ax.set_ylabel("tmp50 (samples at 44100/s)")
    fig.suptitle("Criterion III regression vs Lindau tmp50 (reconciled)")
    fig.tight_layout()
    fig.savefig(figdir / "criterion3_regression.png", dpi=130)
    plt.close(fig)

    # 3. k sweep curves.
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for tk, color in (("pra", "tab:blue"), ("ab", "tab:orange")):
        ks = [e["k"] for e in out["tracks"][tk]["k_sweep"] if "r2" in e]
        r2 = [e["r2"] * 100 for e in out["tracks"][tk]["k_sweep"]
              if "r2" in e]
        ax.plot(ks, r2, "o-", color=color, label=tk, ms=3)
    for crit, kk in (("I", 0), ("II", 1), ("III", 2), ("IV", 3)):
        ax.axvline(kk, color="gray", lw=0.5, ls=":")
        ax.annotate(crit, (kk, 100), fontsize=8, ha="center")
    ax.axhline(93.49, color="black", lw=0.7, ls="--",
               label="2011 Criterion III (93.49)")
    ax.set_xlabel("k (threshold = tail mean - k * tail SD)")
    ax.set_ylabel("R2 (%)")
    ax.set_ylim(-5, 105)
    ax.legend(fontsize=8)
    ax.set_title("Full threshold sweep, R2 vs k (n = 8 rooms)")
    fig.tight_layout()
    fig.savefig(figdir / "k_sweep.png", dpi=130)
    plt.close(fig)

    # 4. LOOCV predicted vs actual, both tracks.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, tk in zip(axes, ("pra", "ab")):
        p = out["tracks"][tk]["primary_criterion_iii"]
        if "loocv_predictions" not in p:
            continue
        vals = out["tracks"][tk]["crossings"]["III"]
        ok = [i for i, v in enumerate(vals) if v is not None]
        y = np.array(out["ground_truth_samples_reconciled"])[ok]
        preds = np.array(p["loocv_predictions"])
        ax.scatter(y, preds, marker="x")
        for j, i in enumerate(ok):
            ax.annotate(str(INCLUDED[i]), (y[j], preds[j]), fontsize=8,
                        xytext=(3, 3), textcoords="offset points")
        lim = [min(y.min(), preds.min()) * 0.9,
               max(y.max(), preds.max()) * 1.1]
        ax.plot(lim, lim, "k:", lw=0.8)
        ax.set_xlabel("actual tmp50 (samples)")
        ax.set_ylabel("LOOCV predicted (samples)")
        ax.set_title(f"{tk}: LOOCV R2 = {p['loocv_r2']*100:.1f}%, "
                     f"RMSE = {p['loocv_rmse_ms']:.1f} ms", fontsize=10)
    fig.suptitle("Leave-one-out cross-validation, Criterion III")
    fig.tight_layout()
    fig.savefig(figdir / "loocv.png", dpi=130)
    plt.close(fig)


def write_summary(out):
    L = []
    L.append("# Analysis summary (auto-generated by analysis/analyze.py)\n")
    L.append("\n## Thesis baseline recomputed (thesis Table 6.1b samples "
             "vs recovered ground truth)\n\n")
    L.append("| Crit | recomputed R2 | published R2 | recomputed slope | "
             "published slope |\n|---|---|---|---|---|\n")
    for crit in ("I", "II", "III", "IV"):
        e = out["thesis_baseline"][crit]
        pub = e["published"]
        L.append(f"| {crit} | {e['r2']*100:.2f}% | {pub['r2_percent']}% | "
                 f"{e['slope']:.4f} | {pub['slope']} |\n")
    for tk in ("pra", "ab"):
        T = out["tracks"][tk]
        p = T["primary_criterion_iii"]
        L.append(f"\n## Track {tk}\n")
        if "r2" in p:
            L.append(f"Primary Criterion III: R2 = {p['r2']*100:.2f}% "
                     f"(2011: 93.49%), slope = {p['slope']:.4f} "
                     f"(2011: 0.3197), intercept = {p['intercept']:.0f} "
                     f"(2011: 325), n = {p['n_detected']}\n")
            L.append(f"LOOCV: R2 = {p['loocv_r2']*100:.2f}%, "
                     f"RMSE = {p['loocv_rmse_ms']:.2f} ms\n")
            b = p["bootstrap"]
            L.append(f"Bootstrap 95% CI: R2 [{b['r2_ci'][0]*100:.1f}, "
                     f"{b['r2_ci'][1]*100:.1f}]%, slope "
                     f"[{b['slope_ci'][0]:.4f}, {b['slope_ci'][1]:.4f}]\n")
            L.append(f"With raw digitized ground truth: R2 = "
                     f"{p['digitized_gt_r2']*100:.2f}%\n")
        L.append("\nCriteria vs 2011:\n\n")
        L.append("| Crit | our R2 | 2011 R2 | our slope | 2011 slope | "
                 "corr(our x, 2011 x) |\n|---|---|---|---|---|---|\n")
        for crit in ("I", "II", "III", "IV"):
            e = T["criteria"][crit]
            if "r2" in e:
                tr = e["thesis_regression"]
                L.append(f"| {crit} | {e['r2']*100:.2f}% | "
                         f"{tr['r2_percent']}% | {e['slope']:.4f} | "
                         f"{tr['slope']} | "
                         f"{e.get('corr_with_thesis_samples', float('nan')):.3f} |\n")
        L.append("\nSensitivity (Criterion III R2 under one-factor "
                 "changes):\n\n| varied | value | R2 | n |\n|---|---|---|---|\n")
        for e in T["sensitivity"]:
            r2s = f"{e['r2']*100:.2f}%" if "r2" in e else "NA"
            L.append(f"| {e['varied']} | {e['value']} | {r2s} | "
                     f"{e['n_detected']} |\n")
        if T["exploratory"]:
            L.append("\nExploratory (k=2, pra RIRs, separate from the "
                     "primary):\n\n| measure | R2 | n |\n|---|---|---|\n")
            for m, e in T["exploratory"].items():
                r2s = f"{e['r2']*100:.2f}%" if "r2" in e else "NA"
                L.append(f"| {m} | {r2s} | {e['n_detected']} |\n")
    with open(ROOT / "results/summary.md", "w") as f:
        f.write("".join(L))


if __name__ == "__main__":
    main()
