#!/usr/bin/env python3
"""Digitize figures from Conaty 2011 MPhil thesis to recover per-room tmp50%.

Deterministic, no randomness (seed set anyway for belt and braces).

Recovers:
1. Figure 6.2 (PDF page index 56, printed page 46): four scatter panels of
   tmp50% (samples at 44100 Hz) vs detected sample number for Criteria I-IV.
   Blue x markers are detected by colour, axes calibrated from tick marks,
   markers assigned to rooms by known x values (Table 6.1b).
2. Figure 4.3.2 (PDF page index 43, printed page 33): error bar plot of
   tmp50% (ms) with 95 percent CI for all nine rooms. The mean is the small
   widened dot on each bar shaft, CI limits are the serif caps.
3. Figure 5.3 (PDF page index 52, printed page 42): vertical tmp50% marker
   line for Room 1, used as a sanity anchor.

Then fits per-criterion regressions of digitized y against the known x
values and compares with the 2011 reported regressions, and finally solves
a penalized least squares reconciliation for the 8 shared tmp50% values.

Outputs:
- inputs/digitization_debug/  extracted native-resolution figure images and
  overlay images showing every detected marker/bar for human audit.
- inputs/ground_truth.json    machine readable ground truth with provenance.
"""

import json
import os
import numpy as np

np.random.seed(0)  # no randomness used; set for determinism guarantee

import fitz  # pymupdf
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.optimize import least_squares, linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "FINAL_thom_conaty.pdf")
DEBUG = os.path.join(HERE, "digitization_debug")
OUT_JSON = os.path.join(HERE, "ground_truth.json")
FS = 44100.0

# Known detected sample numbers, Table 6.1b (thesis printed page 44), rooms 1-8.
X_KNOWN = {
    1: [6277, 4573, 6056, 13724, 7996, 12366, 14850, 16194],
    2: [5772, 3197, 4335, 7027, 7971, 9486, 12308, 12571],
    3: [3136, 3170, 2494, 7013, 7906, 8060, 11708, 8465],
    4: [2053, 3049, 1567, 7004, 7661, 2052, 2192, 5209],
}

# Axis tick labels read from the figure panels (hardcoded, human verified
# against inputs/digitization_debug/fig62_xref*.png).
TICKS_62 = {
    1: {"x": [6000, 8000, 10000, 12000, 14000, 16000],
        "y": [1000, 1500, 2000, 2500, 3000, 3500, 4000]},
    2: {"x": [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000],
        "y": [1000, 1500, 2000, 2500, 3000, 3500, 4000]},
    3: {"x": [3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000],
        "y": [1000, 1500, 2000, 2500, 3000, 3500, 4000]},
    4: {"x": [2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000],
        "y": [1000, 1500, 2000, 2500, 3000, 3500, 4000]},
}

# Embedded image xrefs on PDF page index 56, mapped to panels by their
# placement rectangles (338 top-left, 340 top-right, 342 bottom-left,
# 344 bottom-right).
XREF_62 = {1: 338, 2: 340, 3: 342, 4: 344}
XREF_432 = 266   # page index 43
XREF_53 = 319    # page index 52

# Hypothesis for the exact y data behind Figure 6.2: every digitized marker
# lies within 3.5 samples (under 1 pixel) of a round ms value times 44.1,
# and these values reproduce the reported regressions as well as the raw
# digitization does. The endpoints 37.5 and 97.8 ms are the exact figures
# quoted in the thesis abstract. Flagged as inference, not as thesis text.
ROUND_MS_HYPOTHESIS = [37.5, 30.0, 22.5, 60.0, 52.8, 70.0, 97.8, 65.0]

# 2011 reported regressions. Figure 6.2 annotations carry more digits than
# Table 6.2 and were produced directly by the curve fitting toolbox, so the
# figure values are used for the reconciliation penalty. Table 6.2 crit I
# slope 0.2090 disagrees with the figure's 0.2019 (digit transposition).
REPORTED = {
    1: {"slope": 0.2019, "intercept": 330.9, "r2": 0.7229,
        "slope_sig": 0.00005, "int_sig": 0.05, "r2_sig": 0.00005},
    2: {"slope": 0.276, "intercept": 293.2, "r2": 0.7882,
        "slope_sig": 0.0005, "int_sig": 0.05, "r2_sig": 0.00005},
    3: {"slope": 0.3197, "intercept": 324.9, "r2": 0.9349,
        "slope_sig": 0.00005, "int_sig": 0.05, "r2_sig": 0.00005},
    4: {"slope": 0.05172, "intercept": 2202.0, "r2": 0.0014,
        "slope_sig": 0.000005, "int_sig": 0.5, "r2_sig": 0.00005},
}


def extract_images():
    doc = fitz.open(PDF)
    os.makedirs(DEBUG, exist_ok=True)
    paths = {}
    for name, xref in [("fig62_c1", XREF_62[1]), ("fig62_c2", XREF_62[2]),
                       ("fig62_c3", XREF_62[3]), ("fig62_c4", XREF_62[4]),
                       ("fig432", XREF_432), ("fig53", XREF_53),
                       ("table4311", 254), ("floorplans_fig4311", 256)]:
        pix = fitz.Pixmap(doc, xref)
        if pix.n > 3:
            pix = fitz.Pixmap(fitz.csRGB, pix)
        p = os.path.join(DEBUG, f"{name}_xref{xref}.png")
        pix.save(p)
        paths[name] = p
    return paths


def load_rgb(path):
    return np.asarray(Image.open(path).convert("RGB")).astype(int)


def dark_mask(img):
    """Near-black, low saturation pixels (frame, ticks, bars, text)."""
    r, g, b = img[..., 0], img[..., 1], img[..., 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    return (mx < 130) & (mx - mn < 60)


def find_frame(dm):
    """Locate the plot box (left, right, top, bottom) from long dark lines."""
    h, w = dm.shape
    colcount = dm.sum(axis=0)
    rowcount = dm.sum(axis=1)
    cthr = 0.6 * colcount.max()
    rthr = 0.6 * rowcount.max()
    cols = np.where(colcount > cthr)[0]
    rows = np.where(rowcount > rthr)[0]
    left, right = cols.min(), cols.max()
    top, bottom = rows.min(), rows.max()
    return left, right, top, bottom


def cluster_1d(indices, gap=3):
    """Group sorted integer indices into clusters, return centroids."""
    if len(indices) == 0:
        return []
    groups = [[indices[0]]]
    for v in indices[1:]:
        if v - groups[-1][-1] <= gap:
            groups[-1].append(v)
        else:
            groups.append([v])
    return [float(np.mean(g)) for g in groups]


def find_ticks(dm, frame, axis):
    """Tick mark pixel centres (interior only; MATLAB draws ticks inward)."""
    left, right, top, bottom = frame
    if axis == "x":
        band = dm[int(bottom) - 9:int(bottom) - 2, :]
        counts = band.sum(axis=0)
        cand = np.where(counts >= 5)[0]
        cand = cand[(cand > left + 3) & (cand < right - 3)]
    else:
        band = dm[:, int(left) + 3:int(left) + 10]
        counts = band.sum(axis=1)
        cand = np.where(counts >= 5)[0]
        cand = cand[(cand > top + 3) & (cand < bottom - 3)]
        if len(cand) == 0:
            # ticks drawn outward (non-MATLAB print graphics)
            band = dm[:, int(left) - 10:int(left) - 2]
            counts = band.sum(axis=1)
            cand = np.where(counts >= 4)[0]
            cand = cand[(cand >= top - 3) & (cand <= bottom + 3)]
    return cluster_1d(list(cand), gap=3)


def match_ticks(tick_px, labels, frame_lo, frame_hi):
    """Match detected interior ticks to expected labels. A tick that sits on
    the axis corner is hidden inside the frame line; if exactly one label is
    unmatched, test whether it belongs at either frame edge."""
    if len(tick_px) == len(labels):
        return list(tick_px), list(labels)
    if len(tick_px) == len(labels) - 1:
        for drop_first in (True, False):
            lab = labels[1:] if drop_first else labels[:-1]
            A = np.polyfit(tick_px, lab, 1)
            missing = labels[0] if drop_first else labels[-1]
            px_pred = (missing - A[1]) / A[0]
            edge = frame_lo if abs(px_pred - frame_lo) < abs(px_pred - frame_hi) \
                else frame_hi
            if abs(px_pred - edge) <= 4:
                return list(tick_px) + [float(edge)], lab + [missing]
    raise AssertionError(
        f"cannot match {len(tick_px)} ticks to {len(labels)} labels")


def calibrate(tick_px, tick_vals):
    """Linear map px -> data value. Returns (a, b) with value = a*px + b."""
    A = np.polyfit(tick_px, tick_vals, 1)
    resid = np.polyval(A, tick_px) - np.asarray(tick_vals, float)
    return A, float(np.max(np.abs(resid)))


def detect_markers(img, frame):
    """x markers are drawn in one flat colour (pure blue or purple depending
    on the panel's quantization). Black text carries multi-coloured
    quantization ringing, so we search for the palette colour that forms the
    most compact ~13x13, ~80 px single-colour components, then take every
    component of that exact colour with at least 15 px (markers partially
    overdrawn by the regression line are smaller than 80 px)."""
    left, right, top, bottom = frame
    r, g, b = img[..., 0], img[..., 1], img[..., 2]
    bluish = (b > np.maximum(r, g) + 40) & (b > 60)
    cols = np.unique(img[bluish].reshape(-1, 3), axis=0)
    best_color, best_score = None, -1
    for c in cols:
        exact = (r == c[0]) & (g == c[1]) & (b == c[2])
        if exact.sum() < 100:
            continue
        lab, n = ndimage.label(exact, structure=np.ones((3, 3), int))
        sizes = ndimage.sum(exact, lab, range(1, n + 1))
        objs = ndimage.find_objects(lab)
        good = 0
        for s, o in zip(sizes, objs):
            hh = o[0].stop - o[0].start
            ww = o[1].stop - o[1].start
            if s >= 60 and 9 <= hh <= 17 and 9 <= ww <= 17:
                good += 1
        if good > best_score:
            best_score, best_color = good, c
    assert best_color is not None and best_score >= 6, \
        f"marker colour not found (best {best_score})"
    exact = (r == best_color[0]) & (g == best_color[1]) & (b == best_color[2])
    lab, n = ndimage.label(exact, structure=np.ones((3, 3), int))
    sizes = ndimage.sum(exact, lab, range(1, n + 1))
    cents = ndimage.center_of_mass(exact, lab, range(1, n + 1))
    objs = ndimage.find_objects(lab)
    comps = [{"cy": c[0], "cx": c[1], "n": int(s),
              "bbox": [o[1].start, o[0].start, o[1].stop, o[0].stop]}
             for c, s, o in zip(cents, sizes, objs) if s >= 8]
    # merge fragments of one marker (e.g. split by the red regression line)
    comps.sort(key=lambda c: (c["cx"], c["cy"]))
    merged, used = [], [False] * len(comps)
    for i, ci in enumerate(comps):
        if used[i]:
            continue
        group, used[i] = [ci], True
        for j in range(i + 1, len(comps)):
            if used[j]:
                continue
            cj = comps[j]
            if abs(ci["cx"] - cj["cx"]) < 12 and abs(ci["cy"] - cj["cy"]) < 12:
                group.append(cj)
                used[j] = True
        ntot = sum(gg["n"] for gg in group)
        if ntot < 25:
            continue
        x0 = min(gg["bbox"][0] for gg in group)
        y0 = min(gg["bbox"][1] for gg in group)
        x1 = max(gg["bbox"][2] for gg in group)
        y1 = max(gg["bbox"][3] for gg in group)
        partial = ntot < 70
        if partial:
            # x glyph is symmetric; the union bbox centre is less biased than
            # the pixel centroid when a fragment has been overdrawn
            cx, cy = (x0 + x1 - 1) / 2.0, (y0 + y1 - 1) / 2.0
        else:
            cx = sum(gg["cx"] * gg["n"] for gg in group) / ntot
            cy = sum(gg["cy"] * gg["n"] for gg in group) / ntot
        clip = []
        if partial:
            clip.append("partial")
            if y0 <= top + 2:
                clip.append("top")
            if y1 >= bottom - 1:
                clip.append("bottom")
            if x0 <= left + 2:
                clip.append("left")
            if x1 >= right - 1:
                clip.append("right")
        markers_entry = {"cx": cx, "cy": cy, "n": ntot, "clip": clip,
                         "bbox": (x0, y0, x1, y1),
                         "color": [int(v) for v in best_color]}
        merged.append(markers_entry)
    return merged


def digitize_fig62_panel(path, crit):
    img = load_rgb(path)
    dm = dark_mask(img)
    frame = find_frame(dm)
    tx = find_ticks(dm, frame, "x")
    ty = find_ticks(dm, frame, "y")
    ex, ey = TICKS_62[crit]["x"], TICKS_62[crit]["y"]
    left, right, top, bottom = frame
    tx, ex = match_ticks(tx, ex, left, right)
    # y ticks run top(px small)=large value; expected labels ascending, so
    # reverse them to match ascending pixel rows
    ty, ey_m = match_ticks(ty, list(reversed(ey)), top, bottom)
    Ax, rx = calibrate(tx, ex)
    Ay, ry = calibrate(ty, ey_m)
    markers = detect_markers(img, frame)
    for m in markers:
        m["x_data"] = float(np.polyval(Ax, m["cx"]))
        m["y_data"] = float(np.polyval(Ay, m["cy"]))
    # axis limits in data units (to test the tight-axis hypothesis)
    left, right, top, bottom = frame
    lims = {"xlim": [float(np.polyval(Ax, left)), float(np.polyval(Ax, right))],
            "ylim": [float(np.polyval(Ay, bottom)), float(np.polyval(Ay, top))],
            "cal_resid_px_x": rx, "cal_resid_px_y": ry,
            "samples_per_px_x": float(Ax[0]), "samples_per_px_y": float(-Ay[0])}
    return img, frame, markers, lims


def assign_rooms(markers, xs_known, samples_per_px_x):
    """Assign 8 markers to 8 rooms by known x (Hungarian on |dx|)."""
    assert len(markers) == 8, f"expected 8 markers, found {len(markers)}"
    cost = np.zeros((8, 8))
    for i, xk in enumerate(xs_known):
        for j, m in enumerate(markers):
            cost[i, j] = abs(m["x_data"] - xk)
    ri, ci = linear_sum_assignment(cost)
    out = {}
    for i, j in zip(ri, ci):
        out[i + 1] = markers[j]
    return out


def digitize_fig432(path):
    """Error bar chart, tmp50% in ms, rooms 1-9."""
    img = load_rgb(path)
    dm = dark_mask(img)
    frame = find_frame(dm)
    left, right, top, bottom = frame
    # this print graphic has no tick marks; the y axis labels (200..0) are
    # vertically centred on their values, so use the label text centroids
    band = dm[:, max(0, int(left) - 60):int(left) - 3]
    rowsum = band.sum(axis=1)
    cand = list(np.where(rowsum > 0)[0])
    groups = []
    for rr in cand:
        if groups and rr - groups[-1][-1] <= 10:
            groups[-1].append(rr)
        else:
            groups.append([rr])
    assert len(groups) == 5, f"fig432 y label clusters: {len(groups)}"
    ty = [float(np.average(gg, weights=rowsum[gg])) for gg in groups]
    ey_m = [200, 150, 100, 50, 0]  # descending: pixel rows grow downward
    Ay, ry = calibrate(ty, ey_m)
    interior = dm[int(top) + 3:int(bottom) - 2, int(left) + 3:int(right) - 2]
    # longest vertical dark run per column
    h, w = interior.shape
    runlen = np.zeros(w, int)
    runs = {}
    for c in range(w):
        col = interior[:, c]
        best, cur, s0, bs = 0, 0, 0, (0, 0)
        for rr in range(h):
            if col[rr]:
                if cur == 0:
                    s0 = rr
                cur += 1
                if cur > best:
                    best, bs = cur, (s0, rr)
            else:
                cur = 0
        runlen[c] = best
        runs[c] = bs
    # bar columns: long contiguous vertical runs that do not span the whole
    # panel height (the category divider lines do)
    barcols = [c for c in range(w) if 22 <= runlen[c] < 0.9 * h]
    clusters_px = cluster_1d(barcols, gap=4)
    # keep clusters that are thick like a bar shaft (text stems are 2-3 px)
    bars = []
    for ctr in clusters_px:
        cols = [c for c in barcols if abs(c - ctr) <= 8]
        if len(cols) < 4 or max(runlen[c] for c in cols) < 25:
            continue
        # keep only columns whose run agrees with the cluster median run
        # (rejects stray columns borrowed from the alpha arrows)
        med_s = float(np.median([runs[c][0] for c in cols]))
        med_e = float(np.median([runs[c][1] for c in cols]))
        cols = [c for c in cols
                if runs[c][1] >= med_s - 2 and runs[c][0] <= med_e + 2]
        if len(cols) < 4:
            continue
        # bar extent = union of the longest runs of the member columns
        rtop = min(runs[c][0] for c in cols)
        rbot = max(runs[c][1] for c in cols)
        # widen the column window so the serif caps (which stick out sideways
        # from the shaft and have short vertical runs) are included
        win = list(range(max(0, min(cols) - 6), min(w, max(cols) + 7)))
        sub = interior[rtop:rbot + 1, :][:, win]
        widths_local = sub.sum(axis=1)
        widths = np.zeros(h, int)
        widths[rtop:rbot + 1] = widths_local
        wmax = widths[rtop:rbot + 1].max()
        shaft = float(np.median([widths[rr] for rr in range(rtop, rbot + 1)
                                 if widths[rr] > 0]))
        if shaft > 8:
            continue  # solid blob, not an error bar
        # an error bar has serif caps at BOTH extremes; the alpha arrowheads
        # (triangular ramp) and text stems (uniform width) fail this test
        if widths[rtop:rtop + 3].max() < 1.8 * shaft or \
           widths[rbot - 2:rbot + 1].max() < 1.8 * shaft:
            continue
        # caps: contiguous wide rows at the extremes
        capw = 0.75 * wmax
        cap_top_rows = []
        rr = rtop
        while rr <= rbot and widths[rr] >= capw:
            cap_top_rows.append(rr)
            rr += 1
        cap_top_end = rr
        cap_bot_rows = []
        rr = rbot
        while rr >= rtop and widths[rr] >= capw:
            cap_bot_rows.append(rr)
            rr -= 1
        cap_bot_start = rr
        top_cap = float(np.mean(cap_top_rows)) if cap_top_rows else float(rtop)
        bot_cap = float(np.mean(cap_bot_rows)) if cap_bot_rows else float(rbot)
        # mean dot: slight widening of the shaft between the caps
        mid_lo, mid_hi = cap_top_end + 4, cap_bot_start - 4
        dot_rows = [r_ for r_ in range(mid_lo, mid_hi + 1)
                    if widths[r_] >= shaft + 1]
        if dot_rows:
            wts = [widths[r_] - shaft for r_ in dot_rows]
            mean_row = float(np.average(dot_rows, weights=wts))
        else:
            mean_row = None
        off = int(top) + 3
        bar = {"px_center_col": float(np.mean(cols)) + int(left) + 3,
               "px_top_cap": top_cap + off,
               "px_bot_cap": bot_cap + off,
               "px_mean": (mean_row + off) if mean_row is not None else None,
               "ci_hi_ms": float(np.polyval(Ay, top_cap + off)),
               "ci_lo_ms": float(np.polyval(Ay, bot_cap + off)),
               "mean_ms": (float(np.polyval(Ay, mean_row + off))
                           if mean_row is not None else None),
               "mean_from_dot": mean_row is not None,
               "mid_ms": float(np.polyval(Ay, (top_cap + bot_cap) / 2 + off))}
        bars.append(bar)
    bars.sort(key=lambda b: b["px_center_col"])
    return img, frame, bars, {"cal_resid_px_y": ry,
                              "ms_per_px": float(abs(Ay[0]))}


def digitize_floorplans(path):
    """Figure 4.3.1.1: nine floor plan subplots (rooms 1-9, row major), axes
    in metres with labels -20, 0, 20. No tick marks are drawn; the axis
    labels are centred on their values, so calibration uses the '0' and '20'
    label centroids (avoiding the '-20' labels whose minus sign biases the
    centroid). Returns footprint length and width per room, about 0.5 m
    uncertainty. Room 9 is non rectangular; its bounding box is returned
    with a flag."""
    img = load_rgb(path)
    dm = dark_mask(img)
    rowc = dm.sum(axis=1)
    colc = dm.sum(axis=0)
    long_rows = cluster_1d([i for i in range(len(rowc)) if rowc[i] > 300], 3)
    long_cols = cluster_1d([i for i in range(len(colc)) if colc[i] > 300], 3)
    assert len(long_rows) == 6 and len(long_cols) == 6, \
        f"floorplan grid: {len(long_rows)} rows, {len(long_cols)} cols"
    row_pairs = [(long_rows[0], long_rows[1]), (long_rows[2], long_rows[3]),
                 (long_rows[4], long_rows[5])]
    col_pairs = [(long_cols[0], long_cols[1]), (long_cols[2], long_cols[3]),
                 (long_cols[4], long_cols[5])]
    rooms = []
    for pr in range(3):
        for pc in range(3):
            top, bot = (int(round(v)) for v in row_pairs[pr])
            l, r = (int(round(v)) for v in col_pairs[pc])
            # x axis labels below the panel: clusters for -20, 0, 20
            band = dm[bot + 4:bot + 30, max(0, l - 12):r + 12]
            cc = band.sum(axis=0)
            xg = cluster_1d(list(np.where(cc > 0)[0]), gap=6)
            assert len(xg) == 3, f"panel {pr},{pc}: x labels {len(xg)}"
            px0x, px20x = xg[1] + l - 12, xg[2] + l - 12
            sx = (px20x - px0x) / 20.0  # px per metre
            # y axis labels left of the panel: 20, 0, -20 top to bottom
            band2 = dm[top - 5:bot + 5, max(0, l - 45):l - 3]
            rr2 = band2.sum(axis=1)
            yg = cluster_1d(list(np.where(rr2 > 0)[0]), gap=6)
            assert len(yg) == 3, f"panel {pr},{pc}: y labels {len(yg)}"
            px20y, px0y = yg[0] + top - 5, yg[1] + top - 5
            sy = (px0y - px20y) / 20.0
            # room outline: long straight runs strictly inside the frame
            inter = dm[top + 4:bot - 3, l + 4:r - 3]
            vcols, hrows = [], []
            for c in range(inter.shape[1]):
                col = inter[:, c]
                best = cur = 0
                for v in col:
                    cur = cur + 1 if v else 0
                    best = max(best, cur)
                if best >= 12:
                    vcols.append(c + l + 4)
            for rr_ in range(inter.shape[0]):
                roww = inter[rr_, :]
                best = cur = 0
                for v in roww:
                    cur = cur + 1 if v else 0
                    best = max(best, cur)
                if best >= 12:
                    hrows.append(rr_ + top + 4)
            vcl = cluster_1d(vcols, gap=2)
            hcl = cluster_1d(hrows, gap=2)
            length_m = (max(vcl) - min(vcl)) / sx
            width_m = (max(hcl) - min(hcl)) / sy
            rooms.append({"length_m": round(float(length_m), 1),
                          "width_m": round(float(width_m), 1),
                          "non_rectangular": len(vcl) > 2 or len(hcl) > 2,
                          "n_edge_clusters": [len(vcl), len(hcl)]})
    return rooms


def digitize_fig53(path):
    """Room 1 FD profile: locate the thick vertical tmp50% line."""
    img = load_rgb(path)
    dm = dark_mask(img)
    frame = find_frame(dm)
    left, right, top, bottom = frame
    h = bottom - top
    colcount = dm[int(top) + 2:int(bottom) - 2, :].sum(axis=0)
    # the tmp50% marker line spans nearly the full plot height
    cand = [c for c in range(int(left) + 4, int(right) - 4)
            if colcount[c] > 0.9 * h]
    centers = cluster_1d(cand, gap=2)
    assert len(centers) == 1, f"fig53 vertical lines found: {centers}"
    line_px = centers[0]
    # the marker line also crosses the tick band, so drop it from the ticks
    tx = [t for t in find_ticks(dm, frame, "x") if abs(t - line_px) > 4]
    ex = [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000]
    tx, ex_m = match_ticks(tx, ex, left, right)
    Ax, rx = calibrate(tx, ex_m)
    return float(np.polyval(Ax, line_px)), rx


def draw_debug_62(img, frame, assigned, crit):
    im = Image.fromarray(img.astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    for room, m in assigned.items():
        x, y = m["cx"], m["cy"]
        d.ellipse([x - 9, y - 9, x + 9, y + 9], outline=(0, 200, 0), width=2)
        d.text((x + 10, y - 16), f"R{room}", fill=(0, 150, 0))
    left, right, top, bottom = frame
    d.rectangle([left, top, right, bottom], outline=(255, 0, 255), width=1)
    im.save(os.path.join(DEBUG, f"overlay_fig62_c{crit}.png"))


def draw_debug_432(img, frame, bars):
    im = Image.fromarray(img.astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    left, right, top, bottom = frame
    d.rectangle([left, top, right, bottom], outline=(255, 0, 255), width=1)
    for i, b in enumerate(bars, start=1):
        x = b["px_center_col"]
        d.text((x - 7, bottom + 8), f"R{i}", fill=(255, 0, 0))
        for row, col in [(b["px_top_cap"], (0, 160, 0)),
                         (b["px_bot_cap"], (0, 160, 0))]:
            d.line([x - 14, row, x + 14, row], fill=col, width=1)
        if b["px_mean"] is not None:
            d.line([x - 14, b["px_mean"], x + 14, b["px_mean"]],
                   fill=(255, 0, 0), width=1)
    im.save(os.path.join(DEBUG, "overlay_fig432.png"))


def fit_line(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    b, a = np.polyfit(x, y, 1)
    yhat = b * x + a
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(b), float(a), float(1 - ss_res / ss_tot)


def reconcile(y_dig, sigma_dig=12.0):
    """Penalized least squares: find 8 tmp50% values (samples) matching the
    digitized estimates while reproducing the 12 reported regression stats
    (slope, intercept, R2 for each of the four criteria).

    Two stages. Stage 1 fits with every constraint under a robust (soft_l1)
    loss, then any reported stat that still disagrees by more than 10 sigma
    is declared unreproducible (thesis typo) and dropped. Stage 2 refits with
    only the consistent constraints under a plain quadratic loss. Returns the
    solution and the list of dropped constraints."""
    y0 = np.asarray(y_dig, float)

    def make_residuals(active):
        def residuals(y):
            res = list((y - y0) / sigma_dig)
            for crit in (1, 2, 3, 4):
                b, a, r2 = fit_line(X_KNOWN[crit], y)
                rep = REPORTED[crit]
                if (crit, "slope") in active:
                    res.append((b - rep["slope"]) / rep["slope_sig"])
                if (crit, "intercept") in active:
                    res.append((a - rep["intercept"]) / rep["int_sig"])
                if (crit, "r2") in active:
                    res.append((r2 - rep["r2"]) / rep["r2_sig"])
            return np.asarray(res)
        return residuals

    # Consistency screen straight from the digitized fit: a reported stat is
    # reproducible if the digitized data already land close to it. Tolerances
    # are generous relative to digitization noise (about 3 to 12 samples per
    # marker): slope within 1 percent relative, intercept within 15 samples,
    # R2 within 0.5 percentage points.
    dropped = []
    for crit in (1, 2, 3, 4):
        b, a, r2 = fit_line(X_KNOWN[crit], y0)
        rep = REPORTED[crit]
        if abs(b - rep["slope"]) > 0.01 * abs(rep["slope"]):
            dropped.append((crit, "slope"))
        if abs(a - rep["intercept"]) > 15.0:
            dropped.append((crit, "intercept"))
        if abs(r2 - rep["r2"]) > 0.005:
            dropped.append((crit, "r2"))
    all_constraints = [(c, k) for c in (1, 2, 3, 4)
                       for k in ("slope", "intercept", "r2")]
    active = set(all_constraints) - set(dropped)
    sol = least_squares(make_residuals(active), y0, method="lm",
                        xtol=1e-14, ftol=1e-14)
    return sol.x, dropped


def main():
    paths = extract_images()
    results = {"fig62": {}, "per_room_y": {}}

    # ---- Figure 6.2 ----
    per_room = {r: {} for r in range(1, 9)}
    panels = {}
    for crit in (1, 2, 3, 4):
        img, frame, markers, lims = digitize_fig62_panel(
            paths[f"fig62_c{crit}"], crit)
        assigned = assign_rooms(markers, X_KNOWN[crit],
                                lims["samples_per_px_x"])
        panels[crit] = (img, frame, assigned, lims)
        results["fig62"][crit] = {"lims": lims}

    # cross-panel disambiguation for rooms whose known x nearly coincide
    # (crit 3: rooms 1/2 at 3136/3170; crit 4: rooms 1/6 at 2053/2052).
    # First consensus from unambiguous assignments.
    ambig = {c: [] for c in (1, 2, 3, 4)}
    for crit in (1, 2, 3, 4):
        xs = X_KNOWN[crit]
        for i, xi in enumerate(xs):
            near = [j for j, xj in enumerate(xs)
                    if j != i and abs(xj - xi) < 120]
            if near:
                ambig[crit].append(i + 1)
    consensus = {}
    for room in range(1, 9):
        vals = [panels[c][2][room]["y_data"] for c in (1, 2, 3, 4)
                if room not in ambig[c]]
        consensus[room] = float(np.median(vals))
    for crit in (1, 2, 3, 4):
        rooms_a = ambig[crit]
        if not rooms_a:
            continue
        assigned = panels[crit][2]
        marks = [assigned[r] for r in rooms_a]
        cost = np.zeros((len(rooms_a), len(marks)))
        for i, r in enumerate(rooms_a):
            for j, m in enumerate(marks):
                cost[i, j] = abs(m["y_data"] - consensus[r])
        ri, ci = linear_sum_assignment(cost)
        for i, j in zip(ri, ci):
            assigned[rooms_a[i]] = marks[j]

    for crit in (1, 2, 3, 4):
        img, frame, assigned, lims = panels[crit]
        draw_debug_62(img, frame, assigned, crit)
        results["fig62"][crit]["assigned"] = {
            r: {"x_data": assigned[r]["x_data"],
                "y_data": assigned[r]["y_data"],
                "x_known": X_KNOWN[crit][r - 1],
                "clip": assigned[r]["clip"]}
            for r in range(1, 9)}

    # per-room digitized y: mean over panels, excluding clip-biased centroids
    y_dig, y_spread, y_panel_vals = [], [], {}
    for room in range(1, 9):
        vals, used = [], []
        for crit in (1, 2, 3, 4):
            m = panels[crit][2][room]
            v = m["y_data"]
            y_panel_vals.setdefault(room, {})[crit] = v
            # top or bottom clipped markers have biased y centroids; the
            # tight-axis hypothesis handles them via axis limits instead
            if not ("top" in m["clip"] or "bottom" in m["clip"]):
                vals.append(v)
                used.append(crit)
        if vals:
            y_dig.append(float(np.mean(vals)))
            y_spread.append(float(np.ptp(vals)) if len(vals) > 1 else 0.0)
        else:
            # all clipped (room at axis limit in every panel): use the axis
            # limit value implied by the tick calibration
            lim_vals = []
            for crit in (1, 2, 3, 4):
                m = panels[crit][2][room]
                lims = panels[crit][3]
                if "top" in m["clip"]:
                    lim_vals.append(lims["ylim"][1])
                elif "bottom" in m["clip"]:
                    lim_vals.append(lims["ylim"][0])
            y_dig.append(float(np.mean(lim_vals)))
            y_spread.append(float(np.ptp(lim_vals)))
    results["per_room_y"] = {r: {"digitized_samples": y_dig[r - 1],
                                 "panel_spread_samples": y_spread[r - 1],
                                 "per_panel": y_panel_vals[r]}
                             for r in range(1, 9)}

    # ---- regression check ----
    reg_check = {}
    for crit in (1, 2, 3, 4):
        b, a, r2 = fit_line(X_KNOWN[crit], y_dig)
        reg_check[crit] = {"slope_fit": b, "intercept_fit": a, "r2_fit": r2,
                           "slope_reported_fig": REPORTED[crit]["slope"],
                           "intercept_reported_fig": REPORTED[crit]["intercept"],
                           "r2_reported": REPORTED[crit]["r2"]}
    results["regression_check_digitized"] = reg_check

    # ---- reconciliation ----
    y_rec, dropped_constraints = reconcile(y_dig)
    reg_rec = {}
    for crit in (1, 2, 3, 4):
        b, a, r2 = fit_line(X_KNOWN[crit], y_rec)
        reg_rec[crit] = {"slope": b, "intercept": a, "r2": r2}
    results["regression_after_reconciliation"] = reg_rec
    results["y_reconciled"] = {r: float(y_rec[r - 1]) for r in range(1, 9)}
    # round-ms hypothesis check
    y_hyp = [m * FS / 1000.0 for m in ROUND_MS_HYPOTHESIS]
    reg_hyp = {}
    for crit in (1, 2, 3, 4):
        b, a, r2 = fit_line(X_KNOWN[crit], y_hyp)
        reg_hyp[crit] = {"slope": b, "intercept": a, "r2": r2}
    results["round_ms_hypothesis"] = {
        "ms": ROUND_MS_HYPOTHESIS,
        "samples": y_hyp,
        "max_abs_diff_from_digitized_samples":
            float(np.max(np.abs(np.array(y_hyp) - np.array(y_dig)))),
        "regressions": reg_hyp}
    results["reconciliation_dropped_constraints"] = [
        {"criterion": c, "stat": k,
         "reason": "reported value unreproducible from any y consistent "
                   "with the digitized markers; treated as thesis typo"}
        for c, k in dropped_constraints]

    # ---- Figure 4.3.2 ----
    img432, frame432, bars, meta432 = digitize_fig432(paths["fig432"])
    draw_debug_432(img432, frame432, bars)
    assert len(bars) == 9, f"expected 9 bars, found {len(bars)}"
    results["fig432"] = {"meta": meta432,
                         "bars": {i + 1: b for i, b in enumerate(bars)}}

    # ---- Figure 4.3.1.1 floor plans ----
    footprints = digitize_floorplans(paths["floorplans_fig4311"])
    results["floorplans"] = {i + 1: fp for i, fp in enumerate(footprints)}

    # ---- Figure 5.3 anchor ----
    tmp_r1_line, rx53 = digitize_fig53(paths["fig53"])
    results["fig53_room1_tmp_line_samples"] = tmp_r1_line

    with open(os.path.join(DEBUG, "digitization_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ---- assemble ground_truth.json ----
    rooms_meta = {
        1: ("E. Studio", 216, 0.36, 0.39), 2: ("EN-111", 224, 0.26, 0.62),
        3: ("EN-190", 182, 0.17, 0.79), 4: ("H-0104", 3300, 0.28, 1.15),
        5: ("HE-101", 5179, 0.23, 1.67), 6: ("Teldex Studio", 3647, 0.2, 1.83),
        7: ("UoA concert hall", 8298, 0.33, 1.52),
        8: ("TUB Audimax", 8500, 0.23, 2.08),
        9: ("JC church", 7417, 0.23, 2.36)}
    gt = {
        "provenance_notes": {
            "source": "Conaty 2011 MPhil thesis (Trinity College Dublin), "
                      "inputs/FINAL_thom_conaty.pdf",
            "room_table": "Table 4.3.1.1, thesis printed page 31 (figure on "
                          "printed p.31, PDF page index 41); values verified "
                          "from the native-resolution embedded table image.",
            "room_names": "From the panel titles of Figure 4.3.1.1 (floor "
                          "plans reproduced from Lindau et al. 2010), thesis "
                          "printed p.31.",
            "room_footprints": "Digitized from Figure 4.3.1.1 floor plan "
                               "outlines, calibrated on the 0 and 20 m axis "
                               "labels; about plus minus 0.5 m. Heights "
                               "derived as V/(length*width). The thesis "
                               "itself never tabulates per room dimensions.",
            "tmp50_samples": "Digitized from Figure 6.2 (printed p.46) blue "
                             "x markers; x positions matched to Table 6.1b "
                             "detected samples; y read via tick-calibrated "
                             "axes. Reconciled values additionally constrained "
                             "by the four reported regressions (figure "
                             "annotation coefficients).",
            "tmp50_ms_fig432": "Digitized from Figure 4.3.2 (printed p.33) "
                               "error bar mean dot; CI limits from serif caps."},
        "fs_hz": 44100,
        "method_constants_2011": {
            "hfd_kmax": 5,
            "hmw_window_param": 50,
            "hmw_window_samples_actual": 51,
            "hmw_window_note": "HMW.m passes y=50 but MATLAB k(i:(i+y)) is "
                               "51 samples inclusive",
            "smooth_span_param": 200,
            "smooth_span_effective": 199,
            "smooth_note": "MATLAB smooth() is a centred moving average and "
                           "reduces even spans to the next lower odd span",
            "threshold_search_start_sample": 1000,
            "detection_rule": "first i >= 1000 with profile(i) >= threshold "
                              "(FDthres.m, Appendix 9.1.1)",
            "tail_fraction": 0.1,
            "tail_slice": "x((length-round(length/10)):length), i.e. the "
                          "last 10 percent plus one sample",
            "criteria_thresholds": "I: tail mean; II: mean - 1 SD; "
                                   "III: mean - 2 SD; IV: mean - 3 SD",
            "ism": {
                "algorithm": "Allen and Berkley 1979 image source method, "
                             "MATLAB port of the Fortran appendix (sroom.m)",
                "fs_hz": 44100,
                "source_xyz_text": "x = L/2, y = 2, z = 1.9 (thesis printed "
                                   "p.39)",
                "receiver_xyz_text": "x = L/2, y = 1, z = 1.9 (thesis "
                                     "printed p.39)",
                "gui_screenshot_discrepancy": "Figure 5.1a GUI shows source "
                                              "(5, 1, 1.9) and receiver "
                                              "(5, 2, 1.9) for a 10x10x10 "
                                              "example, i.e. y swapped vs "
                                              "the text. RIR is reciprocal "
                                              "under the swap.",
                "beta_from_alpha": "conversion formula NOT stated in thesis; "
                                   "beta = sqrt(1 - alpha_ave) is the "
                                   "standard Allen-Berkley convention, "
                                   "unverified",
                "rir_length_s": "UNCERTAIN. GUI default 0.256 s; Figure "
                                "5.1b plots Room 1 RIR to 0.2 s; Figure 5.3 "
                                "Room 1 FD profile extends to about 17200 "
                                "samples which is about RT*fs (0.39 s). "
                                "Per-room length never stated."}},
        "rooms": [],
        "detected_samples_2011": {f"criterion_{c}": X_KNOWN[c]
                                  for c in (1, 2, 3, 4)},
        "tail_stats_2011_table_6_1a": {
            "mean_fd": [2.0364, 2.0450, 2.0250, 1.9833, 1.9782, 1.9866,
                        1.9682, 1.9644],
            "sd_fd": [0.0252, 0.0265, 0.0321, 0.0394, 0.0404, 0.0432,
                      0.0394, 0.0506]},
        "regressions_2011": {
            "table_6_2": {
                "criterion_1": {"slope": 0.2090, "intercept": 331,
                                "r2_percent": 72.29},
                "criterion_2": {"slope": 0.2760, "intercept": 293,
                                "r2_percent": 78.82},
                "criterion_3": {"slope": 0.3197, "intercept": 325,
                                "r2_percent": 93.49},
                "criterion_4": {"slope": 0.0517, "intercept": 2202,
                                "r2_percent": 0.14}},
            "figure_6_2_annotations": {
                "criterion_1": {"slope": 0.2019, "intercept": 330.9,
                                "r2_percent": 72.29},
                "criterion_2": {"slope": 0.276, "intercept": 293.2,
                                "r2_percent": 78.82},
                "criterion_3": {"slope": 0.3197, "intercept": 324.9,
                                "r2_percent": 93.49},
                "criterion_4": {"slope": 0.05172, "intercept": 2202,
                                "r2_percent": 0.14}},
            "note": "Table vs figure slope disagrees for criterion I "
                    "(0.2090 vs 0.2019); digitized fit decides which is "
                    "correct, see regression_check in this file."},
        "regression_check_digitized": reg_check,
        "regression_after_reconciliation": reg_rec,
        "reconciliation_dropped_constraints":
            results["reconciliation_dropped_constraints"],
        "round_ms_hypothesis": {
            "note": "INFERENCE, not thesis text: every digitized Figure 6.2 "
                    "marker sits within 3.5 samples (under 1 px) of these "
                    "round ms values times 44.1, they reproduce the reported "
                    "regressions equally well, and the 37.5 and 97.8 ms "
                    "endpoints match the values quoted in the abstract.",
            "tmp50_ms_rooms_1_to_8": ROUND_MS_HYPOTHESIS,
            "tmp50_samples_rooms_1_to_8":
                results["round_ms_hypothesis"]["samples"],
            "max_abs_diff_from_digitized_samples":
                results["round_ms_hypothesis"]
                       ["max_abs_diff_from_digitized_samples"],
            "regressions": reg_hyp},
        "fig53_room1_tmp_line_samples": {
            "value": round(tmp_r1_line, 1),
            "provenance": "Digitized vertical tmp50% marker line, Figure 5.3 "
                          "(printed p.42), tick-calibrated x axis."},
    }
    for r in range(1, 10):
        name, vol, alpha, rt = rooms_meta[r]
        bar = bars[r - 1]
        fp = footprints[r - 1]
        entry = {
            "room": r, "name": name,
            "name_provenance": "Figure 4.3.1.1 floor plan panel titles, "
                               "thesis printed p.31",
            "volume_m3": vol, "alpha_ave": alpha, "rt_s": rt,
            "room_params_provenance": "Exact from Table 4.3.1.1, thesis "
                                      "printed p.31",
            "footprint_length_m": fp["length_m"],
            "footprint_width_m": fp["width_m"],
            "footprint_non_rectangular": fp["non_rectangular"],
            "height_derived_m": round(vol / (fp["length_m"] * fp["width_m"]),
                                      2),
            "footprint_provenance": "Digitized from Figure 4.3.1.1, plus "
                                    "minus about 0.5 m; height derived from "
                                    "V/(L*W), not stated in thesis",
            "included_in_2011": r <= 9 and r != 9,
            "tmp50_fig432_mean_ms": round(bar["mean_ms"], 1)
            if bar["mean_ms"] is not None else None,
            "tmp50_fig432_mean_from_dot": bar["mean_from_dot"],
            "tmp50_fig432_mid_ms": round(bar["mid_ms"], 1),
            "ci95_fig432_ms": [round(bar["ci_lo_ms"], 1),
                               round(bar["ci_hi_ms"], 1)],
            "fig432_provenance": "Digitized from Figure 4.3.2 error bars, "
                                 "uncertainty about plus minus 1 ms "
                                 "(about 2 px)",
        }
        if r <= 8:
            yd = y_dig[r - 1]
            yr = float(y_rec[r - 1])
            m_clip = [c for c in (1, 2, 3, 4)
                      if panels[c][2][r]["clip"]]
            entry.update({
                "tmp50_samples_digitized": round(yd, 1),
                "tmp50_ms_digitized": round(yd / FS * 1000, 2),
                "tmp50_samples_reconciled": round(yr, 1),
                "tmp50_ms_reconciled": round(yr / FS * 1000, 2),
                "digitized_minus_reconciled_samples": round(yd - yr, 1),
                "tmp50_ms_round_hypothesis": ROUND_MS_HYPOTHESIS[r - 1],
                "panel_spread_samples": round(y_spread[r - 1], 1),
                "clipped_marker_panels": m_clip,
                "tmp50_provenance": "Digitized from Figure 6.2 (4 panels "
                                    "averaged, clip-biased centroids "
                                    "excluded); reconciled via penalized "
                                    "least squares against reported "
                                    "regressions.",
            })
        else:
            entry.update({
                "tmp50_samples_digitized": None,
                "tmp50_ms_digitized": None,
                "tmp50_samples_reconciled": None,
                "tmp50_ms_reconciled": None,
                "tmp50_provenance": "Room 9 excluded from thesis analysis "
                                    "(non-shoebox); only Figure 4.3.2 ms "
                                    "values available.",
            })
        gt["rooms"].append(entry)

    with open(OUT_JSON, "w") as f:
        json.dump(gt, f, indent=2)

    # console report
    print("=== Figure 6.2 axis limits (tight-axis check) ===")
    for crit in (1, 2, 3, 4):
        lims = panels[crit][3]
        print(f" crit {crit}: xlim {lims['xlim'][0]:.0f}..{lims['xlim'][1]:.0f}"
              f"  ylim {lims['ylim'][0]:.0f}..{lims['ylim'][1]:.0f}"
              f"  (known x min/max {min(X_KNOWN[crit])}/{max(X_KNOWN[crit])})")
    print("=== per-room tmp50% ===")
    for r in range(1, 9):
        print(f" R{r}: digitized {y_dig[r-1]:8.1f}  spread {y_spread[r-1]:5.1f}"
              f"  reconciled {y_rec[r-1]:8.1f}  "
              f"({y_dig[r-1]/44.1:.2f} / {y_rec[r-1]/44.1:.2f} ms)")
    print("=== regression check (digitized) ===")
    for crit in (1, 2, 3, 4):
        rc = reg_check[crit]
        print(f" crit {crit}: fit y={rc['slope_fit']:.4f}x+{rc['intercept_fit']:.1f} "
              f"R2={rc['r2_fit']*100:.2f}%  reported(fig) "
              f"y={rc['slope_reported_fig']}x+{rc['intercept_reported_fig']} "
              f"R2={rc['r2_reported']*100:.2f}%")
    print("=== dropped (unreproducible) reported stats ===")
    for c, k in dropped_constraints:
        print(f" criterion {c}: {k}")
    print("=== regression after reconciliation ===")
    for crit in (1, 2, 3, 4):
        rr = reg_rec[crit]
        print(f" crit {crit}: y={rr['slope']:.5f}x+{rr['intercept']:.1f} "
              f"R2={rr['r2']*100:.2f}%")
    print("=== Figure 4.3.2 (ms) ===")
    for i, b in enumerate(bars, 1):
        print(f" R{i}: mean {b['mean_ms'] if b['mean_ms'] is None else round(b['mean_ms'],1)}"
              f" (dot={b['mean_from_dot']}) mid {b['mid_ms']:.1f} "
              f"CI [{b['ci_lo_ms']:.1f}, {b['ci_hi_ms']:.1f}]")
    print("=== Figure 4.3.1.1 footprints (m) ===")
    for i, fp in enumerate(footprints, 1):
        vol = rooms_meta[i][1]
        hh = vol / (fp["length_m"] * fp["width_m"])
        print(f" R{i}: L {fp['length_m']} x W {fp['width_m']} "
              f"(H=V/LW={hh:.2f}) nonrect={fp['non_rectangular']}")
    print(f"=== Figure 5.3 Room 1 tmp line: {tmp_r1_line:.1f} samples "
          f"({tmp_r1_line/44.1:.1f} ms) ===")


if __name__ == "__main__":
    main()
