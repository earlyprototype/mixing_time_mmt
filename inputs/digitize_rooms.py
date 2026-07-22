"""Digitization of Figure 4.3.1.1 (thesis page 42, embedded raster 723x755).

This script produced the L and W values hard-coded in src/mixtime/rooms.py.
Calibration: dotted gridlines at data -20, 0, +20 on both axes of each
subplot, detected as near-full-span lines at gray threshold < 210.
Room outlines are solid black strokes (gray < 128, run >= 8 px).
Edge positions are stroke-centre centroids, converted with per-axis scale.
Deterministic: pure image processing, no randomness.

Run:  python inputs/digitize_rooms.py
Writes digitized_rooms.json next to this file.
"""
import io
import json
import os

import fitz
import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation

# Extract the native raster of Figure 4.3.1.1 (image xref 256 on PDF page 42,
# zero-based index 41) straight from the thesis PDF.
_here = os.path.dirname(os.path.abspath(__file__))
_doc = fitz.open(os.path.join(_here, 'FINAL_thom_conaty.pdf'))
_pix = fitz.Pixmap(_doc, 256)
img = np.array(Image.open(io.BytesIO(_pix.tobytes('png'))).convert('L'))
assert img.shape == (755, 723), img.shape
dark = img < 128
gray = img < 210

frame_x = [(80.5, 246.5), (298.5, 464.5), (516.5, 682.5)]
frame_y = [(44.5, 212.5), (277.5, 445.0), (510.5, 677.5)]


def run_lengths(mask, axis):
    m = mask if axis == 1 else mask.T
    out = np.zeros(m.shape, dtype=int)
    for r in range(m.shape[0]):
        row = m[r]
        c = 0
        while c < len(row):
            if row[c]:
                c2 = c
                while c2 < len(row) and row[c2]:
                    c2 += 1
                out[r, c:c2] = c2 - c
                c = c2
            else:
                c += 1
    return out if axis == 1 else out.T


solid = dark & ((run_lengths(dark, 1) >= 8) | (run_lengths(dark, 0) >= 8))


def clusters(idx, weights, gap=4):
    groups = []
    for i in idx:
        if groups and i - groups[-1][-1] <= gap:
            groups[-1].append(i)
        else:
            groups.append([i])
    return [float(np.average(g, weights=weights[g])) for g in groups]


out = []
for k in range(9):
    xl, xr = frame_x[k % 3]
    yt, yb = frame_y[k // 3]
    xl_i, xr_i = int(round(xl)) + 3, int(round(xr)) - 2
    yt_i, yb_i = int(round(yt)) + 3, int(round(yb)) - 2
    halo = binary_dilation(solid[yt_i:yb_i, xl_i:xr_i], iterations=3)
    gsub = gray[yt_i:yb_i, xl_i:xr_i] & ~halo
    h, w = gsub.shape

    cp = gsub.sum(axis=0).astype(float)
    rp = gsub.sum(axis=1).astype(float)
    gx = clusters(np.where(cp > 0.35 * h)[0], cp)
    gy = clusters(np.where(rp > 0.35 * w)[0], rp)
    assert len(gx) == 3 and len(gy) == 3, (k + 1, gx, gy)
    gx = [g + xl_i for g in gx]
    gy = [g + yt_i for g in gy]
    sx = (gx[2] - gx[0]) / 40.0
    sy = (gy[2] - gy[0]) / 40.0

    ssub = solid[yt_i:yb_i, xl_i:xr_i]
    cc = ssub.sum(axis=0).astype(float)
    rc = ssub.sum(axis=1).astype(float)
    vedges = clusters(np.where(cc > 0.5 * cc.max())[0], cc)
    hedges = clusters(np.where(rc > 0.5 * rc.max())[0], rc)
    xe0, xe1 = vedges[0] + xl_i, vedges[-1] + xl_i
    ye0, ye1 = hedges[0] + yt_i, hedges[-1] + yt_i

    L = (xe1 - xe0) / sx
    W = (ye1 - ye0) / sy
    out.append(dict(room=k + 1, L=round(L, 2), W=round(W, 2),
                    sx=round(sx, 4), sy=round(sy, 4),
                    x_m=[round((xe0 - gx[1]) / sx, 2), round((xe1 - gx[1]) / sx, 2)],
                    y_m=[round((gy[1] - ye1) / sy, 2), round((gy[1] - ye0) / sy, 2)],
                    n_vedges=len(vedges), n_hedges=len(hedges)))
    print(out[-1])

vols = [216, 224, 182, 3300, 5179, 3647, 8298, 8500, 7417]
print()
for r, V in zip(out, vols):
    H = V / (r['L'] * r['W'])
    print(f"room {r['room']}: L={r['L']:.2f} W={r['W']:.2f} H={H:.2f}")

with open(os.path.join(_here, 'digitized_rooms.json'), 'w') as f:
    json.dump(out, f, indent=1)
