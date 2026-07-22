"""Thesis-faithful Allen and Berkley (1979) image source model.

Port of the MATLAB sroom.m / lthimage.m reproduced in the thesis appendix
(9.1.4), which is itself a port of the Fortran code in Allen and Berkley's
original paper. Distinctive properties versus a modern ISM such as
pyroomacoustics:

- All geometry is converted to units of samples (metres * fs / c).
- Each image's delay is rounded to an integer sample, so the RIR is a
  sparse train of single-sample spikes (no fractional-delay sinc kernels).
- Amplitude is prod(beta terms) / (distance in samples), with
  beta = sqrt(1 - alpha) the pressure reflection coefficient of each wall.
- The result is highpass filtered with butter(3, 0.01) as in sroom.m
  (cutoff 0.01 of Nyquist, about 220 Hz at fs = 44100).

The thesis fed the same beta to all six walls. Deterministic, no
randomness anywhere.
"""

import numpy as np
from scipy.signal import butter, lfilter

C_SOUND = 343.0


def ab_rir(room, fs=44100, length_samples=None, c=C_SOUND):
    """Generate the Allen-Berkley RIR for a reconstructed room.

    Parameters mirror mixtime.ism.generate_rir. Source (L/2, 2.0, 1.9),
    receiver (L/2, 1.0, 1.9), per the thesis text.

    Returns dict with rir (float64, length_samples), n_images used, and
    the pre-filter sparse RIR for inspection.
    """
    if length_samples is None:
        length_samples = int(round(room.rt_s * fs))
    npts = int(length_samples)
    scale = fs / c
    rl = np.array([room.L, room.W, room.H]) * scale
    r = np.array([room.L / 2.0, 1.0, 1.9]) * scale   # receiver
    r0 = np.array([room.L / 2.0, 2.0, 1.9]) * scale  # source
    beta = np.sqrt(1.0 - room.alpha_ave)

    # Range of image orders, as in sroom.m: n = floor(npts/(2*rl)) + 1.
    n1 = int(np.floor(npts / (rl[0] * 2) + 1))
    n2 = int(np.floor(npts / (rl[1] * 2) + 1))
    n3 = int(np.floor(npts / (rl[2] * 2) + 1))

    nx = np.arange(-n1, n1 + 1)
    ny = np.arange(-n2, n2 + 1)
    nz = np.arange(-n3, n3 + 1)
    NX, NY, NZ = np.meshgrid(nx, ny, nz, indexing="ij")
    NX = NX.ravel()[:, None]  # (M, 1)
    NY = NY.ravel()[:, None]
    NZ = NZ.ravel()[:, None]

    # Eight sign permutations p in {0,1}^3 (l, j, k in the MATLAB code).
    P = np.array([[l, j, k] for l in (0, 1) for j in (0, 1) for k in (0, 1)])
    sign = 1 - 2 * P  # +1 for p=0, -1 for p=1

    # Image position relative to receiver: (r +/- r0) - 2 n rl, per axis.
    dx = (r[0] + sign[None, :, 0] * r0[0]) - 2.0 * NX * rl[0]
    dy = (r[1] + sign[None, :, 1] * r0[1]) - 2.0 * NY * rl[1]
    dz = (r[2] + sign[None, :, 2] * r0[2]) - 2.0 * NZ * rl[2]
    dist = np.sqrt(dx * dx + dy * dy + dz * dz)
    delay = np.round(dist).astype(np.int64)  # integer sample delay

    # Wall bounce counts, as in sroom.m: beta1^|nx-l| beta2^|nx|
    # beta3^|ny-j| beta4^|ny| beta5^|nz-k| beta6^|nz|, all betas equal.
    expo = (np.abs(NX - P[None, :, 0]) + np.abs(NX)
            + np.abs(NY - P[None, :, 1]) + np.abs(NY)
            + np.abs(NZ - P[None, :, 2]) + np.abs(NZ))
    with np.errstate(divide="ignore", invalid="ignore"):
        amp = np.power(beta, expo) / delay

    keep = (delay >= 1) & (delay < npts)
    ht = np.zeros(npts, dtype=np.float64)
    np.add.at(ht, delay[keep], amp[keep])
    n_images = int(keep.sum())

    b, a = butter(3, 0.01, "high")
    rir = lfilter(b, a, ht)
    return {
        "rir": rir,
        "sparse_rir": ht,
        "n_images": n_images,
        "fs": fs,
        "length_samples": npts,
        "beta": float(beta),
    }
