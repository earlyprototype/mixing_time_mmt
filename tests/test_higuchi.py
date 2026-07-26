"""Validation of higuchi_fd against signals of known fractal dimension.

All random signals are generated with seeded numpy Generators, so every
test is deterministic. Tolerances are stated per test.
"""

import numpy as np
import pytest

from mixtime import higuchi_fd

SEED = 20260722


def fgn_davies_harte(n, hurst, rng):
    """Fractional Gaussian noise by circulant spectral synthesis.

    Davies-Harte method: embed the exact fGn autocovariance
    gamma(k) = 0.5*(|k+1|^2H - 2|k|^2H + |k-1|^2H) in a circulant matrix
    of size 2n - 2, whose eigenvalues are the FFT of its first row, and
    color complex Gaussian noise with the eigenvalue square roots. Tiny
    negative eigenvalues from rounding are clipped to zero.
    """
    k = np.arange(n)
    g = 0.5 * (
        np.abs(k + 1) ** (2 * hurst)
        - 2 * np.abs(k) ** (2 * hurst)
        + np.abs(k - 1) ** (2 * hurst)
    )
    row = np.concatenate([g, g[-2:0:-1]])
    lam = np.maximum(np.fft.fft(row).real, 0.0)
    m = row.size
    w = np.zeros(m, dtype=complex)
    w[0] = rng.standard_normal()
    w[m // 2] = rng.standard_normal()
    re = rng.standard_normal(m // 2 - 1)
    im = rng.standard_normal(m // 2 - 1)
    w[1 : m // 2] = (re + 1j * im) / np.sqrt(2.0)
    w[m // 2 + 1 :] = np.conj(w[1 : m // 2][::-1])
    y = np.fft.ifft(np.sqrt(lam) * w) * np.sqrt(m)
    return y.real[:n]


def test_white_noise_full_length():
    # White Gaussian noise has FD 2. Tolerance per task spec: 1.8 to 2.05
    # for the full length estimate (N=5000, kmax=5).
    rng = np.random.default_rng(SEED)
    noise = rng.standard_normal(5000)
    fd = higuchi_fd(noise, kmax=5)
    assert 1.8 <= fd <= 2.05


def test_white_noise_production_windows():
    # Same noise chopped into the 51-sample windows used in production
    # (window=50, kmax=5). Single short windows are noisy estimators, so
    # each window is allowed 1.7 to 2.2 and the mean must be in 1.9 to 2.1
    # (measured: mean 2.003, min 1.851, max 2.135 for this seed).
    rng = np.random.default_rng(SEED)
    noise = rng.standard_normal(5000)
    fds = np.array(
        [higuchi_fd(noise[s : s + 51], kmax=5) for s in range(0, 5000 - 51, 51)]
    )
    assert np.all(fds >= 1.7) and np.all(fds <= 2.2)
    assert 1.9 <= fds.mean() <= 2.1


def test_sine_wave():
    # A densely sampled smooth curve has FD 1. Five cycles over 5000
    # samples, i.e. 1000 samples per period. Tolerance 0.95 to 1.1.
    t = np.linspace(0.0, 1.0, 5000)
    fd = higuchi_fd(np.sin(2 * np.pi * 5 * t), kmax=5)
    assert 0.95 <= fd <= 1.1


@pytest.mark.parametrize("hurst", [0.2, 0.5, 0.8])
def test_fbm_known_dimension(hurst):
    # Fractional Brownian motion (cumsum of exact fGn) has graph FD
    # 2 - H. N=10000. kmax=10 is used because fBm needs a longer range of
    # delays than the production kmax=5 for a stable slope (the task spec
    # allows kmax >= 8). Tolerance 0.15 per spec. Measured for seed 1234:
    # H=0.2 -> 1.801, H=0.5 -> 1.497, H=0.8 -> 1.178.
    rng = np.random.default_rng(1234)
    fbm = np.cumsum(fgn_davies_harte(10000, hurst, rng))
    fd = higuchi_fd(fbm, kmax=10)
    assert abs(fd - (2.0 - hurst)) <= 0.15


def test_weierstrass_known_dimension():
    # W(t) = sum_j b^(-a*j) cos(b^j * pi * t) has box dimension 2 - a
    # (holds for the cosine Weierstrass function whenever b > 1 and
    # b^(1-a) > 1). Here a=0.5, b=3, t in [0,1] with 10000 samples, terms
    # j=0..8 (3^8 = 6561 stays below the 5000 cycles per unit Nyquist
    # limit of this sampling). kmax=10, tolerance 0.15 per spec.
    # Measured: 1.489 vs expected 1.5.
    a, b = 0.5, 3.0
    t = np.linspace(0.0, 1.0, 10000)
    w = np.zeros_like(t)
    for j in range(9):
        w += b ** (-a * j) * np.cos(b**j * np.pi * t)
    fd = higuchi_fd(w, kmax=10)
    assert abs(fd - (2.0 - a)) <= 0.15


def test_constant_signal_returns_nan():
    assert np.isnan(higuchi_fd(np.ones(100), kmax=5))


def test_determinism():
    rng = np.random.default_rng(SEED)
    x = rng.standard_normal(2000)
    assert higuchi_fd(x, kmax=5) == higuchi_fd(x.copy(), kmax=5)
