"""Higuchi fractal dimension, a port of the thesis appendix hfd.m.

MATLAB original (1-based indexing):

    for k=1:kmax
        for m=1:k
            Lmki = sum over i=1..fix((N-m)/k) of abs(x(m+i*k) - x(m+(i-1)*k))
            Ng = (N-1) / (fix((N-m)/k)*k)
            Lmk(m,k) = (Lmki*Ng)/k
    Lk(k) = sum(Lmk(1:k,k))/k
    b = polyfit(log(1./(1:kmax)), log(Lk), 1); FD = b(1)

In 0-based Python the samples x(m + i*k) for i = 0..fix((N-m)/k) are exactly
the strided slice x[m-1::k], which has fix((N-m)/k) + 1 elements, so Lmki is
the sum of absolute first differences of that slice.
"""

import numpy as np


def higuchi_fd(x, kmax=5):
    """Estimate the Higuchi fractal dimension of a 1-D signal.

    Parameters
    ----------
    x : array_like
        Time series (the thesis uses 51-sample RIR segments).
    kmax : int
        Maximum delay k. The thesis default is 5.

    Returns
    -------
    float
        The FD estimate, about 1 for smooth curves and up to 2 for noise.
        Returns nan for degenerate input, i.e. a constant signal (some
        curve length Lk is 0 so its log is undefined) or a signal too
        short to form a single k-spaced difference.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    ks = np.arange(1, kmax + 1)
    lk = np.empty(kmax)
    for k in ks:
        lmk_sum = 0.0
        for m in range(1, k + 1):
            n_i = (n - m) // k
            if n_i < 1:
                return float("nan")
            # Sum over i of abs(x[m-1 + i*k] - x[m-1 + (i-1)*k]), i=1..n_i.
            lmki = np.abs(np.diff(x[m - 1 :: k])).sum()
            ng = (n - 1) / (n_i * k)
            lmk_sum += (lmki * ng) / k
        lk[k - 1] = lmk_sum / k
    if np.any(lk == 0.0):
        return float("nan")
    coeffs = np.polyfit(np.log(1.0 / ks), np.log(lk), 1)
    return float(coeffs[0])
