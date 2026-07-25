"""Moving window HFD profile of an RIR, a port of the thesis appendix HMW.m.

MATLAB original (single channel, window size y, 1-based indexing):

    k = [x; zeros(y,1)];                 % pad y zeros at the end
    for i = 1:length(k)
        TraceSeg = k(i:(i+y));           % y+1 samples
        if length(TraceSeg) + i + 1 > length(k), break, end
        l(i) = hfd(TraceSeg);
    end
    sl = smooth(l, 200);

Length convention. With N = length(x), the padded signal has N + y samples
and TraceSeg has y + 1 samples, so the break condition (y+1) + i + 1 > N + y
first holds at i = N - 1. The break is checked BEFORE the call to hfd, so
the last computed profile value is at 1-based i = N - 2 and the profile has
exactly N - 2 samples. hfd_profile reproduces that: the returned array has
length len(rir) - 2 and its 0-based element s is the Higuchi FD of the
(window + 1)-sample segment padded[s : s + window + 1].

smooth_profile reproduces MATLAB smooth(l, 200). MATLAB reduces an even
span by 1 (so 200 becomes 199) and shrinks the window symmetrically at the
edges: at index i it averages over [i - half, i + half] with
half = min(i, n - 1 - i, (span - 1) / 2).
"""

import numpy as np

from .higuchi import higuchi_fd


def hfd_profile(rir, window=50, kmax=5):
    """Vectorized moving window Higuchi FD profile.

    Computes, for every window position at once, the Higuchi curve lengths
    per (k, m) from a single array of k-spaced absolute differences of the
    padded signal, then fits all the log-log slopes in one batched least
    squares step. Matches hfd_profile_reference to floating point rounding.

    Parameters
    ----------
    rir : array_like
        Room impulse response, 1-D.
    window : int
        MATLAB window parameter y. Segments have window + 1 samples.
    kmax : int
        Maximum delay passed to the Higuchi algorithm.

    Returns
    -------
    np.ndarray
        Profile of length len(rir) - 2 (see module docstring). Positions
        whose segment gives a zero curve length (constant segment) are nan.
    """
    rir = np.asarray(rir, dtype=float)
    n = rir.size
    p = n - 2  # profile length, see module docstring
    if p < 1:
        return np.empty(0)
    seglen = window + 1
    padded = np.concatenate([rir, np.zeros(window)])
    ks = np.arange(1, kmax + 1)
    lk = np.empty((kmax, p))
    for k in ks:
        dk = np.abs(padded[k:] - padded[:-k])
        acc = np.zeros(p)
        for m in range(1, k + 1):
            n_i = (seglen - m) // k
            if n_i < 1:
                lk[:] = np.nan
                return np.full(p, np.nan)
            ng = (seglen - 1) / (n_i * k)
            lmki = np.zeros(p)
            # Within a segment starting at s, term i of Lmki is
            # |padded[s + m-1 + i*k] - padded[s + m-1 + (i-1)*k]|, which is
            # dk[s + m-1 + (i-1)*k]. Accumulate over the n_i offsets for
            # all window positions s at once.
            for j in range(n_i):
                off = m - 1 + j * k
                lmki += dk[off : off + p]
            acc += (lmki * ng) / k
        lk[k - 1] = acc / k
    bad = np.any(lk == 0.0, axis=0)
    with np.errstate(divide="ignore"):
        y = np.log(lk)
    # Batched degree 1 polyfit of y against log(1/k): with the abscissa
    # centered, the slope is xc . y / (xc . xc).
    xk = np.log(1.0 / ks)
    xc = xk - xk.mean()
    slope = (xc @ y) / (xc @ xc)
    slope[bad] = np.nan
    return slope


def hfd_profile_reference(rir, window=50, kmax=5):
    """Slow reference implementation: a plain loop calling higuchi_fd.

    Mirrors the MATLAB HMW.m loop directly. Same output as hfd_profile.
    """
    rir = np.asarray(rir, dtype=float)
    n = rir.size
    p = n - 2
    if p < 1:
        return np.empty(0)
    padded = np.concatenate([rir, np.zeros(window)])
    out = np.empty(p)
    for s in range(p):
        out[s] = higuchi_fd(padded[s : s + window + 1], kmax)
    return out


def smooth_profile(profile, span=200):
    """Centered moving average matching MATLAB smooth(x, span).

    An even span is reduced by 1 (MATLAB behavior), so the default 200
    acts as 199. At the edges the window shrinks symmetrically: index i
    averages profile[i - half : i + half + 1] with
    half = min(i, n - 1 - i, (effective_span - 1) // 2).

    Non-finite profile values (hfd_profile emits nan for constant
    segments, e.g. the silence before the direct sound) are excluded from
    each window's average rather than propagated, so one degenerate window
    only affects outputs whose window overlaps it. An output is nan only
    when every value in its window is non-finite.
    """
    x = np.asarray(profile, dtype=float)
    n = x.size
    if n == 0:
        return x.copy()
    eff = span if span % 2 == 1 else span - 1
    hw = (eff - 1) // 2
    idx = np.arange(n)
    half = np.minimum(np.minimum(idx, n - 1 - idx), hw)
    good = np.isfinite(x)
    csum = np.concatenate([[0.0], np.cumsum(np.where(good, x, 0.0))])
    ccnt = np.concatenate([[0.0], np.cumsum(good.astype(float))])
    sums = csum[idx + half + 1] - csum[idx - half]
    cnts = ccnt[idx + half + 1] - ccnt[idx - half]
    out = np.full(n, np.nan)
    nz = cnts > 0
    out[nz] = sums[nz] / cnts[nz]
    return out
