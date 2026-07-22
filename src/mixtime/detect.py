"""Threshold detection on the HFD profile, a port of the thesis FDthres.m.

MATLAB original:

    xtenpercent = round(length(x)/10);
    tail = x((length(x)-xtenpercent):length(x));   % xtenpercent+1 samples
    thresholds: mean(tail) - n*std(tail) for Criterion n+1, n = 0..3
    for i = 1000:length(x)
        if x(i) >= threshold, break, end

Indexing note. MATLAB is 1-based: the loop starts at 1-based sample 1000,
which is 0-based index 999. detect_crossing takes `start` in the MATLAB
1-based convention (default 1000) and also RETURNS the MATLAB 1-based
sample number, because the thesis tables are 1-based. Subtract 1 to get a
0-based index into the profile array.
"""

import numpy as np

# Criterion number -> k in threshold = tail_mean - k * tail_sd.
CRITERIA = {"I": 0, "II": 1, "III": 2, "IV": 3}


def _matlab_round(v):
    """MATLAB round: half away from zero (v is nonnegative here)."""
    return int(np.floor(v + 0.5))


def tail_stats(profile, tail_frac=0.1):
    """Mean and standard deviation of the profile tail.

    The tail is the last round(n * tail_frac) + 1 samples, matching the
    MATLAB slice x((end - xtenpercent):end). The standard deviation uses
    ddof=1 to match MATLAB std.

    Returns
    -------
    (mean, sd) : tuple of floats
    """
    x = np.asarray(profile, dtype=float)
    n = x.size
    m = _matlab_round(n * tail_frac)
    tail = x[n - m - 1 :]
    return float(tail.mean()), float(tail.std(ddof=1))


def detect_crossing(profile, k, tail_frac=0.1, start=1000):
    """First sample at which the profile reaches the tail threshold.

    threshold = tail_mean - k * tail_sd. Criterion I is k=0, II is k=1,
    III is k=2, IV is k=3 (see CRITERIA).

    Parameters
    ----------
    profile : array_like
        Smoothed HFD profile.
    k : float
        Number of tail standard deviations subtracted from the tail mean.
    tail_frac : float
        Fraction of the profile used for the tail statistics.
    start : int
        First sample to consider, in the MATLAB 1-based convention. The
        thesis uses 1000, i.e. 0-based index 999.

    Returns
    -------
    int or None
        The MATLAB 1-based sample number of the first sample at or after
        `start` with profile >= threshold, or None if never crossed.
        (The MATLAB loop would silently return length(x) in that case,
        this port returns None instead.)
    """
    x = np.asarray(profile, dtype=float)
    mean, sd = tail_stats(x, tail_frac)
    threshold = mean - k * sd
    i0 = start - 1  # 0-based index of the 1-based start sample
    hits = np.nonzero(x[i0:] >= threshold)[0]
    if hits.size == 0:
        return None
    return int(i0 + hits[0] + 1)  # back to 1-based
