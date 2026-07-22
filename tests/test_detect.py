"""Tests for tail statistics and threshold crossing detection.

Indexing reminder: detect_crossing takes and returns MATLAB 1-based sample
numbers (the thesis convention). 0-based array index = returned value - 1.
"""

import numpy as np

from mixtime import tail_stats, detect_crossing, CRITERIA


def test_tail_stats_hand_checked():
    # n=20, tail_frac=0.1: xtenpercent = round(2.0) = 2, so the tail is
    # the last 3 samples (MATLAB slice is xtenpercent + 1 samples).
    profile = np.arange(20, dtype=float)
    mean, sd = tail_stats(profile, tail_frac=0.1)
    assert mean == 18.0  # mean of [17, 18, 19]
    assert abs(sd - 1.0) < 1e-12  # std with ddof=1 of [17, 18, 19]


def test_detect_crossing_first_index():
    # 50 zeros then 50 twos. Tail = last 11 samples, all 2.0, so the
    # threshold is 2.0 for every k. First crossing is 0-based index 50,
    # i.e. 1-based sample 51.
    profile = np.concatenate([np.zeros(50), np.full(50, 2.0)])
    assert detect_crossing(profile, k=0, start=1) == 51
    assert detect_crossing(profile, k=2, start=1) == 51


def test_start_parameter_respected():
    # The profile crosses at sample 51, but scanning must begin at start.
    profile = np.concatenate([np.zeros(50), np.full(50, 2.0)])
    assert detect_crossing(profile, k=0, start=60) == 60
    assert detect_crossing(profile, k=0, start=51) == 51


def _ramp_plateau_profile():
    # Deterministic synthetic profile: linear ramp from 0 to 2 over 500
    # samples, then a noisy plateau around 2 for 500 samples. The tail
    # (last 101 samples) sits inside the plateau.
    rng = np.random.default_rng(42)
    ramp = np.linspace(0.0, 2.0, 500)
    plateau = 2.0 + 0.05 * rng.standard_normal(500)
    return np.concatenate([ramp, plateau])


def test_criteria_ordering_monotonic():
    # A larger k lowers the threshold, so on a rising profile the crossing
    # can only move earlier (or stay put). Criterion I..IV is k=0..3.
    profile = _ramp_plateau_profile()
    crossings = [
        detect_crossing(profile, k=CRITERIA[c], start=1)
        for c in ("I", "II", "III", "IV")
    ]
    assert all(c is not None for c in crossings)
    for earlier_k, later_k in zip(crossings[1:], crossings[:-1]):
        assert earlier_k <= later_k


def test_none_when_threshold_never_reached():
    # A strongly negative k pushes the threshold above every profile
    # value, so no sample ever crosses it.
    profile = _ramp_plateau_profile()
    assert detect_crossing(profile, k=-100, start=1) is None


def test_detect_determinism():
    profile = _ramp_plateau_profile()
    assert detect_crossing(profile, k=2, start=1) == detect_crossing(
        profile.copy(), k=2, start=1
    )
