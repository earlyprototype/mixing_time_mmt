"""Tests for the moving window HFD profile and MATLAB-style smoothing."""

import numpy as np

from mixtime import hfd_profile, hfd_profile_reference, smooth_profile


def test_fast_equals_reference():
    # The vectorized profile must match the plain per-window loop to 1e-9
    # on a short random signal.
    rng = np.random.default_rng(7)
    sig = rng.standard_normal(300)
    fast = hfd_profile(sig, window=50, kmax=5)
    ref = hfd_profile_reference(sig, window=50, kmax=5)
    assert fast.shape == ref.shape
    np.testing.assert_allclose(fast, ref, atol=1e-9, rtol=0)


def test_profile_length_convention():
    # HMW.m checks its break condition before calling hfd, so the last
    # computed 1-based index is N - 2 and the profile has N - 2 samples.
    rng = np.random.default_rng(3)
    sig = rng.standard_normal(500)
    assert hfd_profile(sig).size == 498
    assert hfd_profile_reference(sig).size == 498


def test_noise_profile_settles_near_two():
    # The raw profile of pure white noise should hover around FD 2 away
    # from the zero-padded end. Measured for this seed: mid-region mean
    # 1.987, range 1.77 to 2.15.
    rng = np.random.default_rng(99)
    noise = rng.standard_normal(5000)
    prof = hfd_profile(noise, window=50, kmax=5)
    mid = prof[1000:3000]
    assert 1.9 <= mid.mean() <= 2.1
    assert np.all(mid >= 1.6) and np.all(mid <= 2.3)


def test_smoothing_preserves_length():
    rng = np.random.default_rng(11)
    prof = rng.standard_normal(1234)
    assert smooth_profile(prof, span=200).size == prof.size


def test_smoothing_hand_computed_small_case():
    # span=5 (odd, half width 2) with symmetric edge shrinking, computed
    # by hand.
    x = np.array([1.0, 4.0, 2.0, 8.0, 5.0, 7.0, 3.0])
    expected = np.array(
        [
            1.0,  # half 0
            (1 + 4 + 2) / 3.0,  # half 1
            (1 + 4 + 2 + 8 + 5) / 5.0,  # half 2
            (4 + 2 + 8 + 5 + 7) / 5.0,  # half 2
            (2 + 8 + 5 + 7 + 3) / 5.0,  # half 2
            (5 + 7 + 3) / 3.0,  # half 1
            3.0,  # half 0
        ]
    )
    np.testing.assert_allclose(smooth_profile(x, span=5), expected, atol=1e-12)


def test_even_span_reduced_by_one():
    # MATLAB smooth reduces an even span by 1, so span=6 must equal span=5.
    rng = np.random.default_rng(21)
    x = rng.standard_normal(50)
    np.testing.assert_allclose(
        smooth_profile(x, span=6), smooth_profile(x, span=5), atol=1e-12
    )


def test_profile_determinism():
    rng = np.random.default_rng(5)
    sig = rng.standard_normal(400)
    a = hfd_profile(sig)
    b = hfd_profile(sig.copy())
    np.testing.assert_array_equal(a, b)
