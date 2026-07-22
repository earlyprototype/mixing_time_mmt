"""Higuchi FD mixing time pipeline, a Python port of the 2011 thesis MATLAB code.

Modules:
    higuchi: Higuchi fractal dimension of a single segment (hfd.m).
    profile: moving window HFD profile and smoothing (HMW.m, smooth).
    detect: tail statistics and threshold crossing detection (FDthres.m).
"""

from .higuchi import higuchi_fd
from .profile import hfd_profile, hfd_profile_reference, smooth_profile
from .detect import tail_stats, detect_crossing, CRITERIA

__all__ = [
    "higuchi_fd",
    "hfd_profile",
    "hfd_profile_reference",
    "smooth_profile",
    "tail_stats",
    "detect_crossing",
    "CRITERIA",
]
