"""Higuchi FD mixing time pipeline, a Python port of the 2011 thesis MATLAB code.

Modules:
    higuchi: Higuchi fractal dimension of a single segment (hfd.m).
    profile: moving window HFD profile and smoothing (HMW.m, smooth).
    detect: tail statistics and threshold crossing detection (FDthres.m).
    rooms: the nine reconstructed Lindau et al. 2010 rooms (Table 4.3.1.1).
    rt_check: Schroeder backward-integration RT60 estimation.
    ism: image-source RIR generation (import mixtime.ism explicitly; it is
        not imported here to keep the pyroomacoustics dependency optional).
"""

from .higuchi import higuchi_fd
from .profile import hfd_profile, hfd_profile_reference, smooth_profile
from .detect import tail_stats, detect_crossing, CRITERIA
from .rooms import Room, ROOMS, get_room
from .rt_check import schroeder_rt60, schroeder_edc_db

__all__ = [
    "Room",
    "ROOMS",
    "get_room",
    "schroeder_rt60",
    "schroeder_edc_db",
    "higuchi_fd",
    "hfd_profile",
    "hfd_profile_reference",
    "smooth_profile",
    "tail_stats",
    "detect_crossing",
    "CRITERIA",
]
