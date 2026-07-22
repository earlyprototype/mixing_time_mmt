"""Reverberation time estimation from an RIR by Schroeder backward integration.

schroeder_rt60 fits a straight line to the Schroeder energy decay curve
between -5 dB and -25 dB (a T20 measurement) and extrapolates the slope to
60 dB of decay. Deterministic, no randomness involved.
"""

import numpy as np


def schroeder_edc_db(rir):
    """Schroeder energy decay curve in dB, normalized to 0 dB at t = 0."""
    e = np.asarray(rir, dtype=np.float64) ** 2
    edc = np.cumsum(e[::-1])[::-1]
    if edc[0] <= 0.0:
        raise ValueError("RIR has no energy")
    edc = edc / edc[0]
    with np.errstate(divide="ignore"):
        return 10.0 * np.log10(edc)


def schroeder_rt60(rir, fs, db_start=-5.0, db_end=-25.0):
    """RT60 estimate from backward integration and a T20-style linear fit.

    A least squares line is fitted to the decay curve over the region from
    db_start to db_end (defaults -5 to -25 dB) and the time to fall 60 dB
    is extrapolated from its slope.
    """
    db = schroeder_edc_db(rir)
    below_start = np.nonzero(db <= db_start)[0]
    below_end = np.nonzero(db <= db_end)[0]
    if below_start.size == 0 or below_end.size == 0:
        raise ValueError(
            f"decay never reaches the fit range [{db_start}, {db_end}] dB; "
            "RIR too short or too little decay"
        )
    i0, i1 = below_start[0], below_end[0]
    if i1 <= i0 + 1:
        raise ValueError("fit range too short for a line fit")
    t = np.arange(i0, i1 + 1, dtype=np.float64) / fs
    slope, _ = np.polyfit(t, db[i0:i1 + 1], 1)
    if slope >= 0.0:
        raise ValueError("non-decaying fit region")
    return -60.0 / slope
