"""Room definitions for the nine Lindau et al. 2010 rooms (thesis Table 4.3.1.1).

Tabulated data (volume, average absorption coefficient, reverberation time)
are transcribed exactly from Table 4.3.1.1 of the 2011 thesis, which in turn
transcribed Lindau et al. 2010. Only the volume is tabulated, so footprint
length L and width W were digitized from the floor plans of Figure 4.3.1.1
(thesis page 42). Height is derived as H = V / (L * W).

Digitization method (see inputs/rooms_reconstruction.md for full detail):
the native 723x755 px raster of Figure 4.3.1.1 was extracted from
inputs/FINAL_thom_conaty.pdf (page 42, image xref 256) with pymupdf. Each of
the nine subplots was calibrated on its dotted gridlines, which sit at data
coordinates -20, 0 and +20 m on both axes (x scale 3.95 px/m, y scale about
3.49 to 3.68 px/m depending on subplot row). Room outlines were located as
solid black strokes and edge positions taken as stroke-centre centroids.
Estimated digitization uncertainty is about +/- 1 px, i.e. +/- 0.3 m per edge.

Sanity rule applied: 2 m < H < 25 m. All nine derived heights fall inside
this band, so no proportional adjustment of L and W was needed for any room.

Room 9 (JC church) is non-rectangular in the figure (notched outline). Its
L and W here describe the bounding rectangle of the outline, so H derived
from V / (L * W) underestimates the true height. The thesis excluded Room 9
from modelling as non-shoebox; it is kept here with included=False.

Everything in this module is deterministic and derived at import time.
"""

from dataclasses import dataclass


SABINE_COEFF = 0.161  # s/m, as used in the thesis (Sabine, air at ~20 C)


@dataclass(frozen=True)
class Room:
    number: int
    name: str
    volume_m3: float
    alpha_ave: float
    rt_s: float          # tabulated RT from Table 4.3.1.1
    L: float             # footprint length in m, digitized from Fig. 4.3.1.1
    W: float             # footprint width in m, digitized from Fig. 4.3.1.1
    H: float             # derived, V / (L * W)
    S: float             # total shoebox surface area, 2(LW + LH + WH)
    sabine_rt_s: float   # 0.161 V / (S alpha), from reconstructed geometry
    s_match_m2: float    # surface area that WOULD make Sabine equal rt_s
    included: bool       # False for Room 9, excluded by the thesis
    notes: str


def _make(number, name, volume, alpha, rt, L, W, included, notes):
    H = volume / (L * W)
    if not (2.0 < H < 25.0):
        # By construction this does not trigger for the digitized values.
        # Kept as a hard guard so any future edit of L or W cannot silently
        # produce an absurd height.
        raise ValueError(
            f"room {number}: derived H={H:.2f} m outside (2, 25); "
            "adjust L and W proportionally and record the change"
        )
    S = 2.0 * (L * W + L * H + W * H)
    sabine = SABINE_COEFF * volume / (S * alpha)
    s_match = SABINE_COEFF * volume / (rt * alpha)
    return Room(
        number=number,
        name=name,
        volume_m3=float(volume),
        alpha_ave=float(alpha),
        rt_s=float(rt),
        L=float(L),
        W=float(W),
        H=round(H, 3),
        S=round(S, 2),
        sabine_rt_s=round(sabine, 4),
        s_match_m2=round(s_match, 2),
        included=included,
        notes=notes,
    )


_FIG = "L, W digitized from Fig. 4.3.1.1"

ROOMS = [
    _make(1, "E. Studio", 216, 0.36, 0.39, 9.48, 7.62, True,
          _FIG + "; H derived from V"),
    _make(2, "EN-111", 224, 0.26, 0.62, 6.95, 7.18, True,
          _FIG + "; H derived from V"),
    _make(3, "EN-190", 182, 0.17, 0.79, 8.60, 7.18, True,
          _FIG + "; H derived from V"),
    _make(4, "H-0104", 3300, 0.28, 1.15, 25.69, 20.71, True,
          _FIG + "; H derived from V"),
    _make(5, "HE-101", 5179, 0.23, 1.67, 26.83, 25.86, True,
          _FIG + "; H derived from V"),
    _make(6, "Teldex", 3647, 0.20, 1.83, 28.35, 16.14, True,
          _FIG + "; H derived from V"),
    _make(7, "UoA concert hall", 8298, 0.33, 1.52, 37.72, 16.46, True,
          _FIG + "; H derived from V"),
    _make(8, "TUB Audimax", 8500, 0.23, 2.08, 29.74, 27.89, True,
          _FIG + "; H derived from V"),
    _make(9, "JC church", 7417, 0.23, 2.36, 25.82, 14.28, False,
          "EXCLUDED by thesis as non-shoebox; L, W are the bounding "
          "rectangle of the notched outline in Fig. 4.3.1.1, whose area "
          "exceeds the true footprint, so the derived H underestimates "
          "the true height"),
]


def get_room(number):
    """Return the Room with the given 1-based number."""
    for r in ROOMS:
        if r.number == number:
            return r
    raise KeyError(f"no room number {number}")
