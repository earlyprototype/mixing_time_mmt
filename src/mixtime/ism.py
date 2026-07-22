"""Image-source RIR generation for the reconstructed Lindau rooms.

Uses pyroomacoustics ShoeBox in pure image-source mode: no air absorption,
no ray tracing, no randomized image method. The whole computation is
deterministic, there is nothing random to seed.

Absorption convention. pyroomacoustics takes the ENERGY absorption
coefficient (pra.Material(energy_absorption=alpha)). The thesis instead
converted alpha_ave to a pressure reflection coefficient
beta = sqrt(1 - alpha) and fed beta to an Allen and Berkley (1979) code.
Energy reflection there is beta**2 = 1 - alpha, so feeding alpha_ave to
pyroomacoustics is equivalent in energy terms.

Geometry convention: x = length L, y = width W, z = height H.
Thesis text (Section 5.1): receiver at (L/2, 1, 1.9), source at
(L/2, 2, 1.9). The thesis GUI screenshot (Figure 5.1a) swapped source and
receiver y values, but the ISM transfer path is reciprocal so the RIR is
identical either way; we follow the text. Positions are clipped to stay at
least 0.5 m from every wall (matters only if a future room is very small;
none of the nine reconstructed rooms triggers the clip).

max_order. Chosen so image sources cover the requested RIR length:
ceil(c * t_len / min_dim) + 3 safety, capped at MAX_ORDER_CAP = 100 to keep
runtime and memory tractable for the big rooms. With the reconstructed
geometries and default lengths (round(rt_s * fs) samples) the per-room
required orders are 42 to 95, so the cap never truncates a default request.
If the cap ever binds, the returned metadata flags it and the energy decay
over the window of interest should be re-checked (the Schroeder curve of
the raw RIR must be fully developed, about -60 dB, at t_len).

Default RIR length: length_samples = round(rt_s * fs). Assumption, recorded
here and in inputs/rooms_reconstruction.md. For Room 1 at fs = 44100 this
gives 17199 samples, matching the roughly 17200-sample profile length seen
in thesis Figure 5.3.
"""

import math

import numpy as np
import pyroomacoustics as pra


SPEED_OF_SOUND = 343.0   # m/s, pyroomacoustics default
MAX_ORDER_CAP = 100      # documented cap, see module docstring
POSITION_MARGIN = 0.5    # m, minimum distance of source and mic from walls
SAFETY_ORDERS = 3


def default_length_samples(room, fs=44100):
    """Assumed RIR length: tabulated RT times sample rate, rounded."""
    return int(round(room.rt_s * fs))


def required_max_order(room, fs=44100, length_samples=None,
                       c=SPEED_OF_SOUND):
    """Reflection order needed for image sources to cover the RIR window.

    The most distant image source of order n is at least n * min_dim away,
    so covering time t_len needs roughly c * t_len / min_dim orders. A small
    safety margin is added; the result is capped at MAX_ORDER_CAP.
    """
    if length_samples is None:
        length_samples = default_length_samples(room, fs)
    t_len = length_samples / fs
    min_dim = min(room.L, room.W, room.H)
    order = math.ceil(c * t_len / min_dim) + SAFETY_ORDERS
    return min(order, MAX_ORDER_CAP), order


def _clip_position(pos, dims, margin=POSITION_MARGIN):
    """Clip a position to lie at least margin inside the box; report if so."""
    clipped = []
    out = []
    for p, d in zip(pos, dims):
        lo, hi = margin, d - margin
        if lo > hi:
            raise ValueError(f"room dimension {d} m too small for "
                             f"{margin} m margins")
        q = min(max(p, lo), hi)
        clipped.append(q != p)
        out.append(q)
    return out, any(clipped)


def generate_rir(room, fs=44100, length_samples=None, max_order=None):
    """Generate a pure ISM RIR for one reconstructed room.

    Parameters
    ----------
    room : mixtime.rooms.Room
    fs : int, sample rate (thesis used 44100).
    length_samples : int or None. None means round(room.rt_s * fs).
    max_order : int or None. None means required_max_order(...) (capped).

    Returns
    -------
    dict with keys:
        rir : float64 array of exactly length_samples (trimmed or
            zero-padded). The pyroomacoustics global fractional-delay
            offset (frac_delay_length // 2 samples) is removed so the
            direct sound arrives at distance / c, as in Allen and Berkley.
        raw_length : length in samples of the untrimmed pyroomacoustics RIR
            (after removing the global delay).
        raw_rir : the untrimmed float64 RIR (for decay checks).
        fs, length_samples, max_order : as used.
        order_capped : True if the required order exceeded MAX_ORDER_CAP.
        source, mic : positions actually used, in m.
        positions_clipped : True if either position had to be moved to
            respect the 0.5 m wall margin.
    """
    if length_samples is None:
        length_samples = default_length_samples(room, fs)
    capped_order, needed_order = required_max_order(room, fs, length_samples)
    if max_order is None:
        max_order = capped_order
    order_capped = needed_order > max_order

    dims = [room.L, room.W, room.H]
    source, s_clip = _clip_position([room.L / 2.0, 2.0, 1.9], dims)
    mic, m_clip = _clip_position([room.L / 2.0, 1.0, 1.9], dims)

    shoebox = pra.ShoeBox(
        dims,
        fs=fs,
        materials=pra.Material(energy_absorption=room.alpha_ave),
        max_order=max_order,
        air_absorption=False,
        ray_tracing=False,
        use_rand_ism=False,
    )
    shoebox.add_source(source)
    shoebox.add_microphone(mic)
    shoebox.compute_rir()
    raw = np.asarray(shoebox.rir[0][0], dtype=np.float64)

    # Remove the causal global delay pyroomacoustics prepends for its
    # fractional-delay filters, so t = 0 is the source emission time.
    global_delay = pra.constants.get("frac_delay_length") // 2
    raw = raw[global_delay:]

    rir = np.zeros(length_samples, dtype=np.float64)
    n = min(length_samples, raw.size)
    rir[:n] = raw[:n]

    return {
        "rir": rir,
        "raw_length": int(raw.size),
        "raw_rir": raw,
        "fs": fs,
        "length_samples": int(length_samples),
        "max_order": int(max_order),
        "order_capped": bool(order_capped),
        "source": source,
        "mic": mic,
        "positions_clipped": bool(s_clip or m_clip),
    }
