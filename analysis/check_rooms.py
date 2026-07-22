"""Sanity check of the reconstructed Lindau rooms and one end-to-end RIR.

Prints the geometry table for all nine rooms, including the Sabine RT from
the reconstructed surfaces versus the tabulated RT (their ratio measures how
wrong the reconstructed proportions are, or how non-Sabine the real room
was), and the surface area that WOULD make Sabine match the tabulated RT.

Then generates the Room 1 RIR only (fast) with the pure image-source model,
verifies the energy decay is fully developed over the analysis window, and
reports the Schroeder RT60 against the tabulated 0.39 s.

Run from the repo root:  python analysis/check_rooms.py
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from mixtime.rooms import ROOMS, get_room
from mixtime.ism import generate_rir, required_max_order
from mixtime.rt_check import schroeder_rt60, schroeder_edc_db


def main():
    hdr = (f"{'rm':>2} {'name':<16} {'V':>6} {'alpha':>5} {'L':>6} {'W':>6} "
           f"{'H':>6} {'S':>7} {'RT_tab':>6} {'RT_sab':>6} {'ratio':>5} "
           f"{'S_match':>7} {'inc':>3}")
    print(hdr)
    print("-" * len(hdr))
    for r in ROOMS:
        ratio = r.sabine_rt_s / r.rt_s
        print(f"{r.number:>2} {r.name:<16} {r.volume_m3:>6.0f} "
              f"{r.alpha_ave:>5.2f} {r.L:>6.2f} {r.W:>6.2f} {r.H:>6.2f} "
              f"{r.S:>7.1f} {r.rt_s:>6.2f} {r.sabine_rt_s:>6.3f} "
              f"{ratio:>5.2f} {r.s_match_m2:>7.1f} "
              f"{'yes' if r.included else 'NO':>3}")

    print()
    room1 = get_room(1)
    order, needed = required_max_order(room1)
    print(f"Room 1 end-to-end check: max_order={order} "
          f"(uncapped requirement {needed})")
    t0 = time.time()
    res = generate_rir(room1)
    dt = time.time() - t0
    fs = res["fs"]
    n = res["length_samples"]
    print(f"  generated in {dt:.1f} s, length={n} samples "
          f"({n / fs:.4f} s), raw untrimmed length={res['raw_length']}")
    print(f"  source={res['source']}, mic={res['mic']}, "
          f"positions_clipped={res['positions_clipped']}, "
          f"order_capped={res['order_capped']}")

    # Energy decay development over the window of interest. Two checks:
    # 1. the raw RIR extends well past the analysis window, so no image
    #    source truncation artifact can fall inside the window;
    # 2. the raw Schroeder level at the window end, for information. Since
    #    the window length is one tabulated RT, the expected level there is
    #    about -60 dB times (tabulated RT / actual ISM RT), not -60 dB.
    db_raw = schroeder_edc_db(res["raw_rir"])
    level_at_window_end = db_raw[min(n, db_raw.size) - 1]
    coverage = res["raw_length"] / n
    print(f"  raw RIR covers {coverage:.2f}x the analysis window "
          f"(want > 1, comfortably)")
    print(f"  raw Schroeder level at window end: "
          f"{level_at_window_end:.1f} dB")

    rt60 = schroeder_rt60(res["rir"], fs)
    print(f"  Schroeder RT60 (T20 fit, -5 to -25 dB): {rt60:.3f} s")
    print(f"  tabulated RT: {room1.rt_s:.2f} s, "
          f"Sabine from geometry: {room1.sabine_rt_s:.3f} s, "
          f"measured/tabulated = {rt60 / room1.rt_s:.3f}")


if __name__ == "__main__":
    main()
