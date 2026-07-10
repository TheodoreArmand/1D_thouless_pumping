#!/usr/bin/env python3
"""Parking-acceptance check for the R3a occupancy gate.

The R3a design (vs3_n2_k32_success_roadmap_report.html section 6-R3 / 8.2)
requires that low-amplitude reserve pads stay *parked* until they are populated:
while a term's amplitude |u_k| is below the occupancy threshold u_on, its centre
R_k and its width must not drift. This script verifies that hard acceptance gate
offline from a run's snapshots.csv (N=2 schema).

For every basis term k it tracks over the recorded events:
  - |u_k(t)|
  - R_k(t)                         (2-vector, complex; parking uses the real part)
  - lambda_min[Re(A_k + B_k)]      (square-integrability width floor)
A term is a *reserve* if |u_k(0)| < u_on. For each reserve, over the prefix of
events where it is still parked (|u_k(t)| < u_on) we require:
  - ||Re R_k(t) - Re R_k(0)|| < r_tol      (default 1.0 code units ~ a quarter
                                            short-lattice period)
  - |width(t) - width(0)| / max(width(0), eps) < w_tol   (default 0.25)

Prints a per-term table and an overall PASS/FAIL. Exit code 0 on PASS, 1 on FAIL.

Usage:
  check_occupancy_gate_parking.py <run_dir_or_snapshots.csv> [--u-on 0.03]
                                  [--r-tol 1.0] [--w-tol 0.25]
"""
import argparse
import csv
import math
import os
import sys


def load_snapshots(path):
    """Return (events, ncols) where events is a list of dicts keyed by
    basis_index -> row dict, ordered by first appearance (step order)."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"empty snapshots file: {path}")
    # Group rows by (snapshot, step) event, preserving order.
    events = []
    cur_key = None
    cur = None
    for r in rows:
        key = (r.get("snapshot"), r.get("step"))
        if key != cur_key:
            cur = {}
            events.append((r.get("step"), float(r["t"]), cur))
            cur_key = key
        cur[int(r["basis_index"])] = r
    return events, reader.fieldnames


def cval(row, name):
    return complex(float(row[name + "_re"]), float(row[name + "_im"]))


def width_floor(row):
    """lambda_min of Re(A + B) for the N=2 diagonal-B convention.
    A is symmetric (A00,A01,A11), B is diagonal (B00,B11)."""
    a00 = float(row["A00_re"]); a01 = float(row["A01_re"]); a11 = float(row["A11_re"])
    b00 = float(row["B00_re"]); b11 = float(row["B11_re"])
    m00 = a00 + b00
    m11 = a11 + b11
    m01 = a01  # B off-diagonal is zero by convention
    tr = m00 + m11
    det = m00 * m11 - m01 * m01
    disc = max(0.0, tr * tr / 4.0 - det)
    return tr / 2.0 - math.sqrt(disc)


def r_real_vec(row, N=2):
    return [float(row[f"R{a}_re"]) for a in range(N)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="run directory containing snapshots.csv, or the "
                                 "snapshots.csv path itself")
    ap.add_argument("--u-on", type=float, default=0.03,
                    help="occupancy threshold on |u_k| (default 0.03)")
    ap.add_argument("--r-tol", type=float, default=1.0,
                    help="max ||Re R_k(t)-Re R_k(0)|| while parked (default 1.0)")
    ap.add_argument("--w-tol", type=float, default=0.25,
                    help="max relative width drift while parked (default 0.25)")
    args = ap.parse_args()

    path = args.path
    if os.path.isdir(path):
        path = os.path.join(path, "snapshots.csv")
    if not os.path.exists(path):
        raise SystemExit(f"not found: {path}")

    events, fields = load_snapshots(path)
    if "A01_re" not in fields:
        raise SystemExit("this checker expects the N=2 snapshots schema "
                         "(A00/A01/A11,B00/B11,R0/R1); got: " + ",".join(fields))

    step0, t0, ev0 = events[0]
    indices = sorted(ev0.keys())

    print(f"snapshots : {path}")
    print(f"events    : {len(events)}  (t {events[0][1]:.4f} -> {events[-1][1]:.4f})")
    print(f"terms     : {len(indices)}   u_on={args.u_on}  r_tol={args.r_tol}  "
          f"w_tol={args.w_tol}")
    print()

    reserves = []
    for k in indices:
        u0 = abs(cval(ev0[k], "u"))
        if u0 < args.u_on:
            reserves.append(k)

    print(f"reserves (|u(0)|<u_on): {len(reserves)} of {len(indices)}  -> "
          f"{reserves}")
    print()

    hdr = (f"{'term':>4} {'|u0|':>10} {'|u_end|':>10} {'parked_evts':>11} "
           f"{'max|dR|':>10} {'max_dW_rel':>11} {'verdict':>8}")
    print(hdr)
    print("-" * len(hdr))

    overall_ok = True
    any_reserve_moved = False
    for k in reserves:
        u0 = abs(cval(ev0[k], "u"))
        r0 = r_real_vec(ev0[k])
        w0 = width_floor(ev0[k])
        max_dr = 0.0
        max_dw = 0.0
        parked_evts = 0
        for (_, _, ev) in events:
            if k not in ev:
                continue
            ut = abs(cval(ev[k], "u"))
            if ut >= args.u_on:
                # term has been occupied; parking no longer required from here on
                break
            parked_evts += 1
            rt = r_real_vec(ev[k])
            dr = math.sqrt(sum((rt[a] - r0[a]) ** 2 for a in range(len(r0))))
            wt = width_floor(ev[k])
            dw = abs(wt - w0) / max(abs(w0), 1e-12)
            max_dr = max(max_dr, dr)
            max_dw = max(max_dw, dw)
        u_end = abs(cval(events[-1][2][k], "u")) if k in events[-1][2] else float("nan")
        ok = (max_dr < args.r_tol) and (max_dw < args.w_tol)
        if not ok:
            overall_ok = False
        if max_dr > 1e-9:
            any_reserve_moved = True
        print(f"{k:>4} {u0:>10.3e} {u_end:>10.3e} {parked_evts:>11} "
              f"{max_dr:>10.3e} {max_dw:>11.3e} {'OK' if ok else 'FAIL':>8}")

    print()
    if not reserves:
        print("VERDICT: no reserves detected at t=0 (all terms above u_on) — "
              "nothing to check.")
        return 0
    if not any_reserve_moved:
        print("note: all reserves stayed exactly put (max|dR|~0) — consistent with "
              "a frozen gate, but confirm the run actually stepped.")
    print("VERDICT:", "PASS" if overall_ok else "FAIL")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
