#!/usr/bin/env python3
"""Continuously refresh the N=2 K32 full-depth-schedule live report.

Run from the repository root with no arguments for a 30-second cadence:

    python3 rice_mele_reference/watch_vs3_n2_k32_full_depth_schedule_live_report.py

Press Ctrl-C to stop. Use ``--iterations 1`` for a single refresh.
"""
from __future__ import annotations

import argparse
import time
from datetime import datetime

import make_vs3_n2_k32_full_depth_schedule_live_report as configured


def stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--interval",
        type=float,
        default=30.0,
        help="Target seconds between the start of consecutive refreshes (default: 30).",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=0,
        help="Number of refreshes; 0 means continue until Ctrl-C (default: 0).",
    )
    args = parser.parse_args()
    if args.interval < 1.0:
        parser.error("--interval must be at least 1 second")
    if args.iterations < 0:
        parser.error("--iterations must be nonnegative")

    count = 0
    try:
        while args.iterations == 0 or count < args.iterations:
            count += 1
            cycle_started = time.monotonic()
            output = configured.report.write_report()
            elapsed = time.monotonic() - cycle_started
            print(
                f"[{stamp()}] refresh={count} elapsed_s={elapsed:.2f} output={output}",
                flush=True,
            )
            if args.iterations > 0 and count >= args.iterations:
                break
            time.sleep(max(0.0, args.interval - elapsed))
    except KeyboardInterrupt:
        print(f"\n[{stamp()}] stopped after {count} refreshes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
