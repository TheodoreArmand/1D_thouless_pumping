#!/usr/bin/env python3
"""Refresh the Vs3/Vl3 N=2 K32 live reports on a fixed interval."""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
REPORTS = [
    REPO / "rice_mele_reference" / "make_vs3_n2_k32_sweep_live_report.py",
    REPO / "rice_mele_reference" / "make_vs3_n2_k32_sweep_fine20_live_report.py",
]


def stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def refresh_once() -> int:
    worst = 0
    for script in REPORTS:
        print(f"[{stamp()}] refresh {script.relative_to(REPO)}", flush=True)
        proc = subprocess.run([sys.executable, str(script)], cwd=REPO)
        worst = max(worst, proc.returncode)
        print(f"[{stamp()}] return_code={proc.returncode}", flush=True)
    return worst


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interval", type=float, default=600.0)
    parser.add_argument(
        "--iterations",
        type=int,
        default=0,
        help="Number of refresh cycles; 0 means run until the Slurm time limit.",
    )
    args = parser.parse_args()

    count = 0
    last_rc = 0
    while args.iterations <= 0 or count < args.iterations:
        count += 1
        print(f"[{stamp()}] cycle={count} start", flush=True)
        last_rc = refresh_once()
        print(f"[{stamp()}] cycle={count} done last_rc={last_rc}", flush=True)
        if args.iterations > 0 and count >= args.iterations:
            break
        time.sleep(max(args.interval, 1.0))
    return last_rc


if __name__ == "__main__":
    raise SystemExit(main())
