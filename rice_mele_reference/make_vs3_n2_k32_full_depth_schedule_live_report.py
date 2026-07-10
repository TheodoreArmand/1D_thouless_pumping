#!/usr/bin/env python3
"""Build the live report for the unified full-depth Vs3/Vl3 protocol."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_full_depth_schedule_live_report.html"
report.SCHEDULE_CSV = (
    report.REPO
    / "rice_mele_reference"
    / "Vs3Vl3_3_3"
    / "gap_adaptive_vs3vl3_full_depth_schedule.csv"
)
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep"
report.ECG_TASK_SUFFIX = "_full_depth_schedule"
report.GRID_ROOT = (
    report.REPO
    / "out"
    / "vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule"
    / "grid_ref"
)
report.GRID_FALLBACK_ROOT = None
report.GRID_AUDIT_ROOT = None
report.ALLOW_OLD_GRID_FALLBACK = False
report.GRID_LOG_DIR = (
    report.REPO
    / "slurm"
    / "vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule"
    / "grid_ref"
    / "ecg1d_vs3_n2_grid_k32_sweep_20260709"
)
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 unified full-depth schedule"
report.PROTOCOL_NOTE = (
    "<strong>Current protocol:</strong> full-depth Hamiltonian driven by the "
    "matching full-depth-derived phase schedule."
)
report.REPORT_COMMAND = (
    "python3 rice_mele_reference/"
    "watch_vs3_n2_k32_full_depth_schedule_live_report.py"
)
report.REPORT_NOTE = (
    " This report accepts only the unified protocol: full-depth Hamiltonian and "
    "<code>gap_adaptive_vs3vl3_full_depth_schedule.csv</code>. Historical hybrid "
    "ECG/grid outputs are intentionally excluded."
)
report.GRID_NOTE = (
    "The grid reference must come from the isolated "
    "<code>..._full_depth_schedule/grid_ref</code> root; no legacy fallback is used."
)


if __name__ == "__main__":
    report.main()
