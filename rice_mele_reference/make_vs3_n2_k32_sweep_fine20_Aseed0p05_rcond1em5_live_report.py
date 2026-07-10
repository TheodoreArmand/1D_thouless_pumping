#!/usr/bin/env python3
"""Build the compact live HTML report for K32 fine20 Aseed0p05 rcond=1e-5."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_fine20_Aseed0p05_rcond1em5_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_fine20_Aseed0p05_rcond1em5"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 fine20 Aseed0p05 rcond=1e-5 live sweep"
report.REPORT_COMMAND = (
    "python3 rice_mele_reference/make_vs3_n2_k32_sweep_fine20_Aseed0p05_rcond1em5_live_report.py"
)
report.REPORT_NOTE = (
    " ECG data here use fine20 local dt, <code>rcond=1e-5</code>, and "
    "<code>initial_pathpad_N2_K32_Aseed0p05.csv</code>. "
    "The 16 core terms keep <code>A=0</code>; the 16 path-pad terms use "
    "<code>A=0.05[[1,-1],[-1,1]]</code> with unchanged B. "
    "Only <code>free</code> and <code>gauss sigma=1</code> are shown."
)
report.CASES = [case for case in report.CASES if case["key"] in {"free", "gauss_sigma1p000"}]


if __name__ == "__main__":
    report.main()
