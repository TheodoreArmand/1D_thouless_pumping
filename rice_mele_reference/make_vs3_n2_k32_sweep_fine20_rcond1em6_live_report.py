#!/usr/bin/env python3
"""Build the compact live HTML report for Vs3/Vl3 N=2 K32 fine20 rcond=1e-6."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_fine20_rcond1em6_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em6"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 fine20 rcond=1e-6 live sweep"
report.REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_fine20_rcond1em6_live_report.py"
report.REPORT_NOTE = (
    " ECG data here use fine20 local dt with <code>rcond=1e-6</code>. "
    "Only <code>free</code> and <code>gauss sigma=1</code> are shown."
)
report.CASES = [case for case in report.CASES if case["key"] in {"free", "gauss_sigma1p000"}]


if __name__ == "__main__":
    report.main()
