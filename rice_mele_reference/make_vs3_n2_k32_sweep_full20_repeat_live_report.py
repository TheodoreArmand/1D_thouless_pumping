#!/usr/bin/env python3
"""Build the compact live HTML report for the repeated Vs3/Vl3 N=2 K32 full-cycle dt/20 sweep."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_full20_repeat_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 full-cycle dt/20 repeat live sweep"
report.REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_full20_repeat_live_report.py"
report.REPORT_NOTE = (
    " ECG data here repeat the full-cycle fixed <code>dt=0.0005</code> run "
    "in an independent output tree, without touching the original full20 run."
)


if __name__ == "__main__":
    report.main()
