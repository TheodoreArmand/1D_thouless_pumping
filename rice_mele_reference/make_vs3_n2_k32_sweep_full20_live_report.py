#!/usr/bin/env python3
"""Build the compact live HTML report for the Vs3/Vl3 N=2 K32 full-cycle dt/20 sweep."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_full20_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_full20"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 full-cycle dt/20 live sweep"
report.REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_full20_live_report.py"
report.REPORT_NOTE = (
    " ECG data here use fixed <code>dt=0.0005</code> for the full cycle, "
    "i.e. global <code>dt/20</code> relative to the base <code>dt=0.01</code>."
)


if __name__ == "__main__":
    report.main()
