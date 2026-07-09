#!/usr/bin/env python3
"""Build the compact live HTML report for the Vs3/Vl3 N=2 K32 dt/20 sweep."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_fine20_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_fine20"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 dt/20 live sweep"
report.REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_fine20_live_report.py"
report.REPORT_NOTE = (
    " ECG data here use local fine-step factor <code>0.05</code>, i.e. "
    "<code>dt/20</code> inside <code>s in [0.2,0.3] U [0.7,0.8]</code>."
)


if __name__ == "__main__":
    report.main()
