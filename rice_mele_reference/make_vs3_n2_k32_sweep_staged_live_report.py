#!/usr/bin/env python3
"""Build the compact live HTML report for the Vs3/Vl3 N=2 K32 staged-dt sweep."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_staged_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_staged"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 staged-dt live sweep"
report.REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_staged_live_report.py"
report.REPORT_NOTE = (
    " ECG data here use staged local dt: <code>dt/20</code> in "
    "<code>s in [0.20,0.30] U [0.70,0.80]</code>, "
    "<code>dt/10</code> in shoulders "
    "<code>[0.18,0.20] U [0.30,0.32] U [0.68,0.70] U [0.80,0.82]</code>, "
    "and base <code>dt</code> elsewhere."
)


if __name__ == "__main__":
    report.main()
