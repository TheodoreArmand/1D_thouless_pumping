#!/usr/bin/env python3
"""Build the compact live HTML report for K32 fine20 COMdiag rcond=1e-5."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_fine20_COMdiag_rcond1em5_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em5"
report.GRID_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 fine20 COMdiag rcond=1e-5 live sweep"
report.REPORT_COMMAND = (
    "python3 rice_mele_reference/make_vs3_n2_k32_sweep_fine20_COMdiag_rcond1em5_live_report.py"
)
report.REPORT_NOTE = (
    " ECG data use <code>initial_pathpad_N2_K32_COMdiag.csv</code>, with "
    "<code>dt=0.01</code> and local <code>dt/20</code> windows "
    "<code>s=0.1-0.4</code> and <code>s=0.6-0.9</code>. "
    "Only free and Gaussian <code>sigma=1</code> are shown."
)
report.CASES = [
    {
        "key": "free",
        "title": "COMdiag free",
        "mode": "free",
        "sigma": None,
        "grid_index": 0,
        "color": "#2267a8",
    },
    {
        "key": "gauss_sigma1p000",
        "title": "COMdiag gauss sigma=1",
        "mode": "gauss",
        "sigma": 1.0,
        "grid_index": 3,
        "color": "#dc2626",
    },
]


if __name__ == "__main__":
    report.main()
