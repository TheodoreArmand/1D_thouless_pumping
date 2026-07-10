#!/usr/bin/env python3
"""Build the compact live HTML report for K32 fine20 COMdiag rcond=1e-5."""
from __future__ import annotations

import make_vs3_n2_k32_sweep_live_report as report


report.OUT_HTML = report.REPO / "vs3_n2_k32_sweep_fine20_COMdiag_rcond1em5_live_report.html"
report.ECG_ROOT = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em5"
ordinary_grid_root = report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref"
report.GRID_ROOT = (
    report.REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep" / "grid_ref_comdiag"
)
report.GRID_FALLBACK_ROOT = ordinary_grid_root
report.GRID_AUDIT_ROOT = ordinary_grid_root
report.GRID_LOG_DIR = (
    report.REPO
    / "slurm"
    / "vs3_n2_dt0p01_T160pi_K32_sweep"
    / "grid_ref_comdiag"
    / "ecg1d_vs3_n2_grid_comdiag_20260710"
)
report.GRID_LOG_PREFIX = "vs3n2gcd"
report.GRID_PRIMARY_KIND = "COMdiag-matched K32 grid CSV"
report.GRID_LOG_KIND = "COMdiag-matched grid Slurm log"
report.GRID_FALLBACK_KIND = "ordinary-basis K32 grid fallback — matched pending"
report.GRID_AUDIT_KIND = "ordinary-basis K32 grid CSV"
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 fine20 COMdiag rcond=1e-5 live sweep"
report.REPORT_COMMAND = (
    "python3 rice_mele_reference/make_vs3_n2_k32_sweep_fine20_COMdiag_rcond1em5_live_report.py"
)
report.REPORT_NOTE = (
    " ECG and the primary grid reference use <code>initial_pathpad_N2_K32_COMdiag.csv</code>. "
    "The grid uses uniform <code>dt=0.01</code>; ECG uses <code>dt=0.01</code> and local "
    "<code>dt/20</code> windows "
    "<code>s=0.1-0.4</code> and <code>s=0.6-0.9</code>. "
    r"Only free and Gaussian \(\sigma=1\) are shown. "
    "Initial-state audit: fidelity with the ordinary K32 state is 0.999939422; "
    r"\(P(0)\) changes from 0.999658891 to 0.999685570, "
    r"\(r_{12}(0)\) from 8.045593816 to 8.045508451, and Gaussian "
    r"\(V_g(0)\) from \(7.906531\times10^{-7}\) to \(2.868864\times10^{-11}\)."
)
report.GRID_NOTE = (
    "The selected reference is the COMdiag-matched 1024-grid CSV when complete, "
    "then its live Slurm log while running; the ordinary-basis K32 grid is used only "
    "until matched data exist and remains visible in the reference audit."
)
report.CASES = [
    {
        "key": "free",
        "title": "COMdiag free",
        "mode": "free",
        "sigma": None,
        "grid_index": 0,
        "grid_log_index": 0,
        "color": "#2267a8",
    },
    {
        "key": "gauss_sigma1p000",
        "title": "COMdiag gauss sigma=1",
        "mode": "gauss",
        "sigma": 1.0,
        "grid_index": 3,
        "grid_log_index": 1,
        "color": "#dc2626",
    },
]


if __name__ == "__main__":
    report.main()
