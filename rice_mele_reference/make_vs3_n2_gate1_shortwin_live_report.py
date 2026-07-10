#!/usr/bin/env python3
"""Live report for the R3a Gate-1 short-window occupancy-gate experiment.

Monitors out/vs3_n2_gate1_shortwin_rcond1em5/{gaussFrozenA,gaussA01,free}_{gate,
nogate}/... (driver ecg1d_vs3_n2_gate1_shortwin, full-depth schedule) against the
matching full-depth grid reference. Reuses the sweep report machinery (transport /
delta-Delta P / health) and adds an R3a occupancy-gate panel per case built from
occupancy_gate.csv: the global discarded-RHS fraction, the parked B+R+A share of
that discarded mass (the H1 signature), and the active param-dim / occupancy
census. Compare each (gate vs nogate) pair.

  python3 rice_mele_reference/make_vs3_n2_gate1_shortwin_live_report.py [--watch S]
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import make_vs3_n2_k32_sweep_live_report as report

REPO = report.REPO
RCOND = "rcond1em5"  # the submitted primary set (array 0-5)

# --- retarget the base module at the Gate-1 layout + full-depth grid ref -------
report.OUT_HTML = REPO / "vs3_n2_gate1_shortwin_live_report.html"
report.SCHEDULE_CSV = (
    REPO / "rice_mele_reference" / "Vs3Vl3_3_3"
    / "gap_adaptive_vs3vl3_full_depth_schedule.csv"
)
report.ECG_ROOT = REPO / "out" / ("vs3_n2_gate1_shortwin_" + RCOND)
report.ECG_TASK_SUFFIX = "_full_depth_schedule"
report.GRID_ROOT = (
    REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule" / "grid_ref"
)
report.GRID_FALLBACK_ROOT = None
report.GRID_AUDIT_ROOT = None
report.ALLOW_OLD_GRID_FALLBACK = False
report.GRID_LOG_DIR = (
    REPO / "slurm" / "vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule"
    / "grid_ref" / "ecg1d_vs3_n2_grid_k32_sweep_20260709"
)
report.REPORT_TITLE = "Vs3/Vl3 N=2 K32 - R3a Gate-1 occupancy gate (rcond 1e-5, full-depth)"
report.PROTOCOL_NOTE = (
    "<strong>R3a Gate-1 short window</strong> (s in [0, 0.35], fine dt/10 in "
    "[0.18, 0.32] around the gap minimum at s=0.25). Each mode runs with the "
    "occupancy gate ON and OFF; compare the pairs. Full-depth Hamiltonian + "
    "matching full-depth phase schedule."
)
report.REPORT_COMMAND = (
    "python3 rice_mele_reference/make_vs3_n2_gate1_shortwin_live_report.py"
)
report.REPORT_NOTE = (
    " Gauss runs are frozen-A (R1/E1) and A_01-only (R2). The occupancy panel is "
    "the R3a readout: the gate should collapse the parked B+R+A discarded share "
    "toward 0 and lower the global discarded fraction at the transition."
)
report.GRID_NOTE = (
    "Grid reference from <code>..._full_depth_schedule/grid_ref</code> (free and "
    "gauss_sigma1p000); the gauss/A01 ECG cases are compared against the "
    "gauss_sigma1p000 grid, the free cases against the free grid."
)

TASK_DIR = "a8p000_K32_tmax175p929_VsER3p000_VlER3p000" + report.ECG_TASK_SUFFIX


# --- path overrides: the Gate-1 dirs differ from the full-cycle sweep layout ---
def ecg_progress_path(case: dict) -> Path:
    return report.ECG_ROOT / case["ecg_subdir"] / TASK_DIR / "progress.csv"


def grid_csv_path(case: dict) -> Path:
    # gauss frozen-A and A_01 both compare against the gauss_sigma1p000 grid;
    # free compares against the free grid.
    if case["mode"] == "free":
        gcase = {"mode": "free"}
    else:
        gcase = {"mode": "gauss", "key": "gauss_sigma1p000"}
    return report.grid_csv_path_at(report.GRID_ROOT, gcase)


report.ecg_progress_path = ecg_progress_path
report.grid_csv_path = grid_csv_path


# --- R3a occupancy-gate panel from occupancy_gate.csv --------------------------
def occupancy_plot(case: dict) -> str | None:
    path = report.ECG_ROOT / case["ecg_subdir"] / TASK_DIR / "occupancy_gate.csv"
    if not path.exists():
        return None
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None

    def col(name: str) -> list[float]:
        return [float(r[name]) for r in rows]

    phi = [float(r["phi"]) / math.pi for r in rows]
    gfrac = col("disc_frac_global")
    parked = [
        float(r["share_B_parked"]) + float(r["share_R_parked"]) + float(r["share_A_parked"])
        for r in rows
    ]
    active = col("active_param_dim")
    n_core = col("n_core")
    n_parked = col("n_parked")

    plt = report.plt
    fig, axes = plt.subplots(2, 1, figsize=(3.18, 2.25), sharex=True)
    axes[0].plot(phi, gfrac, color=case["color"], lw=1.15, label="global disc frac")
    axes[0].plot(phi, parked, color="#b91c1c", lw=1.0, ls="--", label="parked B+R+A share")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].set_ylabel("disc", fontsize=7.2)
    axes[0].legend(frameon=False, fontsize=6.0, loc="upper left", handlelength=1.6)
    axes[0].set_title("occupancy gate", fontsize=8.4)

    axes[1].plot(phi, active, color="#334155", lw=1.1, label="active dim")
    ax2 = axes[1].twinx()
    ax2.plot(phi, n_core, color="#2563eb", lw=0.9, label="n_core")
    ax2.plot(phi, n_parked, color="#059669", lw=0.9, ls=":", label="n_parked")
    ax2.tick_params(colors="#475569", labelsize=6.6, width=0.6)
    for spine in ax2.spines.values():
        spine.set_color("#cbd5e1")
        spine.set_linewidth(0.7)
    axes[1].set_ylabel("dim", fontsize=7.2)
    axes[1].set_xlabel("phi / pi", fontsize=7.5)

    for ax in axes:
        ax.set_xlim(0.0, 0.6)
        report.style_axis(ax)
    fig.tight_layout(h_pad=0.2)
    return report.fig_uri(fig)


_orig_build_case = report.build_case
_orig_case_section = report.case_section


def build_case(case: dict) -> dict:
    item = _orig_build_case(case)
    item["occupancy_plot"] = occupancy_plot(case)
    return item


def case_section(item: dict) -> str:
    html = _orig_case_section(item)
    occ = item.get("occupancy_plot")
    if occ:
        fig_html = (
            '<figure><img src="' + occ + '" alt="' + item["title"] + ' occupancy">'
            "<figcaption>R3a occupancy gate: global discarded fraction &amp; parked "
            "B+R+A share (top); active dim &amp; core/parked census (bottom).</figcaption></figure>"
        )
        html = html.replace(
            "</figure>\n      </div>", "</figure>" + fig_html + "\n      </div>", 1
        )
    return html


report.build_case = build_case
report.case_section = case_section


# --- cases: gate vs nogate for each of gauss frozen-A / gauss A01 / free --------
report.CASES = [
    {"key": "gaussFrozenA_nogate_sigma1p000", "ecg_subdir": "gaussFrozenA_nogate_sigma1p000",
     "title": "gauss frozen-A - nogate", "mode": "gauss", "sigma": 1.0, "grid_index": 3, "color": "#f59e0b"},
    {"key": "gaussFrozenA_gate_sigma1p000", "ecg_subdir": "gaussFrozenA_gate_sigma1p000",
     "title": "gauss frozen-A - GATE", "mode": "gauss", "sigma": 1.0, "grid_index": 3, "color": "#c2410c"},
    {"key": "gaussA01_nogate_sigma1p000", "ecg_subdir": "gaussA01_nogate_sigma1p000",
     "title": "gauss A01-only - nogate", "mode": "gauss", "sigma": 1.0, "grid_index": 3, "color": "#60a5fa"},
    {"key": "gaussA01_gate_sigma1p000", "ecg_subdir": "gaussA01_gate_sigma1p000",
     "title": "gauss A01-only - GATE", "mode": "gauss", "sigma": 1.0, "grid_index": 3, "color": "#1d4ed8"},
    {"key": "free_nogate", "ecg_subdir": "free_nogate",
     "title": "free - nogate", "mode": "free", "sigma": None, "grid_index": 0, "color": "#86efac"},
    {"key": "free_gate", "ecg_subdir": "free_gate",
     "title": "free - GATE", "mode": "free", "sigma": None, "grid_index": 0, "color": "#15803d"},
]


if __name__ == "__main__":
    report.main()
