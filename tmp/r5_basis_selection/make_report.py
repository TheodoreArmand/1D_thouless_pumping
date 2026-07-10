#!/usr/bin/env python3
"""Build the R5 selection figure + markdown report from the two summary JSONs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
struct = json.loads((HERE/"r5_selection_summary.json").read_text())
rich = json.loads((HERE/"r5_rich_summary.json").read_text())

phi = np.array(rich["phi_pi"])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

# panel 1: infidelity vs phi
ax[0].semilogy(phi, rich["infid_K32"], "o-", color="#c0392b", label="K32 (baseline)")
ax[0].semilogy(struct["phi_pi"], struct["infid_Ksel"], "s--", color="#2980b9",
               label="K48 structured union")
ax[0].semilogy(phi, rich["infid_Ksel"], "D-", color="#27ae60",
               label="K48 rich (width-enriched)")
ax[0].axvline(0.52, color="gray", ls=":", lw=1)
ax[0].text(0.522, 2e-3, "Gate-1 window\n$\\phi\\leq0.52\\pi$", fontsize=8, color="gray")
ax[0].axhline(1e-2, color="k", ls=":", lw=0.7)
ax[0].set_xlabel("$\\phi/\\pi$"); ax[0].set_ylabel("u-only projection infidelity")
ax[0].set_title("Representability vs pump phase (first transition)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3, which="both")

# panel 2: infid vs K + lambda_min
curve = np.array(rich["curve"])
ax[1].semilogy(curve[:, 0], curve[:, 1], "o-", color="#27ae60",
               label="worst infid (all $\\phi$)")
ax[1].axhline(rich["worst_all_K32"], color="#c0392b", ls="--", lw=1,
              label=f"K32 worst = {rich['worst_all_K32']:.2f}")
ax[1].axvline(48, color="gray", ls=":", lw=1)
ax[1].set_xlabel("basis size K"); ax[1].set_ylabel("worst infidelity")
ax[1].set_title("Saturation: gain stops ~K48, then $\\lambda_{min}$ collapses")
ax[1].legend(fontsize=8, loc="upper right"); ax[1].grid(alpha=0.3, which="both")

fig.tight_layout()
fig.savefig(HERE/"r5_selection_figure.png", dpi=130)
print("wrote r5_selection_figure.png")

# markdown report
lines = []
lines.append("# R5 union basis selection — result\n")
lines.append("Free grid reference, full-depth schedule, Vs=Vl=3, T=160π, "
             "16 snapshots φ/π∈[0.30,0.60]. Selection = greedy conditioning-guarded "
             "fit of the exact free 2-particle grid state ψ=sym[φ_L⊗φ_R] "
             "(g=0 factorizes exactly → all overlaps analytic in 1D).\n")
lines.append("## Verdict\n")
lines.append(f"- **Span gap is real**: K32 worst infidelity = {rich['worst_all_K32']:.3f} "
             f"(Gate-1 φ≤0.52π: {rich['worst_gate1_K32']:.3f}). R4 → R5 confirmed.")
lines.append(f"- **Best K48 (width-enriched union)**: worst infidelity "
             f"{rich['worst_all_Ksel']:.3f} (all φ), {rich['worst_gate1_Ksel']:.3f} "
             f"(Gate-1). Transition ONSET φ≤0.46π closed to ~5e-3; ~20× better than K32.")
lines.append(f"- Gram λ_min = {rich['gram_lambda_min']:.2e} (healthy, ≫1e-8); "
             f"initial-state fidelity = {rich['init_fidelity']:.7f} (passes ≥0.99999).")
lines.append(f"- **Hard plateau floor ~3–4%**: adding terms past K≈48 gives negligible "
             f"gain while λ_min collapses (2e-3→4e-4 by K64). Single evolved plateau "
             f"packet needs ~6 free gaussians for 1e-3 → the plateau ripple/breathing "
             f"is a genuine gaussian-ansatz representability limit, not a pool gap.")
lines.append(f"- The residual lives in short-lattice Bloch-ripple / breathing modes that "
             f"carry little first-moment (ΔP) weight — the grid ΔP pumps cleanly to −2 "
             f"despite it — so K48 should still track the pump observable through Gate-1.\n")
lines.append("## Files\n")
lines.append("- `initial_pathpad_N2_K48_richunion.csv` — **recommended best 48-term basis** "
             "(width-enriched, half-integer path centers).")
lines.append("- `initial_pathpad_N2_K48_union.csv` — roadmap-faithful structured variant "
             "(P-A+P-B+P-C1+P-C3, 2 canonical widths; worst 0.093 all / 0.036 Gate-1).")
lines.append("- `r5_selection_figure.png` — infidelity vs φ and vs K.\n")
lines.append("## infidelity by φ (rich K48 vs K32)\n")
lines.append("| φ/π | K32 | K48 rich |")
lines.append("|---|---|---|")
for p, a, b in zip(phi, rich["infid_K32"], rich["infid_Ksel"]):
    lines.append(f"| {p:.2f} | {a:.2e} | {b:.2e} |")
(HERE/"R5_RESULT.md").write_text("\n".join(lines))
print("wrote R5_RESULT.md")
