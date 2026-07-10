# R5 union basis selection — result

Free grid reference, full-depth schedule, Vs=Vl=3, T=160π, 16 snapshots φ/π∈[0.30,0.60]. Selection = greedy conditioning-guarded fit of the exact free 2-particle grid state ψ=sym[φ_L⊗φ_R] (g=0 factorizes exactly → all overlaps analytic in 1D).

## Verdict

- **Span gap is real**: K32 worst infidelity = 0.968 (Gate-1 φ≤0.52π: 0.668). R4 → R5 confirmed.
- **Best K48 (width-enriched union)**: worst infidelity 0.045 (all φ), 0.032 (Gate-1). Transition ONSET φ≤0.46π closed to ~5e-3; ~20× better than K32.
- Gram λ_min = 2.08e-03 (healthy, ≫1e-8); initial-state fidelity = 0.9999931 (passes ≥0.99999).
- **Hard plateau floor ~3–4%**: adding terms past K≈48 gives negligible gain while λ_min collapses (2e-3→4e-4 by K64). Single evolved plateau packet needs ~6 free gaussians for 1e-3 → the plateau ripple/breathing is a genuine gaussian-ansatz representability limit, not a pool gap.
- The residual lives in short-lattice Bloch-ripple / breathing modes that carry little first-moment (ΔP) weight — the grid ΔP pumps cleanly to −2 despite it — so K48 should still track the pump observable through Gate-1.

## Files

- `initial_pathpad_N2_K48_richunion.csv` — **recommended best 48-term basis** (width-enriched, half-integer path centers).
- `initial_pathpad_N2_K48_union.csv` — roadmap-faithful structured variant (P-A+P-B+P-C1+P-C3, 2 canonical widths; worst 0.093 all / 0.036 Gate-1).
- `r5_selection_figure.png` — infidelity vs φ and vs K.

## infidelity by φ (rich K48 vs K32)

| φ/π | K32 | K48 rich |
|---|---|---|
| 0.30 | 8.01e-03 | 3.93e-03 |
| 0.32 | 1.02e-02 | 4.55e-03 |
| 0.34 | 1.05e-02 | 4.30e-03 |
| 0.36 | 1.43e-02 | 4.77e-03 |
| 0.38 | 1.36e-02 | 4.66e-03 |
| 0.40 | 1.79e-02 | 5.07e-03 |
| 0.42 | 2.47e-02 | 5.15e-03 |
| 0.44 | 2.87e-02 | 5.22e-03 |
| 0.46 | 5.51e-02 | 5.88e-03 |
| 0.48 | 1.33e-01 | 8.09e-03 |
| 0.50 | 3.62e-01 | 1.68e-02 |
| 0.52 | 6.68e-01 | 3.17e-02 |
| 0.54 | 8.59e-01 | 3.97e-02 |
| 0.56 | 9.23e-01 | 4.29e-02 |
| 0.58 | 9.54e-01 | 4.41e-02 |
| 0.60 | 9.68e-01 | 4.47e-02 |