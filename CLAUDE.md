# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

C++17 simulation of 1D Thouless pumping in an optical superlattice using Explicitly Correlated Gaussians (ECG) evolved by real-time TDVP. The wavefunction is a sum of K shifted complex Gaussians ψ_k = u_k·exp(−xᵀA_k x − (x−R_k)ᵀB_k(x−R_k)), (anti)symmetrized over the N! particle permutations; the parameters z = (u, B, R, A) evolve by i·C·ż = ∂⟨H⟩ with a Tikhonov + SVD-truncated linear solve and RK4 integration. Current physics target: the legacy-PRB pump H(t) = p²/2m − Vs·cos(4πx/a) − Vl·cos(2πx/a + φ(t)) with N=1, K=16, Vs=Vl=3 Er, T=160π, gap-adaptive φ(t) schedule.

## Build and run

Requires Eigen3 and OpenMP. There are no tests and no lint.

```bash
cmake -S . -B build              # configure once
cmake --build build -j           # → build/ecg1d_thouless_pumping
./build/ecg1d_thouless_pumping   # run from the repo root — config paths are relative
```

- The executable is deliberately **not a CLI**. Every parameter is hard-coded in `pumpconfig/legacy_prb_3_3.cpp` (registered by name in `pumpconfig/pump_config.cpp`). To change a run: edit that config or add a new named config, then rebuild. For a smoke test, shorten `total_time`/`pump_period` there.
- **A full default run takes ~41 hours** (251k fixed-dt RK4 steps, ~0.6 s/step). Never launch one casually.
- **Re-running the default config overwrites `out/pump_vs3pad_gapadaptive_T160pi/a8p000_K16_.../`** — the archived output of a successful 41-hour run. `out/` is gitignored, so that data exists nowhere else. Change `out_root` (or move the old dir) before experimental runs.
- Thread count: `OMP_NUM_THREADS`. `EIGEN_DONT_PARALLELIZE` is set intentionally — OpenMP parallelizes the metric/gradient assembly loops in `src/tdvp_solver.cpp`, so Eigen must not nest threads.
- Precedent: the binary used for a specific run is archived as a timestamped copy in `build/` (e.g. `ecg1d_thouless_pumping_short_snapshot_20260707_173458`).

## Architecture

One executable = ECG/TDVP engine (`src/`, namespace `ecg1d`) + run configs (`pumpconfig/`) + driver (`main.cpp`).

`src/` is layered bottom-up:

1. **Matrix elements**: `permutation` (N! permutation set with signs) → `pair_cache` (per ⟨bra i | ket j, permutation p⟩ cached K, K⁻¹, μ, det, overlap M_G) → `interaction_kernels` (kinetic P_Mij, delta-contact G, Gaussian H, one-body cosine kernels) → `hamiltonian` (sums kernels over i,j,p into ⟨S⟩ and ⟨H⟩ functionals).
2. **Analytic gradients**: `derivatives` (12-case derivative tables of log-overlap and K⁻¹ w.r.t. u/A/B/R — the largest and most delicate file), `observable_derivatives`, `hamiltonian_gradient` (∂(⟨H⟩/S)/∂z per term).
3. **TDVP**: `tdvp_solver` — assembles metric C and gradient g over an `AlphaIndex` list (the OpenMP hot loops), `HamiltonianTerms` toggles kinetic/delta/gaussian + arbitrary cosine terms. `realtime_tdvp` — `compute_rhs_dz` solves (C + λ_C·I)dz = −i·g via SVD with relative `rcond` truncation, records per-step linear-algebra diagnostics (`RealtimeStepResult`); Euler/RK4/RK45 integrators plus the time-dependent-H RK4 (`realtime_tdvp_step_rk4_time_dependent`) that the pump driver uses.
4. **I/O**: `snapshot_io` (in-memory-buffered snapshots.csv; writer currently assumes N=1), `trace_io` (observables + solver-health trace), `run_report` (config.txt / summary.txt), `csv_utils`.

`main.cpp` flow: load gap-adaptive phase schedule CSV and initial basis CSV → fixed-dt RK4 TDVP loop with per-step norm re-enforcement and finiteness guards → write trace/snapshots/basis/summary into a parameter-tagged output dir.

Key convention: TDVP parameters are addressed by `AlphaIndex{a1,..}` with a1 = 1:u, 2:B, 3:R, 4:A. `make_alpha_list` in main.cpp deliberately **freezes A** (for N=1, A and B are both width-like directions and together make C ill-conditioned), so param_dim = K·(1+2N).

## Units and conventions (gotchas)

- Code units ħ = 1, m = 1 (`src/physical_constants.hpp`); `recoil_energy_for_lattice` defines Er with k = 2π/lattice_a (long-lattice period a = 8), while the Rice–Mele reference notes use short-lattice units (d_s = 1, E_{r,s} = 1).
- **Factor-of-2 depth trap**: "Vs=Vl=3" labels two different potentials — the full-depth "grid" convention vs the half-depth "main.cpp" convention that the ECG binary actually integrates. Read `rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_Vs3Vl3_convention_note.md` before comparing any ECG result to reference curves or regenerating schedules.
- Phase schedule CSV: columns `s,phi` with s = t/T strictly increasing and spanning exactly [0,1]; linearly interpolated; `#` comments allowed.
- Basis CSV format (`initial_state/`, `basis_initial/final.csv`): per Gaussian, one `re,im` line each for u, then A (N² entries), B (N²), R (N), then one integer name line.

## Outputs and how a run is judged

Each run writes `out/<out_root>/a<a>_K<K>_tmax<T>_VsER<vs>_VlER<vl>/` containing `config.txt`, `trace.csv`, `snapshots.csv`, `basis_initial.csv`, `basis_final.csv`, `summary.txt`. Success criteria for a pump run: Δpolarization_cell ≈ −1 per cycle, small norm/energy drift, and healthy solver diagnostics — max actual-solve condition number, effective rank vs param_dim, discarded-rhs fraction, min Re(B) and min Re(A+B) staying positive. These are printed per sampled step and aggregated in summary.txt.

Reference/analysis material: `rice_mele_reference/Vs3Vl3_3_3/` (phase schedules, Rice–Mele band-model benchmark data + README), `successful_case/Vs3Vl3_3_3/` (Jupyter density-movie notebook and HTML report over the good run's snapshots.csv), `reference/` (papers).

## Notes

- Project documentation and notebook narration mix Chinese and English.
- `slurm/` and `txt_task/` are empty placeholders (cluster scripts, task notes).
