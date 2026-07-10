# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

C++17 simulation of 1D Thouless pumping in an optical superlattice using Explicitly Correlated Gaussians (ECG) evolved by real-time TDVP. The wavefunction is a sum of K shifted complex Gaussians ψ_k = u_k·exp(−xᵀA_k x − (x−R_k)ᵀB_k(x−R_k)), (anti)symmetrized over the N! particle permutations; the parameters z = (u, B, R, A) evolve by i·C·ż = ∂⟨H⟩ with a Tikhonov + SVD-truncated linear solve and RK4 integration. Physics: the legacy-PRB pump H(t) = p²/2m − Vs·cos(4πx/a) − Vl·cos(2πx/a + φ(t)) with Vs=Vl=3 Er, T=160π, gap-adaptive φ(t) schedule.

**Current campaign: full-cycle N=2 pumping at K=32** (two bosons, path-pad basis, optionally with Gaussian pair repulsion g=0.3 E_rs), run as Slurm sweeps over sigma_gauss, `rcond`, dt refinement, and basis variants, monitored through auto-refreshed HTML reports. The archived N=1 K=16 run is the baseline. PairCache memoization (see Architecture) made full N=2 cycles practical: ~0.4–0.6 s/step at K=32 on 32 cores → roughly 1–2 days per full T=160π run (~240–280k steps).

## Build and run

Requires Eigen3 and OpenMP. There are no tests and no lint (only standalone validation targets, below).

```bash
cmake -S . -B build              # configure once
cmake --build build -j           # builds all targets
./build/ecg1d_thouless_pumping   # run from the repo root — config paths are relative
```

Targets fall into three groups (all share `ECG_SRC` + `PUMPCONFIG_SRC` in CMakeLists.txt):

- **Core**: `ecg1d_thouless_pumping` (N=1, `main.cpp`, deliberately not a CLI — every parameter is hard-coded in `pumpconfig/legacy_prb_3_3.cpp`, registered by name in `pumpconfig/pump_config.cpp`); `ecg1d_2ninteraction` / `ecg1d_2gaussian` (short N=2 validation windows, `total_time=12`).
- **Task drivers**: one `txt_task/<name>.cpp` per experiment, each registered as its own CMake target (copy an existing `add_executable` block — they are all identical boilerplate). A driver builds a `PumpConfig` in code (usually starting from `make_pump_config("legacy_prb_3_3")` and overriding N/K/dt/rcond/basis/out_root), then calls `pump2::run_pump(cfg, opt)`. Unlike `main.cpp`, drivers may take positional args for the sweep case, e.g. `./build/ecg1d_vs3_n2_dt001_fine20_aseed_rcond rcond1em5 gauss 1.0`. Naming decode: `fine5`/`fine20` = dt divided by 5/20 inside the fine windows, `dt0005` = uniform dt=5e-4, `staged` = per-window `ManualFineDtWindow` schedule, `Aseed0p05` = A-seeded basis, `rcond1emX` = rcond=1e−X.
- **Validation/bench**: `ecg1d_cache_hermitian_validate`, `ecg1d_gauss_gradient_check`, `ecg1d_vs3_n2_k32_fullstep_threads_bench`, `ecg1d_legacy_prb_3_3_300step` (N=1 smoke/bench, ~0.045 s/step on 8 threads post-memoization).

Every driver appends `key=value` provenance lines (`opt.config_appendix`) to the run's `config.txt` — sweep parameter, basis variant, fine-dt settings, rcond case, wall-time planning estimate. Keep doing this for new drivers.

- **Re-running the N=1 default config overwrites `out/pump_vs3pad_gapadaptive_T160pi/a8p000_K16_.../`** — the archived output of a successful long run. `out/` is gitignored, so that data exists nowhere else. Change `out_root` (or move the old dir) before experimental runs.
- Thread count: `OMP_NUM_THREADS`. `EIGEN_DONT_PARALLELIZE` is set intentionally — OpenMP parallelizes the metric/gradient assembly loops in `src/tdvp_solver.cpp`, so Eigen must not nest threads.
- Precedent: the binary used for a specific run is archived as a timestamped copy in `build/`.

## Running on the cluster (Slurm)

Anything longer than a smoke test goes through Slurm — never launch a multi-hour run in the foreground.

- sbatch scripts live in `txt_task/<task>_<YYYYMMDD>.sbatch`. Each starts with a comment block of copy-paste **Submit / Monitor / After completion / Cancel** commands (submit uses `sbatch --parsable ... | tee <meta_dir>/<task>.jobid`).
- Script conventions: `--chdir` to the repo root, `--exclude=bar`, explicit `--cpus-per-task/--mem/--time`, `--output/--error` and the `.jobid` file all under a mirrored metadata dir `slurm/<sweep_name>/<task_name>/`; `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` with `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores`; sweeps as `--array` over case tables (`case_modes`/`case_sigmas`/`case_labels`) inside the script; `echo` all job metadata (binary, basis, rcond, planning estimate) at start.
- Partitions: `epyc`/`zen5` for compute; the 1-cpu HTML auto-refresh jobs run on `wild`.
- `slurm/` holds only per-job metadata/log mirrors and bench results, never code. `tmp/time_summary.md` (+ `time_summary_all.csv`) is a generated progress/ETA digest across all run families.

## Live monitoring and reports

- Every run writes `progress.csv` into its task dir (`N1ProgressWriter`/`N2ProgressWriter`, `src/run_report.hpp` + `src/progress_output_n{1,2}.cpp`) and prints per-step progress lines to stdout — i.e. into the Slurm `.out`. `tail -f` the `.out` or watch `progress.csv` to monitor a run in real time.
- Live HTML reports: `rice_mele_reference/make_vs3_n2_k32_sweep_live_report.py` is the base module; each variant script (`make_vs3_n2_k32_sweep_<variant>_live_report.py`) imports it and overrides `OUT_HTML`/`ECG_ROOT`/`REPORT_TITLE`/`CASES`. `rice_mele_reference/refresh_vs3_n2_k32_reports.py --reports <keys>` regenerates a set of reports on an interval and is itself run as a long-lived Slurm job (`txt_task/vs3_n2_k32_reports_auto_refresh_*.sbatch`). Register new variants in its `REPORT_SCRIPTS` dict.
- Reports land at the repo root as `*_live_report.html` (some archived under `html/`); they compare each ECG sweep case against the 2D-grid split-step reference in `out/vs3_n2_dt0p01_T160pi_K32_sweep/grid_ref`.

## Architecture

One engine (`src/`, namespace `ecg1d`) + run configs (`pumpconfig/`) + drivers. The N=2 drivers keep `main.cpp` untouched and share their copied pump loop in `nointeraction_src/pump_common.cpp` (`pump2::run_pump` + `RunOptions`); Gaussian-specific setup is isolated in `2gaussian_src/`.

`src/` is layered bottom-up:

1. **Matrix elements**: `permutation` (N! permutation set with signs) → `pair_cache` (per ⟨bra i | ket j, permutation p⟩ cached K, K⁻¹, μ, det, overlap M_G) → `interaction_kernels` (kinetic P_Mij, delta-contact G, Gaussian H, one-body cosine kernels) → `hamiltonian` (sums kernels over i,j,p into ⟨S⟩ and ⟨H⟩ functionals). **Memoization**: `PairCacheTable` precomputes the full (i,j,p) table once per RHS evaluation; `compute_rhs_dz` installs it via `PairCacheTableScope` so every downstream kernel (overlap, Hamiltonian, metric, gradient, trace) reuses it instead of rebuilding. This is the accelerator that unlocked full-cycle N=2 (the "cost warning" in `txt_task/n2_notes.md` predates it and is stale).
2. **Analytic gradients**: `derivatives` (12-case derivative tables of log-overlap and K⁻¹ w.r.t. u/A/B/R — the largest and most delicate file), `observable_derivatives`, `hamiltonian_gradient` (∂(⟨H⟩/S)/∂z per term).
3. **TDVP**: `tdvp_solver` — assembles metric C and gradient g over an `AlphaIndex` list (the OpenMP hot loops), `HamiltonianTerms` toggles kinetic/delta/gaussian + arbitrary cosine terms. `realtime_tdvp` — `compute_rhs_dz` solves (C + λ_C·I)dz = −i·g via SVD with relative `rcond` truncation; `lambda_C` and `rcond` come from `PumpConfig` (defaults 1e-8 / 1e-4) so drivers can sweep them. Per-step diagnostics (`RealtimeStepResult`): `raw_cond` is deprecated (kept as NaN for CSV compatibility); `sv_max` is the largest singular value of the solve matrix — the active cutoff is `rcond·sv_max`. Euler/RK4/RK45 integrators plus the time-dependent-H RK4 (`realtime_tdvp_step_rk4_time_dependent`) that the pump drivers use.
4. **I/O**: `snapshot_io` (in-memory-buffered snapshots.csv; legacy N=1 schema, extended N>=2 schema), `trace_io` (observables + solver-health trace incl. `sv_max`; N>=2 trace appends `V_gauss`), `run_report` (config.txt / summary.txt / progress writers / dt-schedule helpers `time_step_at`/`next_time_step`), `progress_output_n{1,2}`, `csv_utils`.

**Manual fine-dt windows** (gap-adaptive time step): `PumpConfig.manual_fine_dt_{enabled,factor,windows_s,schedule_s}` shrink dt by `factor` while s = t/T is inside the given windows (the campaign uses 0.05× in s∈[0.2,0.3]∪[0.7,0.8], where the band gap is smallest); `schedule_s` allows per-window factors ("staged"). Wired into the N=2 loop via `next_time_step`; the N=1 `main.cpp` loop is still fixed-dt. `estimate_time_steps` accounts for the windows when planning wall time.

**N=2 `RunOptions`** (`nointeraction_src/pump_common.hpp`): `base_terms` (kinetic_only for free; gaussian driver adds `.gaussian`), `trace_vint` (sample ⟨V_gauss⟩), `evolve_A` (include A in the variational parameters — off by default, preserving frozen-A), `config_appendix`.

Driver flow (`main.cpp` and `pump2::run_pump` alike): load gap-adaptive phase schedule CSV and initial basis CSV → fixed/fine-dt RK4 TDVP loop with per-step norm re-enforcement and finiteness guards → write trace/snapshots/basis/progress/summary into a parameter-tagged output dir.

Key convention: TDVP parameters are addressed by `AlphaIndex{a1,..}` with a1 = 1:u, 2:B, 3:R, 4:A. A is frozen by default (for N=1, A and B are both width-like directions and together make C ill-conditioned), so param_dim = K·(1+2N); N=2 drivers can opt in via `evolve_A`. The `Aseed0p05` basis variant instead bakes a fixed A = 0.05·[[1,−1],[−1,1]] into the 16 path-pad terms of `initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32_Aseed0p05.csv` (generator: `rice_mele_reference/Vs3Vl3_3_3/make_vs3vl3_initial_basis_n2_pathpad_Aseed.py`) while keeping A frozen during evolution.

## Units and conventions (gotchas)

- Code units ħ = 1, m = 1 (`src/physical_constants.hpp`); `recoil_energy_for_lattice` uses k = 2π/lattice_a = π/d_s, i.e. the **short-lattice** recoil E_r,s (= π²/32 ≈ 0.3084 for a = 8) — the same energy unit as the Rice–Mele reference data (which sets d_s = 1, E_r,s = 1). Time conversion: T_code = (32/π²) · T[ħ/E_r,s] ≈ 3.24 · T.
- **Factor-of-2 depth trap**: "Vs=Vl=3" labels two different potentials — the full-depth "grid" convention vs the half-depth "main.cpp" convention that the ECG binary actually integrates. Read `rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_Vs3Vl3_convention_note.md` before comparing any ECG result to reference curves or regenerating schedules.
- Phase schedule CSV: columns `s,phi` with s = t/T strictly increasing and spanning exactly [0,1]; linearly interpolated; `#` comments allowed.
- Basis CSV format (`initial_state/`, `basis_initial/final.csv`): per Gaussian, one `re,im` line each for u, then A (N² entries), B (N²), R (N), then one integer name line.

## Outputs and how a run is judged

Each run writes `out/<out_root>/a<a>_K<K>_tmax<T>_VsER<vs>_VlER<vl>/` containing `config.txt`, `progress.csv`, `trace.csv`, `snapshots.csv`, `basis_initial.csv`, `basis_final.csv`, `summary.txt`. Sweep drivers nest per-case roots: `out/<sweep_name>/{free,gauss_sigma<tag>}/a8p000_K32_tmax502p655_.../`. Success criteria for a pump run: Δpolarization_cell ≈ −1 per cycle, small norm/energy drift, and healthy solver diagnostics — max actual-solve condition number, `sv_max`, effective rank vs param_dim, discarded-rhs fraction, min Re(B) and min Re(A+B) staying positive. These are printed per sampled step and aggregated in summary.txt.

Reference/analysis material: `rice_mele_reference/Vs3Vl3_3_3/` (phase schedules, Rice–Mele band-model benchmark data + README, N=2 path-pad basis generators), `rice_mele_reference/lohes_experience/` (quantitative reference for the Lohse 2016 paper lattice: J₁/J₂/Δ, gaps, adiabatic times, gap-adaptive schedules in both conventions + HTML report, N=1/N=2 basis generators, N=2 split-step reference), `successful_case/Vs3Vl3_3_3/` (Jupyter density-movie notebook and HTML report over the good N=1 run's snapshots.csv), `reference/` (papers), `vs3_n2_k32_failure_analysis_report.html` (root-level post-mortem of failed sweep cases).

## Notes

- Project documentation and notebook narration mix Chinese and English.
- `txt_task/` = task drivers (.cpp), sbatch scripts, and working notes (`n2_notes.md` — its full-N=2 cost warning is stale, see Architecture; `rk45_adaptive_dt_proposal.md` — adaptive-dt design proposal, not implemented).
