#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <functional>
#include <vector>
#include <limits>

namespace ecg1d {

// Real-time TDVP step for Schrodinger equation i*d|psi>/dt = H|psi>:
//   i * C_{alpha,beta} * z_dot_beta = g_alpha
//   -> dz = -i * C^{-1} * g    (C here is the tdvp metric with P_perp)
// vs imaginary-time version  (C * z_dot = -g  ->  dz = -C^{-1} * g)
//
// `dt` is a real positive step; integrator options (order 1 Euler or order 4 RK4).
//
// Returns dz as VectorXcd of size = alpha_z_list.size(); caller applies to basis.
enum class RtIntegrator { Euler, RK4, RK45 };

// Optional per-parameter grouping used to decompose the rcond-discarded RHS by
// arbitrary caller-defined groups (e.g. parameter block u/B/R/A x occupancy tier
// core/waking/parked). group_of_param[p] is the group id in [0, n_groups) of the
// p-th entry of alpha_z_list. When supplied to the k1 solve the step result gets
// discarded_sq_by_group[g] = ||P_disc rhs restricted to group g||^2 and
// total_sq_by_group[g] = ||rhs restricted to group g||^2, so per-group
// sqrt(disc/total) is that group's truncated fraction of the gradient drive.
// This is the H1 discriminator from vs3_n2_k32_success_roadmap_report.html section 8.2.
struct RhsGrouping {
    std::vector<int> group_of_param;  // size must equal alpha_z_list.size()
    int n_groups = 0;
};

struct RealtimeStepResult {
    std::vector<BasisParams> basis;   // updated basis
    double used_dt;
    double cond_C;                    // == actual_solve_cond (matrix actually SVD'd)
    int    effective_rank;            // rank of the actual solve matrix under rcond
    // Diagnostic: smallest 3 singular values of the actual solve matrix at step
    // entry (k1 stage). sv_small[0] = smallest, [1] = second, [2] = third.
    // NaN if fewer than that many singular values exist.
    double sv_small[3] = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};
    // Diagnostic: norm of dz returned by the full step (combined RK4 dz).
    double dz_norm = 0.0;

    // Extended k1-stage linear-algebra diagnostics (see compute_rhs_dz).
    // raw_cond:        deprecated; kept as NaN for legacy CSV compatibility.
    // actual_solve_cond: cond of the matrix actually fed to the SVD solve
    //                  (= C + lambda_C I). Equal to cond_C; kept under an
    //                  explicit name for clarity.
    // sv_max:          largest eigenvalue/singular value of that Hermitian solve
    //                  matrix. The active relative cutoff is rcond * sv_max.
    // relative_raw_residual: ||C_raw dz - rhs|| / ||rhs|| for the k1 solution,
    //                  i.e. how well dz solves the original TDVP equation.
    // discarded_rhs_fraction: fraction of the rhs (in the solve/U basis) removed
    //                  by the rcond truncation.
    // metric_norm:     sqrt(real(dz^dagger C_raw dz)) for the k1 dz.
    double raw_cond                  = std::numeric_limits<double>::quiet_NaN();
    double actual_solve_cond         = std::numeric_limits<double>::quiet_NaN();
    double sv_max                    = std::numeric_limits<double>::quiet_NaN();
    double relative_raw_residual     = std::numeric_limits<double>::quiet_NaN();
    double discarded_rhs_fraction    = std::numeric_limits<double>::quiet_NaN();
    double metric_norm               = std::numeric_limits<double>::quiet_NaN();

    // Grouped decomposition of the discarded/total RHS (filled only when a
    // RhsGrouping is passed to the step). Both are size n_groups; entry g holds
    // the squared RHS norm in group g that was discarded / present in total.
    // Empty when no grouping was requested.
    std::vector<double> discarded_sq_by_group;
    std::vector<double> total_sq_by_group;
};

// Apply one real-time TDVP step. Does NOT adaptively rescale dt.
RealtimeStepResult realtime_tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                                      const std::vector<BasisParams>& basis,
                                      double dt,
                                      const HamiltonianTerms& terms,
                                      const SolverConfig& config,
                                      RtIntegrator integrator = RtIntegrator::RK4);

// Classical RK4 for explicitly time-dependent Hamiltonians:
//   k1 uses H(t), k2/k3 use H(t + dt/2), k4 uses H(t + dt).
// This keeps the old frozen-H interface intact while allowing drivers to use
// true stage-time Hamiltonians.
// When `grouping` is non-null it is applied to the k1 solve to fill the
// grouped discarded/total RHS breakdown in the returned RealtimeStepResult.
RealtimeStepResult realtime_tdvp_step_rk4_time_dependent(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    double t,
    double dt,
    const std::function<HamiltonianTerms(double)>& terms_at,
    const SolverConfig& config,
    const RhsGrouping* grouping = nullptr);

struct RealtimeTrace {
    std::vector<double> t;
    std::vector<double> E;        // <H> / <psi|psi>  (real, conserved)
    std::vector<double> norm;     // <psi|psi> (real, conserved)
    std::vector<double> x_mean;   // <X> = sum_a <x_a>  (any N; normalized)
    std::vector<double> p_mean;   // <P> = sum_a <p_a>  (any N; normalized)
    std::vector<double> x2;       // <x^2>                          normalized
    std::vector<double> p2;       // <p^2>                          normalized
    // Raw (un-normalized) moments: <psi|A|psi> without dividing by <psi|psi>.
    // Useful as numerical witness of the Q9 structural conservation laws under
    // SVD-truncated TDVP (norm leaks but raw <psi|H|psi> stays flat to round-off).
    std::vector<double> x_mean_raw; // <psi|X|psi>   (any N)
    std::vector<double> p_mean_raw; // <psi|P|psi>   (any N)
    std::vector<double> x2_raw;     // <psi|x^2|psi>
    std::vector<double> p2_raw;     // <psi|p^2|psi>
    std::vector<Cd>     overlap0; // <psi(0)|psi(t)> — fidelity |<psi(0)|psi(t)>|^2/(n0*n(t))
};

struct RealtimeEvolutionConfig {
    double dt          = 1e-3;
    RtIntegrator integrator = RtIntegrator::RK4;
    int    sample_every = 10;        // record trace every N steps
    bool   verbose      = true;
    int    print_every  = 100;

    // RK45 (Dormand-Prince 5(4)) adaptive stepping — ignored unless integrator=RK45
    //   rk45_abs_tol: absolute tolerance
    //   rk45_rel_tol: relative tolerance
    //   rk45_dt_min/max: clamp on step size
    //   The initial dt is taken from rt_cfg.dt.
    double rk45_abs_tol = 1e-6;
    double rk45_rel_tol = 1e-6;
    double rk45_dt_min  = 1e-8;
    double rk45_dt_max  = 1.0;

    // If true, after each step rescale u so that <psi|psi> == initial_norm
    // (where initial_norm is computed at t=0 from basis_init). This is a
    // pragmatic post-correction that enforces physical unitarity when the raw
    // TDVP+SVD step drifts norm due to rcond truncation or Tikhonov λ not
    // being exactly consistent with ψ̇ ⊥ ψ. E conservation is NOT enforced,
    // only norm.
    bool enforce_norm = false;
};

struct RealtimeEvolutionResult {
    std::vector<BasisParams> basis_final;
    RealtimeTrace            trace;
    int    n_steps;
};

// Evolve `basis_init` from t=0 to t=T_total using real-time TDVP.
// At every sample_every'th step, record observables into RealtimeTrace.
RealtimeEvolutionResult realtime_tdvp_evolution(
    const std::vector<AlphaIndex>& alpha_z_list,
    std::vector<BasisParams> basis_init,
    double T_total,
    const HamiltonianTerms& terms,
    const SolverConfig& solver_cfg = SolverConfig{},
    const RealtimeEvolutionConfig& rt_cfg = {});

} // namespace ecg1d
