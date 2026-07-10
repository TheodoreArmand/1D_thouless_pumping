#pragma once

// Shared pump-driver plumbing for the N>=2 executables
// (main_2ninteraction.cpp and main_2gaussian.cpp).
//
// This is a faithful port of the anonymous-namespace helpers and the RK4
// evolution loop of main.cpp (which stays untouched and owns its private
// copies; see the provenance note in pump_common.cpp). Everything here is
// N-general: the same code drives the N=1 mini validation configs and the
// N=2 runs.
//
// Statistics note: the engine's permutation signs are hardcoded +1
// (src/permutation.cpp), i.e. SYMMETRIC (bosonic) states — correct for the
// 87Rb Lohse system these drivers target.

#include "basis_params.hpp"
#include "pumpconfig/pump_config.hpp"
#include "tdvp_solver.hpp"

#include <string>
#include <vector>

namespace pump2 {

struct RunOptions {
    // Base Hamiltonian terms; cosine terms are filled per time step from the
    // pump schedule. kinetic_only() for the free driver; the gaussian driver
    // additionally sets .gaussian = true.
    ecg1d::HamiltonianTerms base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    // Sample <V_gauss>(t) into n2_trace.csv (gaussian driver).
    bool trace_vint = false;
    // Include A upper-triangular entries in the TDVP variational parameters.
    // Kept off by default so the free/no-interaction driver preserves the
    // current frozen-A behavior.
    bool evolve_A = false;
    // When evolve_A is set, restrict A to its strict off-diagonal entries
    // (A_01 for N=2) and freeze diag A at 0 (run_pump pins it before the
    // evolution). Diagonal B + R + u already span any symmetric quadratic form
    // W = A + B, so this is the SAME primitive function family without the
    // exact gauge redundancy of full A (param_dim 256 -> 192 at N=2, K=32).
    // See vs3_n2_k32_failure_analysis_report.html section 5 / experiment E2.
    bool evolve_A_offdiag_only = false;
    // R3a occupancy gate (vs3_n2_k32_success_roadmap_report.html section 6-R3 /
    // 8.2). When enabled, a term k's B/R/A coordinates enter the TDVP evolution
    // parameter list only while |u_k| exceeds the hysteresis band: a parked term
    // wakes at |u_k| > occupancy_u_on and re-parks at |u_k| < occupancy_u_off.
    // The u channel is never gated. Parking becomes an explicit, rcond-independent
    // rule instead of a truncation side-effect. Off by default: the alpha list is
    // then the static make_alpha_list output and behavior is unchanged.
    bool occupancy_gate_enabled = false;
    // Emit occupancy_gate.csv: per-sampled-step census (core/waking/parked term
    // counts, active param dim) plus the discarded-RHS fraction decomposed by
    // parameter block (u/B/R/A) x occupancy tier (core/waking/parked). Works with
    // the gate on or off; when off it records the baseline H1 fingerprint.
    bool occupancy_diag_enabled = false;
    double occupancy_u_on = 0.03;   // wake / core threshold on |u_k|
    double occupancy_u_off = 0.02;  // re-park threshold on |u_k| (<= u_on)
    // Extra "key=value\n" lines appended to config.txt after write_config_txt
    // (keeps run_report untouched and the N=1 config.txt byte-stable).
    std::string config_appendix;
};

// Phase schedule (ports of main.cpp helpers).
void load_phase_schedule(pumpconfig::PumpConfig& cfg);
double phi_at(double t, const pumpconfig::PumpConfig& cfg);

// Basis CSV I/O (N-general; same format as main.cpp).
std::vector<ecg1d::BasisParams> load_basis_csv(const std::string& path, int N);
void write_basis_csv(const std::string& path, const std::vector<ecg1d::BasisParams>& basis);

// TDVP parameter list: u (K), B diagonal (K*N), R (K*N), optionally A entries.
// With evolve_A the A block is the upper triangle including the diagonal
// (K*N*(N+1)/2); with evolve_A_offdiag_only it is the strict upper triangle
// only (K*N*(N-1)/2), freezing diag A.
std::vector<ecg1d::AlphaIndex> make_alpha_list(int N, int K, bool evolve_A = false,
                                               bool evolve_A_offdiag_only = false);

// Full pump run: load basis, evolve, write trace/snapshots/summary (and
// n2_trace.csv with pair-separation / <V_int> when N >= 2).
// Returns the process exit code (0 ok).
int run_pump(pumpconfig::PumpConfig& cfg, const RunOptions& opt);

}  // namespace pump2
