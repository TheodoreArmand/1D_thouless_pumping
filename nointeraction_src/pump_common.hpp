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

// TDVP parameter list: u (K), B diagonal (K*N), R (K*N), optionally A upper
// triangle (K*N*(N+1)/2).
std::vector<ecg1d::AlphaIndex> make_alpha_list(int N, int K, bool evolve_A = false);

// Full pump run: load basis, evolve, write trace/snapshots/summary (and
// n2_trace.csv with pair-separation / <V_int> when N >= 2).
// Returns the process exit code (0 ok).
int run_pump(pumpconfig::PumpConfig& cfg, const RunOptions& opt);

}  // namespace pump2
