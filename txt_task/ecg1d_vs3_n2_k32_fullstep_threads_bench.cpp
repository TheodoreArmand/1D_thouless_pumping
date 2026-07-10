#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include "hamiltonian.hpp"
#include "pumpconfig/pump_config.hpp"
#include "realtime_tdvp.hpp"

#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <omp.h>

namespace {

pumpconfig::PumpConfig make_cfg() {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
    cfg.N = 2;
    cfg.K = 32;
    cfg.dt = 5.0e-4;
    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv";
    cfg.g_gauss_over_Er = 0.3;
    cfg.sigma_gauss = 1.0;
    return cfg;
}

void normalize_basis(std::vector<ecg1d::BasisParams>& basis) {
    const double nrm = ecg1d::overlap(basis).real();
    if (!(nrm > 0.0) || !std::isfinite(nrm)) {
        throw std::runtime_error("cannot normalize basis");
    }
    const double scale = std::sqrt(1.0 / nrm);
    for (auto& b : basis) b.u *= scale;
}

ecg1d::HamiltonianTerms terms_at_time(double t,
                                      const pumpconfig::PumpConfig& cfg,
                                      const ecg1d::HamiltonianTerms& base_terms) {
    ecg1d::HamiltonianTerms terms = base_terms;
    const double phi = pump2::phi_at(t, cfg);
    terms.cosine_terms.push_back({-cfg.Vs, 4.0 * pumpconfig::pi / cfg.lattice_a, 0.0});
    terms.cosine_terms.push_back({-cfg.Vl, 2.0 * pumpconfig::pi / cfg.lattice_a, phi});
    return terms;
}

int run_bench(const std::string& mode, int steps, double s0) {
    pumpconfig::PumpConfig cfg = make_cfg();
    pump2::load_phase_schedule(cfg);

    std::vector<ecg1d::BasisParams> basis = pump2::load_basis_csv(cfg.basis_from, cfg.N);
    if (static_cast<int>(basis.size()) != cfg.K) {
        throw std::runtime_error("basis size mismatch");
    }
    normalize_basis(basis);

    ecg1d::HamiltonianTerms base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    bool evolve_A = false;
    bool evolve_A_offdiag_only = false;
    if (mode == "gauss") {
        pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
        base_terms = opt.base_terms;
        evolve_A = opt.evolve_A;
        evolve_A_offdiag_only = opt.evolve_A_offdiag_only;
    } else if (mode != "free") {
        throw std::runtime_error("mode must be free or gauss");
    }

    std::vector<ecg1d::AlphaIndex> alpha =
        pump2::make_alpha_list(cfg.N, cfg.K, evolve_A, evolve_A_offdiag_only);
    ecg1d::SolverConfig solver_cfg;
    solver_cfg.lambda_C = cfg.lambda_C;
    solver_cfg.rcond = cfg.rcond;

    double t = s0 * cfg.pump_period;
    const auto wall0 = std::chrono::steady_clock::now();
    for (int step = 0; step < steps; ++step) {
        auto terms_at = [&](double stage_t) {
            return terms_at_time(stage_t, cfg, base_terms);
        };
        ecg1d::RealtimeStepResult r =
            ecg1d::realtime_tdvp_step_rk4_time_dependent(alpha, basis, t, cfg.dt, terms_at, solver_cfg);
        basis = std::move(r.basis);
        t += cfg.dt;
        const double norm_now = ecg1d::overlap(basis).real();
        if (norm_now > 1e-15 && std::isfinite(norm_now)) {
            const double scale = std::sqrt(1.0 / norm_now);
            for (auto& bp : basis) bp.u *= scale;
        }
    }
    const auto wall1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(wall1 - wall0).count();
    const double seconds_per_step = seconds / static_cast<double>(steps);
    const double final_norm = ecg1d::overlap(basis).real();

    std::cout << std::setprecision(10)
              << "mode=" << mode << "\n"
              << "K=" << cfg.K << "\n"
              << "param_dim=" << alpha.size() << "\n"
              << "steps=" << steps << "\n"
              << "s0=" << s0 << "\n"
              << "dt=" << cfg.dt << "\n"
              << "OMP_threads=" << omp_get_max_threads() << "\n"
              << "wall_seconds=" << std::scientific << std::setprecision(6) << seconds << "\n"
              << "seconds_per_full_rk4_step=" << seconds_per_step << "\n"
              << "final_norm=" << std::setprecision(12) << final_norm << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = (argc >= 2) ? argv[1] : "gauss";
        const int steps = (argc >= 3) ? std::stoi(argv[2]) : 20;
        const double s0 = (argc >= 4) ? std::stod(argv[3]) : 0.25;
        if (steps <= 0) throw std::runtime_error("steps must be positive");
        if (!(s0 >= 0.0 && s0 <= 1.0) || !std::isfinite(s0)) {
            throw std::runtime_error("s0 must be finite and in [0,1]");
        }
        return run_bench(mode, steps, s0);
    } catch (const std::exception& e) {
        std::cerr << "FATAL K32 full-step thread bench: " << e.what() << "\n";
        return 2;
    }
}
