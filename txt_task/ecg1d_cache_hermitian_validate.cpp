#include "derivatives.hpp"
#include "hamiltonian.hpp"
#include "nointeraction_src/pump_common.hpp"
#include "pair_cache.hpp"
#include "physical_constants.hpp"
#include "pumpconfig/pump_config.hpp"
#include "tdvp_solver.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

using ecg1d::Cd;
using ecg1d::HamiltonianTerms;
using ecg1d::MatrixXcd;
using ecg1d::VectorXcd;

struct TimedMatrixVector {
    MatrixXcd C;
    VectorXcd g;
    double seconds = 0.0;
};

double safe_rel(double abs_err, double scale) {
    return abs_err / std::max(scale, 1.0e-300);
}

double max_abs_matrix(const MatrixXcd& m) {
    return (m.size() > 0) ? m.cwiseAbs().maxCoeff() : 0.0;
}

double max_abs_vector(const VectorXcd& v) {
    return (v.size() > 0) ? v.cwiseAbs().maxCoeff() : 0.0;
}

void normalize_basis(std::vector<ecg1d::BasisParams>& basis) {
    const double nrm = ecg1d::overlap(basis).real();
    if (!(nrm > 0.0) || !std::isfinite(nrm)) {
        throw std::runtime_error("cannot normalize basis");
    }
    const double scale = std::sqrt(1.0 / nrm);
    for (auto& b : basis) b.u *= scale;
}

MatrixXcd assemble_C_full_reference(const std::vector<ecg1d::AlphaIndex>& alpha_z_list,
                                    const std::vector<ecg1d::BasisParams>& basis) {
    const int d = static_cast<int>(alpha_z_list.size());
    const Cd S = ecg1d::overlap(basis);
    VectorXcd first_real = VectorXcd::Zero(d);
    VectorXcd first_false = VectorXcd::Zero(d);

    #pragma omp parallel for schedule(static)
    for (int a = 0; a < d; ++a) {
        const auto& alpha = alpha_z_list[a];
        first_real(a) = ecg1d::partial_z_first(alpha.a1, true, basis,
                                               alpha.a2, alpha.a3, alpha.a4);
        first_false(a) = ecg1d::partial_z_first(alpha.a1, false, basis,
                                                alpha.a2, alpha.a3, alpha.a4);
    }

    MatrixXcd C = MatrixXcd::Zero(d, d);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int a = 0; a < d; ++a) {
        for (int b = 0; b < d; ++b) {
            const auto& alpha = alpha_z_list[a];
            const auto& beta = alpha_z_list[b];
            const Cd second = ecg1d::partial_z_second(
                alpha.a1, beta.a1, basis,
                alpha.a2, alpha.a3, alpha.a4,
                beta.a2, beta.a3, beta.a4);
            C(a, b) = -first_real(a) * first_false(b) / (S * S) + second / S;
        }
    }
    return C;
}

VectorXcd solve_tdvp_dz(const MatrixXcd& C,
                        const VectorXcd& g,
                        double lambda_C,
                        double rcond) {
    const int d = static_cast<int>(g.size());
    const MatrixXcd C_raw = C.conjugate();
    const VectorXcd rhs = Cd(0.0, -1.0) * g;
    MatrixXcd M = C_raw;
    if (lambda_C > 0.0) {
        M += lambda_C * MatrixXcd::Identity(d, d);
    }

    Eigen::SelfAdjointEigenSolver<MatrixXcd> eig(M);
    const Eigen::VectorXd ev = eig.eigenvalues();
    const MatrixXcd V = eig.eigenvectors();
    const double lambda_max = ev(ev.size() - 1);
    const double thr = rcond * lambda_max;
    const VectorXcd Vtb = V.adjoint() * rhs;
    VectorXcd coeff(ev.size());
    for (int i = 0; i < ev.size(); ++i) {
        coeff(i) = (ev(i) > thr) ? (Vtb(i) / ev(i)) : Cd(0.0, 0.0);
    }
    return V * coeff;
}

TimedMatrixVector assemble_C_grad_uncached(const std::vector<ecg1d::AlphaIndex>& alpha,
                                           const std::vector<ecg1d::BasisParams>& basis,
                                           const HamiltonianTerms& terms,
                                           int repeats) {
    TimedMatrixVector out;
    const auto start = std::chrono::steady_clock::now();
    for (int r = 0; r < repeats; ++r) {
        out.C = ecg1d::assemble_C(alpha, basis);
        out.g = ecg1d::assemble_grad(alpha, basis, terms);
    }
    const auto end = std::chrono::steady_clock::now();
    out.seconds = std::chrono::duration<double>(end - start).count() / repeats;
    return out;
}

TimedMatrixVector assemble_C_grad_cached(const std::vector<ecg1d::AlphaIndex>& alpha,
                                         const std::vector<ecg1d::BasisParams>& basis,
                                         const HamiltonianTerms& terms,
                                         int repeats) {
    TimedMatrixVector out;
    const auto start = std::chrono::steady_clock::now();
    for (int r = 0; r < repeats; ++r) {
        ecg1d::PairCacheTable pair_cache_table(basis);
        ecg1d::PairCacheTableScope pair_cache_scope(pair_cache_table);
        out.C = ecg1d::assemble_C(alpha, basis);
        out.g = ecg1d::assemble_grad(alpha, basis, terms);
    }
    const auto end = std::chrono::steady_clock::now();
    out.seconds = std::chrono::duration<double>(end - start).count() / repeats;
    return out;
}

HamiltonianTerms make_terms(const pumpconfig::PumpConfig& cfg,
                            const std::string& mode,
                            double s) {
    const double phi = pump2::phi_at(s * cfg.pump_period, cfg);
    HamiltonianTerms terms = HamiltonianTerms::kinetic_only();
    terms.cosine_terms.push_back({-cfg.Vs, 4.0 * pumpconfig::pi / cfg.lattice_a, 0.0});
    terms.cosine_terms.push_back({-cfg.Vl, 2.0 * pumpconfig::pi / cfg.lattice_a, phi});
    if (mode == "gauss") {
        terms.gaussian = true;
        ecg1d::g_gauss = 0.3 * cfg.recoil_energy;
        ecg1d::sigma_gauss = 1.0;
    }
    return terms;
}

void print_metric(const std::string& name, double abs_err, double scale) {
    std::cout << std::left << std::setw(34) << name
              << " abs=" << std::scientific << std::setprecision(6) << abs_err
              << " rel=" << safe_rel(abs_err, scale) << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = (argc >= 2) ? argv[1] : "gauss";
        if (mode != "free" && mode != "gauss") {
            throw std::runtime_error("usage: ecg1d_cache_hermitian_validate [free|gauss] [K] [repeats] [s]");
        }
        const int K = (argc >= 3) ? std::stoi(argv[2]) : 8;
        const int repeats = std::max(1, (argc >= 4) ? std::stoi(argv[3]) : 1);
        const double s = (argc >= 5) ? std::stod(argv[4]) : 0.25;
        if (K <= 0 || K > 32) throw std::runtime_error("K must be in 1..32");
        if (!(s >= 0.0 && s <= 1.0) || !std::isfinite(s)) {
            throw std::runtime_error("s must be finite and in [0,1]");
        }

        pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
        cfg.N = 2;
        cfg.K = K;
        cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv";
        pump2::load_phase_schedule(cfg);

        std::vector<ecg1d::BasisParams> basis = pump2::load_basis_csv(cfg.basis_from, cfg.N);
        if (static_cast<int>(basis.size()) < K) {
            throw std::runtime_error("basis file has fewer functions than requested K");
        }
        basis.resize(static_cast<std::size_t>(K));
        for (int i = 0; i < K; ++i) basis[static_cast<std::size_t>(i)].name = i;
        normalize_basis(basis);

        const bool evolve_A = (mode == "gauss");
        const std::vector<ecg1d::AlphaIndex> alpha = pump2::make_alpha_list(cfg.N, cfg.K, evolve_A);
        const HamiltonianTerms terms = make_terms(cfg, mode, s);

        std::cout << std::setprecision(10)
                  << "mode=" << mode
                  << " K=" << K
                  << " param_dim=" << alpha.size()
                  << " repeats=" << repeats
                  << " s=" << s
                  << " OMP_threads=" << omp_get_max_threads()
                  << "\n";

        const TimedMatrixVector uncached = assemble_C_grad_uncached(alpha, basis, terms, repeats);
        const TimedMatrixVector cached = assemble_C_grad_cached(alpha, basis, terms, repeats);

        MatrixXcd C_half_timed;
        double half_C_seconds = 0.0;
        {
            const auto start = std::chrono::steady_clock::now();
            for (int r = 0; r < repeats; ++r) {
                ecg1d::PairCacheTable pair_cache_table(basis);
                ecg1d::PairCacheTableScope pair_cache_scope(pair_cache_table);
                C_half_timed = ecg1d::assemble_C(alpha, basis);
            }
            const auto end = std::chrono::steady_clock::now();
            half_C_seconds = std::chrono::duration<double>(end - start).count() / repeats;
        }

        MatrixXcd C_full;
        double full_seconds = 0.0;
        {
            const auto start = std::chrono::steady_clock::now();
            ecg1d::PairCacheTable pair_cache_table(basis);
            ecg1d::PairCacheTableScope pair_cache_scope(pair_cache_table);
            C_full = assemble_C_full_reference(alpha, basis);
            const auto end = std::chrono::steady_clock::now();
            full_seconds = std::chrono::duration<double>(end - start).count();
        }

        const VectorXcd dz_cached = solve_tdvp_dz(cached.C, cached.g, cfg.lambda_C, cfg.rcond);
        const VectorXcd dz_uncached = solve_tdvp_dz(uncached.C, uncached.g, cfg.lambda_C, cfg.rcond);
        const VectorXcd dz_full = solve_tdvp_dz(C_full, cached.g, cfg.lambda_C, cfg.rcond);

        const MatrixXcd dC_cache = cached.C - uncached.C;
        const VectorXcd dg_cache = cached.g - uncached.g;
        const MatrixXcd dC_full = cached.C - C_full;
        const MatrixXcd herm_res = cached.C - cached.C.adjoint();
        const VectorXcd ddz_cache = dz_cached - dz_uncached;
        const VectorXcd ddz_full = dz_cached - dz_full;

        std::cout << std::scientific << std::setprecision(6)
                  << "seconds_per_C_plus_grad_uncached=" << uncached.seconds << "\n"
                  << "seconds_per_C_plus_grad_pair_cached=" << cached.seconds << "\n"
                  << "speedup_cached_vs_uncached=" << (uncached.seconds / cached.seconds) << "\n"
                  << "seconds_half_C_pair_cached=" << half_C_seconds << "\n"
                  << "seconds_full_C_reference_once=" << full_seconds << "\n";
        if (half_C_seconds > 0.0) {
            std::cout << "speedup_half_C_vs_full_C=" << (full_seconds / half_C_seconds) << "\n";
        }

        print_metric("C cached - uncached", max_abs_matrix(dC_cache), max_abs_matrix(uncached.C));
        print_metric("C cached - timed half", max_abs_matrix(cached.C - C_half_timed),
                     max_abs_matrix(cached.C));
        print_metric("g cached - uncached", max_abs_vector(dg_cache), max_abs_vector(uncached.g));
        print_metric("C half - explicit full", max_abs_matrix(dC_full), max_abs_matrix(C_full));
        print_metric("C hermiticity residual", max_abs_matrix(herm_res), max_abs_matrix(cached.C));
        print_metric("dz cached - uncached", max_abs_vector(ddz_cache), max_abs_vector(dz_uncached));
        print_metric("dz half - fullC", max_abs_vector(ddz_full), max_abs_vector(dz_full));

        const double pass_tol = 1.0e-9;
        const bool pass =
            safe_rel(max_abs_matrix(dC_cache), max_abs_matrix(uncached.C)) < pass_tol &&
            safe_rel(max_abs_vector(dg_cache), max_abs_vector(uncached.g)) < pass_tol &&
            safe_rel(max_abs_matrix(dC_full), max_abs_matrix(C_full)) < pass_tol &&
            safe_rel(max_abs_matrix(herm_res), max_abs_matrix(cached.C)) < pass_tol &&
            safe_rel(max_abs_vector(ddz_cache), max_abs_vector(dz_uncached)) < pass_tol &&
            safe_rel(max_abs_vector(ddz_full), max_abs_vector(dz_full)) < pass_tol;

        std::cout << "validation=" << (pass ? "PASS" : "FAIL") << "\n";
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "FATAL cache/hermitian validation: " << e.what() << "\n";
        return 2;
    }
}
