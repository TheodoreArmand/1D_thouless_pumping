#include "nointeraction_src/pump_common.hpp"

#include "csv_utils.hpp"
#include "hamiltonian.hpp"
#include "interaction_kernels.hpp"
#include "pair_cache.hpp"
#include "permutation.hpp"
#include "physical_constants.hpp"
#include "realtime_tdvp.hpp"
#include "run_report.hpp"
#include "snapshot_io.hpp"
#include "trace_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace pump2 {
namespace {

using ecg1d::BasisParams;
using ecg1d::Cd;
using ecg1d::CosineTerm;
using ecg1d::HamiltonianTerms;
using ecg1d::MatrixXcd;
using ecg1d::VectorXcd;
using pumpconfig::PumpConfig;
using pumpconfig::pi;

void validate_config(const PumpConfig& cfg) {
    if (!(cfg.N >= 1) || !(cfg.K >= 1)) {
        throw std::runtime_error("N and K must be positive");
    }
    if (!(cfg.lattice_a > 0.0) || !std::isfinite(cfg.lattice_a)) {
        throw std::runtime_error("lattice_a must be positive and finite");
    }
    if (!(cfg.recoil_energy > 0.0) || !std::isfinite(cfg.recoil_energy)) {
        throw std::runtime_error("recoil_energy must be positive and finite");
    }
    if (!(cfg.Vs_depth_over_Er > 0.0) || !(cfg.Vl_depth_over_Er > 0.0)
        || !std::isfinite(cfg.Vs_depth_over_Er) || !std::isfinite(cfg.Vl_depth_over_Er)) {
        throw std::runtime_error("lattice depths over Er must be positive and finite");
    }
    if (!(cfg.total_time > 0.0) || !(cfg.pump_period > 0.0)) {
        throw std::runtime_error("total_time and pump_period must be positive");
    }
    if (!(cfg.dt > 0.0) || !std::isfinite(cfg.dt)) {
        throw std::runtime_error("fixed dt must be positive and finite");
    }
    if (cfg.phase_s.empty()) {
        throw std::runtime_error("phase schedule CSV was not loaded");
    }
}

std::vector<CosineTerm> pump_terms_at(double t, const PumpConfig& cfg) {
    const double phi = phi_at(t, cfg);
    std::vector<CosineTerm> terms;
    terms.push_back({-cfg.Vs, 4.0 * pi / cfg.lattice_a, 0.0});
    terms.push_back({-cfg.Vl, 2.0 * pi / cfg.lattice_a, +phi});
    return terms;
}

HamiltonianTerms hamiltonian_at(double t, const PumpConfig& cfg, const RunOptions& opt) {
    HamiltonianTerms terms = opt.base_terms;
    terms.cosine_terms = pump_terms_at(t, cfg);
    return terms;
}

void write_cd(std::ofstream& f, Cd z) {
    f << std::setprecision(17) << z.real() << "," << z.imag() << "\n";
}

double normalize_basis(std::vector<BasisParams>& basis, double target_norm = 1.0) {
    const double nrm = ecg1d::overlap(basis).real();
    if (!(nrm > 0.0) || !std::isfinite(nrm)) {
        throw std::runtime_error("cannot normalize basis with nonpositive/nonfinite norm");
    }
    const double scale = std::sqrt(target_norm / nrm);
    for (auto& b : basis) b.u *= scale;
    return nrm;
}

bool basis_all_finite(const std::vector<BasisParams>& basis) {
    for (const auto& b : basis) {
        if (!std::isfinite(b.u.real()) || !std::isfinite(b.u.imag())) return false;
        for (int i = 0; i < b.N(); ++i) {
            if (!std::isfinite(b.R(i).real()) || !std::isfinite(b.R(i).imag())) return false;
            for (int j = 0; j < b.N(); ++j) {
                if (!std::isfinite(b.A(i, j).real()) || !std::isfinite(b.A(i, j).imag())) return false;
                if (!std::isfinite(b.B(i, j).real()) || !std::isfinite(b.B(i, j).imag())) return false;
            }
        }
    }
    return true;
}

void append_config_txt(const std::string& path, const std::string& text) {
    if (text.empty()) return;
    std::ofstream f(path, std::ios::app);
    if (!f.is_open()) throw std::runtime_error("cannot append " + path);
    f << text;
    if (text.back() != '\n') f << "\n";
}

struct N2Diagnostics {
    std::vector<double> t;
    std::vector<double> phi;
    std::vector<double> r12_sq;
    std::vector<double> V_gauss;
};

void sample_n2_diagnostics(const std::vector<BasisParams>& basis,
                           double t,
                           double phi,
                           const HamiltonianTerms& terms,
                           N2Diagnostics& diag) {
    if (basis.empty() || basis[0].N() != 2) return;

    const Cd S = ecg1d::overlap(basis);
    const int K = static_cast<int>(basis.size());
    const ecg1d::PermutationSet perms = ecg1d::PermutationSet::generate(2);
    Cd amp_r12_sq(0.0, 0.0);
    for (int i = 0; i < K; ++i) {
        const Cd conj_ui = std::conj(basis[i].u);
        for (int j = 0; j < K; ++j) {
            Cd term(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                ecg1d::PairCache c = ecg1d::PairCache::build(basis[i], basis[j], perms.matrices[p]);
                const Cd h = c.K_inv(0, 0) + c.K_inv(1, 1) - 2.0 * c.K_inv(0, 1);
                const Cd d = c.mu(0) - c.mu(1);
                const Cd kernel = 0.5 * h + d * d;
                term += static_cast<double>(perms.signs[p]) * c.M_G * kernel;
            }
            amp_r12_sq += conj_ui * basis[j].u * term;
        }
    }

    const double v_gauss = terms.gaussian
        ? (ecg1d::Gaussian_interaction_functional(basis) / S).real()
        : 0.0;
    diag.t.push_back(t);
    diag.phi.push_back(phi);
    diag.r12_sq.push_back((amp_r12_sq / S).real());
    diag.V_gauss.push_back(v_gauss);
}

void write_n2_diagnostics(const std::string& path, const N2Diagnostics& diag) {
    if (diag.t.empty()) return;
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    f << "t,phi,r12_sq,V_gauss\n";
    f << std::setprecision(17);
    for (size_t i = 0; i < diag.t.size(); ++i) {
        f << diag.t[i] << ","
          << diag.phi[i] << ","
          << diag.r12_sq[i] << ","
          << diag.V_gauss[i] << "\n";
    }
}

}  // namespace

void load_phase_schedule(PumpConfig& cfg) {
    if (cfg.phase_schedule_csv.empty()) {
        throw std::runtime_error("phase_schedule_csv must point to a phase schedule CSV");
    }

    cfg.phase_s.clear();
    cfg.phase_phi.clear();
    std::ifstream f(cfg.phase_schedule_csv);
    if (!f.is_open()) throw std::runtime_error("cannot open phase schedule CSV: " + cfg.phase_schedule_csv);

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("s,", 0) == 0) continue;
        const std::vector<std::string> cols = ecg1d::split_simple_csv_line(line);
        if (cols.size() < 2) continue;
        const double s = std::stod(cols[0]);
        const double phi = std::stod(cols[1]);
        if (std::isfinite(s) && std::isfinite(phi)) {
            cfg.phase_s.push_back(s);
            cfg.phase_phi.push_back(phi);
        }
    }

    if (cfg.phase_s.size() < 2) {
        throw std::runtime_error("phase schedule CSV has fewer than two rows: " + cfg.phase_schedule_csv);
    }
    for (size_t i = 1; i < cfg.phase_s.size(); ++i) {
        if (!(cfg.phase_s[i] > cfg.phase_s[i - 1])) {
            throw std::runtime_error("phase schedule s column must be strictly increasing");
        }
    }
    if (std::abs(cfg.phase_s.front()) > 1e-10 || std::abs(cfg.phase_s.back() - 1.0) > 1e-10) {
        throw std::runtime_error("phase schedule s column must run from 0 to 1");
    }
}

double phi_at(double t, const PumpConfig& cfg) {
    const double fraction = t / cfg.pump_period;
    if (fraction <= 0.0) return cfg.phase_phi.front();
    if (fraction >= 1.0) return cfg.phase_phi.back();
    auto it = std::lower_bound(cfg.phase_s.begin(), cfg.phase_s.end(), fraction);
    const size_t hi = static_cast<size_t>(std::distance(cfg.phase_s.begin(), it));
    const size_t lo = hi - 1;
    const double denom = cfg.phase_s[hi] - cfg.phase_s[lo];
    const double w = (fraction - cfg.phase_s[lo]) / denom;
    return cfg.phase_phi[lo] + w * (cfg.phase_phi[hi] - cfg.phase_phi[lo]);
}

std::vector<BasisParams> load_basis_csv(const std::string& path, int N) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);

    auto read_cd = [&]() -> Cd {
        std::string line;
        do {
            if (!std::getline(f, line)) throw std::runtime_error("short CSV: " + path);
        } while (line.empty());
        std::stringstream ss(line);
        double re = 0.0, im = 0.0;
        char comma = ',';
        ss >> re >> comma >> im;
        return Cd(re, im);
    };
    auto read_int = [&]() -> int {
        std::string line;
        if (!std::getline(f, line)) throw std::runtime_error("short CSV: " + path);
        return std::stoi(line);
    };

    std::vector<BasisParams> basis;
    while (f.peek() != EOF) {
        BasisParams b;
        b.u = read_cd();
        b.A = MatrixXcd(N, N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) b.A(i, j) = read_cd();
        b.B = MatrixXcd(N, N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) b.B(i, j) = read_cd();
        b.R = VectorXcd(N);
        for (int i = 0; i < N; ++i) b.R(i) = read_cd();
        b.name = read_int();
        basis.push_back(b);
        while (f.peek() == '\n' || f.peek() == '\r') f.get();
    }
    return basis;
}

void write_basis_csv(const std::string& path, const std::vector<BasisParams>& basis) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    for (const auto& b : basis) {
        const int N = b.N();
        write_cd(f, b.u);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) write_cd(f, b.A(i, j));
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) write_cd(f, b.B(i, j));
        for (int i = 0; i < N; ++i) write_cd(f, b.R(i));
        f << b.name << "\n";
    }
}

std::vector<ecg1d::AlphaIndex> make_alpha_list(int N, int K, bool evolve_A) {
    std::vector<ecg1d::AlphaIndex> alpha;
    for (int i = 0; i < K; ++i) alpha.push_back({1, i, 0, 0});
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < N; ++j) alpha.push_back({2, i, j, 0});
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < N; ++j) alpha.push_back({3, i, j, 0});
    if (evolve_A) {
        for (int i = 0; i < K; ++i)
            for (int m = 0; m < N; ++m)
                for (int n = m; n < N; ++n) alpha.push_back({4, i, m, n});
    }
    return alpha;
}

int run_pump(PumpConfig& cfg, const RunOptions& opt) {
    try {
        load_phase_schedule(cfg);
        validate_config(cfg);

        std::vector<BasisParams> basis = load_basis_csv(cfg.basis_from, cfg.N);
        if (static_cast<int>(basis.size()) != cfg.K) {
            throw std::runtime_error("basis CSV size does not match locked K");
        }
        normalize_basis(basis);

        std::string task_dir = cfg.out_root + "/a" + ecg1d::format_output_tag(cfg.lattice_a)
                             + "_K" + std::to_string(cfg.K)
                             + "_tmax" + ecg1d::format_output_tag(cfg.total_time);
        task_dir += "_VsER" + ecg1d::format_output_tag(cfg.Vs / cfg.recoil_energy)
                  + "_VlER" + ecg1d::format_output_tag(cfg.Vl / cfg.recoil_energy);
        std::filesystem::create_directories(task_dir);

        const std::string config_path = task_dir + "/config.txt";
        ecg1d::write_config_txt(config_path, cfg);
        append_config_txt(config_path, opt.config_appendix);
        write_basis_csv(task_dir + "/basis_initial.csv", basis);

        ecg1d::SolverConfig solver_cfg;
        solver_cfg.lambda_C = cfg.lambda_C;
        solver_cfg.rcond = cfg.rcond;

        std::vector<ecg1d::AlphaIndex> alpha = make_alpha_list(cfg.N, cfg.K, opt.evolve_A);
        const int param_dim = static_cast<int>(alpha.size());
        if (param_dim == 0) throw std::runtime_error("TDVP parameter list is empty");

        std::cout << "config=" << cfg.config_name
                  << " lambda_C=" << std::scientific << std::setprecision(2) << cfg.lambda_C
                  << " rcond=" << cfg.rcond << std::defaultfloat
                  << " | TDVP param dimension=" << param_dim << "\n";

        ecg1d::SnapshotSaver snapshots(task_dir);
        ecg1d::Trace trace;
        N2Diagnostics n2_diag;
        double max_cond_C = 0.0;

        const double dnan = std::numeric_limits<double>::quiet_NaN();
        double max_raw_cond = 0.0;
        double max_actual_solve_cond = 0.0;
        int min_actual_solve_rank = std::numeric_limits<int>::max();
        double max_relative_raw_residual = 0.0;
        double max_discarded_rhs_fraction = 0.0;
        double min_re_B_run = std::numeric_limits<double>::infinity();
        double min_re_AplusB_run = std::numeric_limits<double>::infinity();
        double last_sv_small[3] = {dnan, dnan, dnan};
        double min_accepted_dt = std::numeric_limits<double>::infinity();
        double max_accepted_dt = 0.0;

        int step = 0;
        double t = 0.0;
        const double initial_norm = ecg1d::overlap(basis).real();
        const double phi0 = phi_at(t, cfg);
        HamiltonianTerms terms0 = hamiltonian_at(t, cfg, opt);
        ecg1d::sample_observables(basis, t, phi0, terms0, cfg.lattice_a, trace);
        {
            ecg1d::WidthMonitors w0 = ecg1d::compute_width_monitors(basis);
            min_re_B_run = std::min(min_re_B_run, w0.min_re_B);
            min_re_AplusB_run = std::min(min_re_AplusB_run, w0.min_re_AplusB);
            ecg1d::append_trace_diagnostics(trace, w0, nullptr);
        }
        sample_n2_diagnostics(basis, t, phi0, terms0, n2_diag);
        snapshots.save("initial", step, t, phi0, basis);
        std::cout << "initial"
                  << " t=" << std::fixed << std::setprecision(5) << trace.t.front()
                  << " tau=" << phi0
                  << " E=" << std::setprecision(10) << trace.E_total.front()
                  << " <X>/a=" << trace.polarization_cell.front()
                  << " norm=" << trace.norm.front();
        if (!trace.V_gauss.empty()) {
            std::cout << " V_gauss=" << std::setprecision(10) << trace.V_gauss.front();
        }
        std::cout << std::defaultfloat << "\n";

        auto ensure_finite_basis = [&](const char* where, int next_step, double next_t) {
            if (!basis_all_finite(basis)) {
                std::ostringstream oss;
                oss << "nonfinite basis " << where
                    << " at step " << next_step << " t=" << next_t;
                throw std::runtime_error(oss.str());
            }
        };

        const auto evolution_wall_start = std::chrono::steady_clock::now();
        while (t < cfg.total_time - 1e-15) {
            const double dt_step = std::min(cfg.dt, cfg.total_time - t);
            auto terms_at_stage = [&](double stage_t) {
                return hamiltonian_at(stage_t, cfg, opt);
            };

            ecg1d::RealtimeStepResult r = ecg1d::realtime_tdvp_step_rk4_time_dependent(
                alpha, basis, t, dt_step, terms_at_stage, solver_cfg);
            basis = std::move(r.basis);
            ensure_finite_basis("after TDVP RK4", step + 1, t + dt_step);

            min_accepted_dt = std::min(min_accepted_dt, r.used_dt);
            max_accepted_dt = std::max(max_accepted_dt, r.used_dt);
            max_cond_C = std::max(max_cond_C, r.cond_C);
            if (std::isfinite(r.raw_cond)) max_raw_cond = std::max(max_raw_cond, r.raw_cond);
            if (std::isfinite(r.actual_solve_cond)) {
                max_actual_solve_cond = std::max(max_actual_solve_cond, r.actual_solve_cond);
            }
            min_actual_solve_rank = std::min(min_actual_solve_rank, r.effective_rank);
            if (std::isfinite(r.relative_raw_residual)) {
                max_relative_raw_residual = std::max(max_relative_raw_residual, r.relative_raw_residual);
            }
            if (std::isfinite(r.discarded_rhs_fraction)) {
                max_discarded_rhs_fraction = std::max(max_discarded_rhs_fraction, r.discarded_rhs_fraction);
            }
            for (int s = 0; s < 3; ++s) last_sv_small[s] = r.sv_small[s];

            const double norm_now = ecg1d::overlap(basis).real();
            if (norm_now > 1e-15 && std::isfinite(norm_now)) {
                const double scale = std::sqrt(initial_norm / norm_now);
                for (auto& bp : basis) bp.u *= scale;
            }
            ensure_finite_basis("after norm enforcement", step + 1, t + dt_step);

            const ecg1d::WidthMonitors wmon = ecg1d::compute_width_monitors(basis);
            min_re_B_run = std::min(min_re_B_run, wmon.min_re_B);
            min_re_AplusB_run = std::min(min_re_AplusB_run, wmon.min_re_AplusB);

            t += dt_step;
            step++;
            const double phi = phi_at(t, cfg);
            const bool is_final = (t >= cfg.total_time - 1e-15);
            const bool do_snapshot = is_final ||
                (cfg.snapshot_every > 0 && step % cfg.snapshot_every == 0);
            const bool do_sample = is_final ||
                (cfg.trace_every > 0 && step % cfg.trace_every == 0);

            if (do_snapshot) snapshots.save("tdvp_step", step, t, phi, basis);
            if (do_sample) {
                HamiltonianTerms terms_sample = hamiltonian_at(t, cfg, opt);
                ecg1d::sample_observables(basis, t, phi, terms_sample, cfg.lattice_a, trace);
                ecg1d::append_trace_diagnostics(trace, wmon, &r);
                sample_n2_diagnostics(basis, t, phi, terms_sample, n2_diag);
                const size_t idx = trace.t.size() - 1;
                std::cout << "step " << std::setw(5) << step
                          << " t=" << std::fixed << std::setprecision(5) << t
                          << " tau=" << phi
                          << " E=" << std::setprecision(8) << trace.E_total[idx]
                          << " <X>/a=" << trace.polarization_cell[idx]
                          << " norm=" << trace.norm[idx]
                          << " raw_cond=" << std::scientific << std::setprecision(2) << r.raw_cond
                          << " solve_cond=" << r.actual_solve_cond
                          << " rank=" << r.effective_rank << "/" << param_dim
                          << " raw_resid=" << r.relative_raw_residual
                          << " dt=" << r.used_dt
                          << " |dz|=" << r.dz_norm
                          << " minReB=" << wmon.min_re_B << "@" << wmon.argmin_re_B
                          << " minRe(A+B)=" << wmon.min_re_AplusB << "@" << wmon.argmin_re_AplusB;
                if (!trace.V_gauss.empty()) {
                    std::cout << " V_gauss=" << trace.V_gauss[idx];
                }
                std::cout << std::defaultfloat << "\n";
            }
        }

        const auto evolution_wall_end = std::chrono::steady_clock::now();
        const double evolution_wall_seconds =
            std::chrono::duration<double>(evolution_wall_end - evolution_wall_start).count();
        const double evolution_seconds_per_step =
            (step > 0) ? evolution_wall_seconds / static_cast<double>(step) : dnan;

        snapshots.write();
        ecg1d::write_trace_csv(task_dir + "/trace.csv", trace);
        write_n2_diagnostics(task_dir + "/n2_trace.csv", n2_diag);
        write_basis_csv(task_dir + "/basis_final.csv", basis);

        const double E0 = trace.E_total.front();
        const double E1 = trace.E_total.back();
        const double energy_drift_abs = E1 - E0;
        const double energy_drift_rel =
            (std::abs(E0) > 0.0) ? std::abs(energy_drift_abs) / std::abs(E0) : dnan;
        const double norm0 = trace.norm.front();
        const double norm1 = trace.norm.back();
        const double norm_drift_abs = norm1 - norm0;
        const double norm_drift_rel =
            (std::abs(norm0) > 0.0) ? std::abs(norm_drift_abs) / std::abs(norm0) : dnan;
        const double x2_start = trace.x2.front();
        const double x2_end = trace.x2.back();
        const double x2_growth_ratio = (std::abs(x2_start) > 0.0) ? x2_end / x2_start : dnan;
        if (min_actual_solve_rank == std::numeric_limits<int>::max()) min_actual_solve_rank = -1;

        ecg1d::RunSummaryStats summary_stats;
        summary_stats.param_dim = param_dim;
        summary_stats.steps_total = step;
        summary_stats.min_accepted_dt = min_accepted_dt;
        summary_stats.max_accepted_dt = max_accepted_dt;
        summary_stats.evolution_wall_seconds = evolution_wall_seconds;
        summary_stats.evolution_seconds_per_step = evolution_seconds_per_step;
        summary_stats.max_raw_cond = max_raw_cond;
        summary_stats.max_actual_solve_cond = max_actual_solve_cond;
        summary_stats.max_cond_C = max_cond_C;
        summary_stats.min_actual_solve_rank = min_actual_solve_rank;
        summary_stats.final_sv_small[0] = last_sv_small[0];
        summary_stats.final_sv_small[1] = last_sv_small[1];
        summary_stats.final_sv_small[2] = last_sv_small[2];
        summary_stats.max_relative_raw_residual = max_relative_raw_residual;
        summary_stats.max_discarded_rhs_fraction = max_discarded_rhs_fraction;
        summary_stats.min_re_B_run = min_re_B_run;
        summary_stats.min_re_AplusB_run = min_re_AplusB_run;
        ecg1d::write_run_summary_txt(
            task_dir + "/summary.txt", cfg, trace, summary_stats, snapshots.event_idx, snapshots.row_count);

        std::cout << "[write] " << task_dir << "/trace.csv"
                  << " (" << trace.t.size() << " samples)\n";
        if (!n2_diag.t.empty()) {
            std::cout << "[write] " << task_dir << "/n2_trace.csv"
                      << " (" << n2_diag.t.size() << " samples)\n";
        }
        std::cout << "[write] " << task_dir << "/snapshots.csv"
                  << " (" << snapshots.event_idx << " events, "
                  << snapshots.row_count << " rows)\n";
        std::cout << "[write] " << task_dir << "/basis_initial.csv\n";
        std::cout << "[write] " << task_dir << "/basis_final.csv\n";
        std::cout << "[write] " << task_dir << "/summary.txt\n";

        if (!trace.polarization_cell.empty()) {
            std::cout << "=== summary ===\n"
                      << "param_dim               = " << param_dim << "\n"
                      << "polarization_cell start = " << trace.polarization_cell.front() << "\n"
                      << "polarization_cell end   = " << trace.polarization_cell.back() << "\n"
                      << "delta polarization      = "
                      << trace.polarization_cell.back() - trace.polarization_cell.front() << "\n"
                      << std::scientific << std::setprecision(3)
                      << "max raw cond            = " << max_raw_cond << "\n"
                      << "max actual-solve cond   = " << max_actual_solve_cond << "\n"
                      << "min actual-solve rank   = " << min_actual_solve_rank
                      << " / " << param_dim << "\n"
                      << "max raw residual        = " << max_relative_raw_residual << "\n"
                      << "max discarded rhs frac  = " << max_discarded_rhs_fraction << "\n"
                      << "min Re(B)   over run    = " << min_re_B_run << "\n"
                      << "min Re(A+B) over run    = " << min_re_AplusB_run << "\n"
                      << "norm drift (rel)        = " << norm_drift_rel << "\n"
                      << "energy drift (rel)      = " << energy_drift_rel << "\n"
                      << "x2 growth ratio         = " << x2_growth_ratio << "\n"
                      << std::defaultfloat;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 2;
    }
}

}  // namespace pump2
