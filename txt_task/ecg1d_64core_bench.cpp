#define main ecg1d_locked_default_main
#include "../main.cpp"
#undef main

#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::string n1_task_dir(const PumpConfig& cfg) {
    std::string task_dir = cfg.out_root + "/a" + format_output_tag(cfg.lattice_a)
                         + "_K" + std::to_string(cfg.K)
                         + "_tmax" + format_output_tag(cfg.total_time);
    task_dir += "_VsER" + format_output_tag(cfg.Vs / cfg.recoil_energy)
              + "_VlER" + format_output_tag(cfg.Vl / cfg.recoil_energy);
    if (!cfg.output_protocol_tag.empty()) {
        task_dir += "_" + cfg.output_protocol_tag;
    }
    return task_dir;
}

int run_n1_case(PumpConfig cfg, const std::string& out_root, int steps) {
    cfg.total_time = static_cast<double>(steps) * cfg.dt;
    cfg.trace_every = 1;
    cfg.snapshot_every = 0;
    cfg.out_root = out_root;

    load_phase_schedule(cfg);
    validate_config(cfg);

    std::vector<BasisParams> basis = load_basis_csv(cfg.basis_from, cfg.N);
    if (static_cast<int>(basis.size()) != cfg.K) {
        throw std::runtime_error("basis CSV size does not match locked K");
    }
    normalize_basis(basis);

    const std::string task_dir = n1_task_dir(cfg);
    std::filesystem::create_directories(task_dir);
    write_config_txt(task_dir + "/config.txt", cfg);
    write_basis_csv(task_dir + "/basis_initial.csv", basis);

    SolverConfig solver_cfg;
    solver_cfg.lambda_C = cfg.lambda_C;
    solver_cfg.rcond = cfg.rcond;

    std::vector<AlphaIndex> alpha = make_alpha_list(cfg.N, cfg.K);
    const int param_dim = static_cast<int>(alpha.size());
    const long long steps_total_est = estimate_time_steps(cfg);
    if (param_dim == 0) throw std::runtime_error("TDVP parameter list is empty");

    SnapshotSaver snapshots(task_dir);
    Trace trace;
    N1ProgressWriter progress(task_dir);

    const double dnan = std::numeric_limits<double>::quiet_NaN();
    double max_cond_C = 0.0;
    double max_actual_solve_cond = 0.0;
    double max_solve_sv_max = 0.0;
    double last_sv_max = dnan;
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
    const double initial_norm = overlap(basis).real();
    sample_observables(
        basis, t, experimental_phi_at(t, cfg), pump_hamiltonian_at(t, cfg), cfg.lattice_a, trace);
    {
        WidthMonitors w0 = compute_width_monitors(basis);
        min_re_B_run = std::min(min_re_B_run, w0.min_re_B);
        min_re_AplusB_run = std::min(min_re_AplusB_run, w0.min_re_AplusB);
        append_trace_diagnostics(trace, w0, nullptr);
    }
    snapshots.save("initial", step, t, 0.0, basis);
    progress.write(step, steps_total_est, cfg.total_time, 0.0, 0.0, param_dim, trace, 0);

    const auto evolution_wall_start = std::chrono::steady_clock::now();
    while (t < cfg.total_time - 1e-15) {
        const double dt_step = std::min(cfg.dt, cfg.total_time - t);
        auto terms_at = [&](double stage_t) {
            return pump_hamiltonian_at(stage_t, cfg);
        };

        RealtimeStepResult r = realtime_tdvp_step_rk4_time_dependent(
            alpha, basis, t, dt_step, terms_at, solver_cfg);
        basis = std::move(r.basis);

        min_accepted_dt = std::min(min_accepted_dt, r.used_dt);
        max_accepted_dt = std::max(max_accepted_dt, r.used_dt);
        max_cond_C = std::max(max_cond_C, r.cond_C);
        if (std::isfinite(r.actual_solve_cond)) {
            max_actual_solve_cond = std::max(max_actual_solve_cond, r.actual_solve_cond);
        }
        if (std::isfinite(r.sv_max)) {
            max_solve_sv_max = std::max(max_solve_sv_max, r.sv_max);
            last_sv_max = r.sv_max;
        }
        min_actual_solve_rank = std::min(min_actual_solve_rank, r.effective_rank);
        if (std::isfinite(r.relative_raw_residual)) {
            max_relative_raw_residual = std::max(max_relative_raw_residual, r.relative_raw_residual);
        }
        if (std::isfinite(r.discarded_rhs_fraction)) {
            max_discarded_rhs_fraction =
                std::max(max_discarded_rhs_fraction, r.discarded_rhs_fraction);
        }
        for (int s = 0; s < 3; ++s) last_sv_small[s] = r.sv_small[s];

        const double norm_now = overlap(basis).real();
        if (norm_now > 1e-15 && std::isfinite(norm_now)) {
            const double scale = std::sqrt(initial_norm / norm_now);
            for (auto& bp : basis) bp.u *= scale;
        }

        const WidthMonitors wmon = compute_width_monitors(basis);
        min_re_B_run = std::min(min_re_B_run, wmon.min_re_B);
        min_re_AplusB_run = std::min(min_re_AplusB_run, wmon.min_re_AplusB);

        t += dt_step;
        step++;
        const double phi = experimental_phi_at(t, cfg);
        if (t >= cfg.total_time - 1e-15) snapshots.save("tdvp_step", step, t, phi, basis);

        sample_observables(basis, t, phi, pump_hamiltonian_at(t, cfg), cfg.lattice_a, trace);
        append_trace_diagnostics(trace, wmon, &r);
        const size_t idx = trace.t.size() - 1;
        const auto now = std::chrono::steady_clock::now();
        const double wall_seconds =
            std::chrono::duration<double>(now - evolution_wall_start).count();
        progress.write(
            step, steps_total_est, cfg.total_time, r.used_dt, wall_seconds, param_dim, trace, idx);
    }

    const auto evolution_wall_end = std::chrono::steady_clock::now();
    const double evolution_wall_seconds =
        std::chrono::duration<double>(evolution_wall_end - evolution_wall_start).count();
    const double evolution_seconds_per_step =
        (step > 0) ? evolution_wall_seconds / static_cast<double>(step) : dnan;

    snapshots.write();
    write_trace_csv(task_dir + "/trace.csv", trace);
    write_basis_csv(task_dir + "/basis_final.csv", basis);

    if (min_actual_solve_rank == std::numeric_limits<int>::max()) min_actual_solve_rank = -1;
    RunSummaryStats summary_stats;
    summary_stats.param_dim = param_dim;
    summary_stats.steps_total = step;
    summary_stats.min_accepted_dt = min_accepted_dt;
    summary_stats.max_accepted_dt = max_accepted_dt;
    summary_stats.evolution_wall_seconds = evolution_wall_seconds;
    summary_stats.evolution_seconds_per_step = evolution_seconds_per_step;
    summary_stats.max_actual_solve_cond = max_actual_solve_cond;
    summary_stats.max_cond_C = max_cond_C;
    summary_stats.max_solve_sv_max = max_solve_sv_max;
    summary_stats.final_sv_max = last_sv_max;
    summary_stats.min_actual_solve_rank = min_actual_solve_rank;
    summary_stats.final_sv_small[0] = last_sv_small[0];
    summary_stats.final_sv_small[1] = last_sv_small[1];
    summary_stats.final_sv_small[2] = last_sv_small[2];
    summary_stats.max_relative_raw_residual = max_relative_raw_residual;
    summary_stats.max_discarded_rhs_fraction = max_discarded_rhs_fraction;
    summary_stats.min_re_B_run = min_re_B_run;
    summary_stats.min_re_AplusB_run = min_re_AplusB_run;
    write_run_summary_txt(
        task_dir + "/summary.txt", cfg, trace, summary_stats, snapshots.event_idx, snapshots.row_count);

    std::cout << "[bench] completed " << cfg.config_name << " -> " << task_dir << "\n";
    return 0;
}

int run_n2_free_case(int steps) {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("lohse_n2_free");
    cfg.total_time = static_cast<double>(steps) * cfg.dt;
    cfg.trace_every = 1;
    cfg.snapshot_every = 0;
    cfg.out_root = "out/bench_64core_10step/n2_lohse_free";

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    opt.config_appendix =
        "driver=ecg1d_64core_bench_n2_lohse_free\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n";
    return pump2::run_pump(cfg, opt);
}

int run_n2_gauss_case(int steps) {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("lohse_n2_gauss");
    cfg.total_time = static_cast<double>(steps) * cfg.dt;
    cfg.trace_every = 1;
    cfg.snapshot_every = 0;
    cfg.out_root = "out/bench_64core_10step/n2_lohse_gauss";

    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    opt.config_appendix += "driver=ecg1d_64core_bench_n2_lohse_gauss\n";
    return pump2::run_pump(cfg, opt);
}

}  // namespace

int main() {
    try {
        constexpr int steps = 10;

        PumpConfig cfg_33 = pumpconfig::make_pump_config("legacy_prb_3_3");
        run_n1_case(cfg_33, "out/bench_64core_10step/n1_3_3", steps);

        PumpConfig cfg_lohse = pumpconfig::make_pump_config("lohse_10_5");
        run_n1_case(cfg_lohse, "out/bench_64core_10step/n1_lohse", steps);

        run_n2_free_case(steps);
        run_n2_gauss_case(steps);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FATAL bench: " << e.what() << "\n";
        return 2;
    }
}
