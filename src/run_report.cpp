#include "run_report.hpp"

#include "physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ecg1d {
namespace {

double lattice_constant_offset_per_particle() {
    return 0.0;
}

double short_cosine_coeff(const pumpconfig::PumpConfig& cfg) {
    return -cfg.Vs;
}

double long_cosine_coeff(const pumpconfig::PumpConfig& cfg) {
    return -cfg.Vl;
}

double manual_fine_dt_factor_at_s(const pumpconfig::PumpConfig& cfg, double s) {
    if (!cfg.manual_fine_dt_enabled) return 1.0;
    double factor = 1.0;
    for (const auto& w : cfg.manual_fine_dt_windows_s) {
        if (s >= w.first && s < w.second) {
            factor = std::min(factor, cfg.manual_fine_dt_factor);
        }
    }
    for (const auto& w : cfg.manual_fine_dt_schedule_s) {
        if (s >= w.start_s && s < w.end_s) {
            factor = std::min(factor, w.factor);
        }
    }
    return factor;
}

double manual_fine_dt_min_factor(const pumpconfig::PumpConfig& cfg) {
    if (!cfg.manual_fine_dt_enabled) return 1.0;
    double factor = 1.0;
    if (!cfg.manual_fine_dt_windows_s.empty()) {
        factor = std::min(factor, cfg.manual_fine_dt_factor);
    }
    for (const auto& w : cfg.manual_fine_dt_schedule_s) {
        factor = std::min(factor, w.factor);
    }
    return factor;
}

std::string time_step_mode_string(const pumpconfig::PumpConfig& cfg) {
    if (!cfg.manual_fine_dt_enabled) return "fixed";
    if (!cfg.manual_fine_dt_schedule_s.empty()) return "manual_fine_schedule";
    return "manual_fine_window";
}

std::string manual_fine_windows_string(const pumpconfig::PumpConfig& cfg) {
    std::ostringstream oss;
    oss << std::setprecision(17);
    for (size_t i = 0; i < cfg.manual_fine_dt_windows_s.size(); ++i) {
        if (i) oss << ";";
        oss << cfg.manual_fine_dt_windows_s[i].first
            << "-" << cfg.manual_fine_dt_windows_s[i].second;
    }
    return oss.str();
}

std::string manual_fine_schedule_string(const pumpconfig::PumpConfig& cfg) {
    std::ostringstream oss;
    oss << std::setprecision(17);
    for (size_t i = 0; i < cfg.manual_fine_dt_schedule_s.size(); ++i) {
        if (i) oss << ";";
        const auto& w = cfg.manual_fine_dt_schedule_s[i];
        oss << w.start_s << "-" << w.end_s << "@" << w.factor;
    }
    return oss.str();
}

}  // namespace

std::string format_output_tag(double x) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << x;
    std::string s = oss.str();
    for (char& c : s) {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

long long estimate_time_steps(const pumpconfig::PumpConfig& cfg) {
    if (!cfg.manual_fine_dt_enabled) {
        return static_cast<long long>(std::ceil(cfg.total_time / cfg.dt));
    }

    long long steps = 0;
    double t = 0.0;
    constexpr long long max_steps = 1000000000LL;
    while (t < cfg.total_time - 1e-15) {
        const double dt_step = next_time_step(cfg, t);
        if (!(dt_step > 0.0) || !std::isfinite(dt_step)) break;
        t += dt_step;
        steps++;
        if (steps > max_steps) break;
    }
    return steps;
}

double time_step_at(const pumpconfig::PumpConfig& cfg, double t) {
    if (!cfg.manual_fine_dt_enabled) return cfg.dt;
    const double s = (cfg.pump_period > 0.0) ? (t / cfg.pump_period) : 0.0;
    return cfg.dt * manual_fine_dt_factor_at_s(cfg, s);
}

double next_time_step(const pumpconfig::PumpConfig& cfg, double t) {
    double dt_step = time_step_at(cfg, t);
    if (cfg.manual_fine_dt_enabled) {
        auto clamp_to_boundary = [&](double boundary_s) {
            const double boundary_t = boundary_s * cfg.pump_period;
            if (t < boundary_t && t + dt_step > boundary_t) {
                dt_step = boundary_t - t;
            }
        };
        for (const auto& w : cfg.manual_fine_dt_windows_s) {
            clamp_to_boundary(w.first);
            clamp_to_boundary(w.second);
        }
        for (const auto& w : cfg.manual_fine_dt_schedule_s) {
            clamp_to_boundary(w.start_s);
            clamp_to_boundary(w.end_s);
        }
    }
    return std::min(dt_step, cfg.total_time - t);
}

void write_config_txt(const std::string& path, const pumpconfig::PumpConfig& cfg) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    f << std::setprecision(17)
      << "config_name=" << cfg.config_name << "\n"
      << "N=" << cfg.N << "\n"
      << "K=" << cfg.K << "\n"
      << "model=" << cfg.model << "\n"
      << "potential=" << cfg.potential_label << "\n"
      << "reference_units=" << cfg.reference_units_label << "\n"
      << "hbar=" << hbar << "\n"
      << "mass=" << mass << "\n"
      << "hbar2_over_2m=" << hbar * hbar / (2.0 * mass) << "\n"
      << "d_s_code=" << 0.5 * cfg.lattice_a << "\n"
      << "a_code_long_period=d_l_code=" << cfg.lattice_a << "\n"
      << "pump_period=" << cfg.pump_period << "\n"
      << "total_time=" << cfg.total_time << "\n"
      << "dt=" << cfg.dt << "\n"
      << "time_step_mode=" << time_step_mode_string(cfg) << "\n"
      << "manual_fine_dt_enabled=" << (cfg.manual_fine_dt_enabled ? 1 : 0) << "\n"
      << "manual_fine_dt_factor=" << cfg.manual_fine_dt_factor << "\n"
      << "manual_fine_dt_windows_s=" << manual_fine_windows_string(cfg) << "\n"
      << "manual_fine_dt_schedule_s=" << manual_fine_schedule_string(cfg) << "\n"
      << "manual_fine_dt_min_dt=" << cfg.dt * manual_fine_dt_min_factor(cfg) << "\n"
      << "phase_schedule_csv=" << cfg.phase_schedule_csv << "\n"
      << "phase_schedule_convention=" << cfg.phase_schedule_convention << "\n"
      << "output_protocol_tag=" << cfg.output_protocol_tag << "\n"
      << "phase_schedule_points=" << cfg.phase_s.size() << "\n"
      << "estimated_time_steps=" << estimate_time_steps(cfg) << "\n"
      << "recoil_energy_code_Ers=" << cfg.recoil_energy << "\n"
      << "Er_l_over_Er_s=1\n"
      << "Vs_depth_over_Er_s=" << cfg.Vs_depth_over_Er << "\n"
      << "Vl_depth_over_Er_l=" << cfg.Vl_depth_over_Er << "\n"
      << "Vl_depth_over_Er_s=" << cfg.Vl_depth_over_Er << "\n"
      << "Vs_depth_code=" << cfg.Vs << "\n"
      << "Vl_depth_code=" << cfg.Vl << "\n"
      << "cos_short_coeff_over_Er_s=" << short_cosine_coeff(cfg) / cfg.recoil_energy << "\n"
      << "cos_long_coeff_over_Er_s=" << long_cosine_coeff(cfg) / cfg.recoil_energy << "\n"
      << "constant_offset_per_particle_over_Er_s="
      << lattice_constant_offset_per_particle() / cfg.recoil_energy << "\n"
      << "phase=phi=interpolate_phase_schedule(t/pump_period)\n"
      << "expanded_potential_without_constant=" << cfg.expanded_potential_label << "\n"
      << "lambda_C=" << cfg.lambda_C << "\n"
      << "rcond=" << cfg.rcond << "\n";
    if (cfg.N >= 2) {
        f << "g_gauss_over_Er_s=" << cfg.g_gauss_over_Er << "\n"
          << "sigma_gauss_code=" << cfg.sigma_gauss << "\n";
    }
    f << "snapshot_every=" << cfg.snapshot_every << "\n"
      << "trace_every=" << cfg.trace_every << "\n"
      << "snapshot_buffered=1"
      << "  # snapshots buffered in memory, snapshots.csv written once at end\n"
      << "snapshot_mode="
      << (cfg.snapshot_every == 0 ? "initial_final_only" : "every_N_steps") << "\n"
      << "basis_from=" << cfg.basis_from << "\n";
}

void write_run_summary_txt(const std::string& path,
                           const pumpconfig::PumpConfig& cfg,
                           const Trace& trace,
                           const RunSummaryStats& stats,
                           int snapshot_events,
                           int snapshot_rows) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);

    const double dnan = std::numeric_limits<double>::quiet_NaN();
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

    f << std::setprecision(17)
      << "config_name=" << cfg.config_name << "\n"
      << "model=" << cfg.model << "\n"
      << "potential=" << cfg.potential_label << "\n"
      << "hbar2_over_2m=" << hbar * hbar / (2.0 * mass) << "\n"
      << "recoil_energy=" << cfg.recoil_energy << "\n"
      << "lambda_C=" << cfg.lambda_C << "\n"
      << "rcond=" << cfg.rcond << "\n"
      << "param_dim=" << stats.param_dim << "\n"
      << "time_step_mode=" << time_step_mode_string(cfg) << "\n"
      << "dt=" << cfg.dt << "\n"
      << "manual_fine_dt_enabled=" << (cfg.manual_fine_dt_enabled ? 1 : 0) << "\n"
      << "manual_fine_dt_factor=" << cfg.manual_fine_dt_factor << "\n"
      << "manual_fine_dt_windows_s=" << manual_fine_windows_string(cfg) << "\n"
      << "manual_fine_dt_schedule_s=" << manual_fine_schedule_string(cfg) << "\n"
      << "manual_fine_dt_min_dt=" << cfg.dt * manual_fine_dt_min_factor(cfg) << "\n"
      << "phase_schedule_csv=" << cfg.phase_schedule_csv << "\n"
      << "phase_schedule_convention=" << cfg.phase_schedule_convention << "\n"
      << "output_protocol_tag=" << cfg.output_protocol_tag << "\n"
      << "phase_schedule_points=" << cfg.phase_s.size() << "\n"
      << "estimated_time_steps=" << estimate_time_steps(cfg) << "\n"
      << "steps_total=" << stats.steps_total << "\n"
      << "min_accepted_dt="
      << (std::isfinite(stats.min_accepted_dt) ? stats.min_accepted_dt : dnan) << "\n"
      << "max_accepted_dt=" << stats.max_accepted_dt << "\n"
      << "evolution_wall_seconds=" << stats.evolution_wall_seconds << "\n"
      << "evolution_seconds_per_step=" << stats.evolution_seconds_per_step << "\n"
      << "snapshot_every=" << cfg.snapshot_every << "\n"
      << "trace_every=" << cfg.trace_every << "\n"
      << "snapshot_buffered=1\n"
      << "snapshot_mode="
      << (cfg.snapshot_every == 0 ? "initial_final_only" : "every_N_steps") << "\n"
      << "snapshot_events=" << snapshot_events << "\n"
      << "snapshot_rows=" << snapshot_rows << "\n"
      << "trace_samples=" << trace.t.size() << "\n"
      << "polarization_cell_start=" << trace.polarization_cell.front() << "\n"
      << "polarization_cell_end=" << trace.polarization_cell.back() << "\n"
      << "delta_polarization="
      << trace.polarization_cell.back() - trace.polarization_cell.front() << "\n"
      << "max_raw_cond=nan  # disabled; raw C condition is not computed\n"
      << "max_actual_solve_cond=" << stats.max_actual_solve_cond << "\n"
      << "max_cond_C=" << stats.max_cond_C << "  # == max_actual_solve_cond (legacy name)\n"
      << "max_solve_sv_max=" << stats.max_solve_sv_max << "\n"
      << "final_sv_max=" << stats.final_sv_max << "\n"
      << "min_actual_solve_rank=" << stats.min_actual_solve_rank << "\n"
      << "final_sv_small_0=" << stats.final_sv_small[0] << "\n"
      << "final_sv_small_1=" << stats.final_sv_small[1] << "\n"
      << "final_sv_small_2=" << stats.final_sv_small[2] << "\n"
      << "max_relative_raw_residual=" << stats.max_relative_raw_residual << "\n"
      << "max_discarded_rhs_fraction=" << stats.max_discarded_rhs_fraction << "\n"
      << "width_monitor_AplusB=min_lambda_Re_AplusB\n"
      << "min_re_B_run=" << stats.min_re_B_run << "\n"
      << "min_re_B_final=" << trace.min_re_B.back()
      << " @kN_plus_particle " << trace.argmin_re_B.back() << "\n"
      << "min_re_AplusB_run=" << stats.min_re_AplusB_run << "\n"
      << "min_re_AplusB_final=" << trace.min_re_AplusB.back()
      << " @gaussian " << trace.argmin_re_AplusB.back() << "\n"
      << "norm_start=" << norm0 << "\n"
      << "norm_end=" << norm1 << "\n"
      << "norm_drift_abs=" << norm_drift_abs << "\n"
      << "norm_drift_rel=" << norm_drift_rel << "\n"
      << "energy_start=" << E0 << "\n"
      << "energy_end=" << E1 << "\n"
      << "energy_drift_abs=" << energy_drift_abs << "\n"
      << "energy_drift_rel=" << energy_drift_rel << "\n"
      << "x2_start=" << x2_start << "\n"
      << "x2_end=" << x2_end << "\n"
      << "x2_growth_ratio=" << x2_growth_ratio << "\n";
    if (!trace.V_gauss.empty()) {
        f << "V_gauss_start=" << trace.V_gauss.front() << "\n"
          << "V_gauss_end=" << trace.V_gauss.back() << "\n";
    }
}

}  // namespace ecg1d
