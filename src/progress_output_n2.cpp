#include "run_report.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace ecg1d {
namespace {

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

double progress_fraction(int step, long long steps_total_est, double t, double total_time) {
    if (steps_total_est > 0) {
        return static_cast<double>(step) / static_cast<double>(steps_total_est);
    }
    return (total_time > 0.0) ? t / total_time : quiet_nan();
}

double seconds_per_step(int step, double wall_seconds) {
    return (step > 0) ? wall_seconds / static_cast<double>(step) : quiet_nan();
}

double eta_seconds(int step, long long steps_total_est, double wall_seconds) {
    if (step <= 0 || steps_total_est <= step) return quiet_nan();
    return seconds_per_step(step, wall_seconds) *
           static_cast<double>(steps_total_est - step);
}

double r12_rms(double r12_sq) {
    return (r12_sq >= 0.0 && std::isfinite(r12_sq)) ? std::sqrt(r12_sq) : quiet_nan();
}

void write_header(std::ofstream& f) {
    f << "step,steps_total_est,progress_frac,t,total_time,phi,accepted_dt,"
         "wall_seconds,seconds_per_step,eta_seconds,param_dim,"
         "norm,E_total,T_kin,V_cos,V_const,V_lattice,V_gauss,"
         "x_mean,p_mean,x2,p2,polarization_cell,delta_polarization,"
         "r12_sq,r12_rms,"
         "raw_cond,actual_solve_cond,actual_solve_rank,sv_max,sv_min,"
         "relative_raw_residual,discarded_rhs_fraction,"
         "dz_norm,metric_norm,min_re_B,min_re_AplusB\n";
}

void validate_idx(const Trace& trace, size_t idx) {
    if (idx >= trace.t.size()) {
        throw std::runtime_error("N2ProgressWriter: trace index out of range");
    }
}

}  // namespace

N2ProgressWriter::N2ProgressWriter(const std::string& task_dir)
    : csv_(task_dir + "/progress.csv") {
    if (!csv_.is_open()) {
        throw std::runtime_error("cannot open " + task_dir + "/progress.csv");
    }
    write_header(csv_);
    csv_ << std::setprecision(17);
    csv_.flush();
}

void N2ProgressWriter::write(int step,
                             long long steps_total_est,
                             double total_time,
                             double accepted_dt,
                             double wall_seconds,
                             int param_dim,
                             const Trace& trace,
                             size_t idx,
                             double r12_sq,
                             double v_gauss) {
    validate_idx(trace, idx);

    const double t = trace.t[idx];
    const double frac = progress_fraction(step, steps_total_est, t, total_time);
    const double s_per_step = seconds_per_step(step, wall_seconds);
    const double eta = eta_seconds(step, steps_total_est, wall_seconds);
    const double delta_p = trace.polarization_cell[idx] - trace.polarization_cell.front();
    const double r12 = r12_rms(r12_sq);

    csv_ << step << ","
         << steps_total_est << ","
         << frac << ","
         << t << ","
         << total_time << ","
         << trace.tau[idx] << ","
         << accepted_dt << ","
         << wall_seconds << ","
         << s_per_step << ","
         << eta << ","
         << param_dim << ","
         << trace.norm[idx] << ","
         << trace.E_total[idx] << ","
         << trace.T_kin[idx] << ","
         << trace.V_cos[idx] << ","
         << trace.V_const[idx] << ","
         << trace.V_lattice[idx] << ","
         << v_gauss << ","
         << trace.x_mean[idx] << ","
         << trace.p_mean[idx] << ","
         << trace.x2[idx] << ","
         << trace.p2[idx] << ","
         << trace.polarization_cell[idx] << ","
         << delta_p << ","
         << r12_sq << ","
         << r12 << ","
         << trace.raw_cond[idx] << ","
         << trace.actual_solve_cond[idx] << ","
         << trace.actual_solve_rank[idx] << ","
         << trace.sv_max[idx] << ","
         << trace.sv_min[idx] << ","
         << trace.relative_raw_residual[idx] << ","
         << trace.discarded_rhs_fraction[idx] << ","
         << trace.dz_norm[idx] << ","
         << trace.metric_norm[idx] << ","
         << trace.min_re_B[idx] << ","
         << trace.min_re_AplusB[idx] << "\n";
    csv_.flush();

    std::cout << "progress N=2 step " << step << "/" << steps_total_est
              << " " << std::fixed << std::setprecision(1) << (100.0 * frac) << "%"
              << " t=" << std::setprecision(5) << t
              << " eta_s=" << std::setprecision(1) << eta
              << " E=" << std::setprecision(8) << trace.E_total[idx]
              << " P=" << trace.polarization_cell[idx]
              << " dP=" << delta_p
              << " r12=" << r12
              << " Vg=" << v_gauss
              << " norm=" << trace.norm[idx]
              << " rank=" << trace.actual_solve_rank[idx] << "/" << param_dim
              << " sv_max=" << std::scientific << std::setprecision(2)
              << trace.sv_max[idx]
              << " resid=" << trace.relative_raw_residual[idx]
              << std::defaultfloat << "\n" << std::flush;
}

}  // namespace ecg1d
