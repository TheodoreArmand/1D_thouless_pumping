#pragma once

#include "pumpconfig/pump_config.hpp"
#include "trace_io.hpp"

#include <cstddef>
#include <fstream>
#include <string>

namespace ecg1d {

struct RunSummaryStats {
    int param_dim = 0;
    int steps_total = 0;
    double min_accepted_dt = 0.0;
    double max_accepted_dt = 0.0;
    double evolution_wall_seconds = 0.0;
    double evolution_seconds_per_step = 0.0;
    double max_actual_solve_cond = 0.0;
    double max_cond_C = 0.0;
    double max_solve_sv_max = 0.0;
    double final_sv_max = 0.0;
    int min_actual_solve_rank = -1;
    double final_sv_small[3] = {0.0, 0.0, 0.0};
    double max_relative_raw_residual = 0.0;
    double max_discarded_rhs_fraction = 0.0;
    double min_re_B_run = 0.0;
    double min_re_AplusB_run = 0.0;
};

std::string format_output_tag(double x);

long long estimate_time_steps(const pumpconfig::PumpConfig& cfg);

double time_step_at(const pumpconfig::PumpConfig& cfg, double t);

double next_time_step(const pumpconfig::PumpConfig& cfg, double t);

void write_config_txt(const std::string& path, const pumpconfig::PumpConfig& cfg);

void write_run_summary_txt(const std::string& path,
                           const pumpconfig::PumpConfig& cfg,
                           const Trace& trace,
                           const RunSummaryStats& stats,
                           int snapshot_events,
                           int snapshot_rows);

class N1ProgressWriter {
public:
    explicit N1ProgressWriter(const std::string& task_dir);

    void write(int step,
               long long steps_total_est,
               double total_time,
               double accepted_dt,
               double wall_seconds,
               int param_dim,
               const Trace& trace,
               size_t idx);

private:
    std::ofstream csv_;
};

class N2ProgressWriter {
public:
    explicit N2ProgressWriter(const std::string& task_dir);

    void write(int step,
               long long steps_total_est,
               double total_time,
               double accepted_dt,
               double wall_seconds,
               int param_dim,
               const Trace& trace,
               size_t idx,
               double r12_sq,
               double v_gauss);

private:
    std::ofstream csv_;
};

}  // namespace ecg1d
