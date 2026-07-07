#pragma once

#include "pumpconfig/pump_config.hpp"
#include "trace_io.hpp"

#include <string>

namespace ecg1d {

struct RunSummaryStats {
    int param_dim = 0;
    int steps_total = 0;
    double min_accepted_dt = 0.0;
    double max_accepted_dt = 0.0;
    double evolution_wall_seconds = 0.0;
    double evolution_seconds_per_step = 0.0;
    double max_raw_cond = 0.0;
    double max_actual_solve_cond = 0.0;
    double max_cond_C = 0.0;
    int min_actual_solve_rank = -1;
    double final_sv_small[3] = {0.0, 0.0, 0.0};
    double max_relative_raw_residual = 0.0;
    double max_discarded_rhs_fraction = 0.0;
    double min_re_B_run = 0.0;
    double min_re_AplusB_run = 0.0;
};

std::string format_output_tag(double x);

long long estimate_time_steps(const pumpconfig::PumpConfig& cfg);

void write_config_txt(const std::string& path, const pumpconfig::PumpConfig& cfg);

void write_run_summary_txt(const std::string& path,
                           const pumpconfig::PumpConfig& cfg,
                           const Trace& trace,
                           const RunSummaryStats& stats,
                           int snapshot_events,
                           int snapshot_rows);

}  // namespace ecg1d
