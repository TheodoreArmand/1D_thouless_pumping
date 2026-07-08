#pragma once

#include "basis_params.hpp"
#include "realtime_tdvp.hpp"
#include "tdvp_solver.hpp"

#include <limits>
#include <string>
#include <vector>

namespace ecg1d {

struct Trace {
    std::vector<double> t;
    std::vector<double> tau;
    std::vector<double> norm;
    std::vector<double> E_total;
    std::vector<double> T_kin;
    std::vector<double> V_cos;
    std::vector<double> V_const;
    std::vector<double> V_lattice;
    // Present only for N >= 2 samples. write_trace_csv appends the column
    // when this vector is populated, preserving the legacy N=1 trace schema.
    std::vector<double> V_gauss;
    std::vector<double> x_mean;
    std::vector<double> p_mean;
    std::vector<double> x2;
    std::vector<double> p2;
    std::vector<double> polarization_cell;
    std::vector<int> basis_size;

    std::vector<double> raw_cond;
    std::vector<double> actual_solve_cond;
    std::vector<int> actual_solve_rank;
    std::vector<double> sv_min;
    std::vector<double> relative_raw_residual;
    std::vector<double> discarded_rhs_fraction;
    std::vector<double> dz_norm;
    std::vector<double> metric_norm;
    std::vector<double> min_re_B;
    std::vector<int> argmin_re_B;
    std::vector<double> min_re_AplusB;
    std::vector<int> argmin_re_AplusB;
};

struct WidthMonitors {
    double min_re_B = std::numeric_limits<double>::quiet_NaN();
    int argmin_re_B = -1;
    double min_re_AplusB = std::numeric_limits<double>::quiet_NaN();
    int argmin_re_AplusB = -1;
};

WidthMonitors compute_width_monitors(const std::vector<BasisParams>& basis);

void sample_observables(const std::vector<BasisParams>& basis,
                        double t,
                        double phi,
                        const HamiltonianTerms& terms,
                        double lattice_a,
                        Trace& trace);

void append_trace_diagnostics(Trace& trace,
                              const WidthMonitors& width,
                              const RealtimeStepResult* step_result);

void write_trace_csv(const std::string& path, const Trace& trace);

}  // namespace ecg1d
