#pragma once

#include <string>
#include <vector>

namespace pumpconfig {

constexpr double pi = 3.141592653589793238462643383279502884;

struct ReferenceGap12 {
    double phi0 = 0.0;
    double phipi = 0.0;
    double max = 0.0;
    bool enabled = false;
};

struct PumpConfig {
    std::string config_name;
    std::string model;

    int N = 1;
    int K = 16;

    double mass_code_units = 1.0;
    double lattice_a = 8.0;  // long lattice period d_l; short spacing d_s=a/2
    double recoil_energy = 0.0;
    double Vs_depth_over_Er = 0.0;
    double Vl_depth_over_Er = 0.0;
    double Vs = 0.0;
    double Vl = 0.0;
    double omega = 0.0;  // harmonic trap frequency (0 = off)

    double pump_period = 0.0;
    double total_time = 0.0;
    double dt = 0.0;
    bool dynamic_dt = false;
    double dt_min = 1.0;
    double dt_max = 50.0;
    double dt_gap_power = 2.0;
    std::string integrator = "rk4";
    double adaptive_abs_tol = 1.0e-7;
    double adaptive_rel_tol = 1.0e-4;
    double adaptive_dt_min = 1.0e-6;
    int adaptive_max_reject = 20;

    // Optional physical phase schedule. "linear" keeps phi=2*pi*t/T.
    // "local_csv" loads a CSV with columns s=t/T and phi(s).
    std::string phase_schedule = "linear";
    std::string phase_schedule_csv;
    double dphi_max = 0.0;  // optional numerical cap on phase advance per TDVP step
    std::vector<double> phase_s;
    std::vector<double> phase_phi;

    // snapshot_every: write a per-Gaussian snapshot every N steps.
    //   N > 0 : intermediate snapshots every N steps (plus initial + final).
    //   N == 0: no intermediate snapshots, only the initial and final events.
    int snapshot_every = 25;
    // trace_every: record observables (trace.csv) + stdout line every N steps.
    // Decoupled from snapshot_every so snapshots can be made sparse (or off)
    // without thinning the trace, which is the primary health output.
    int trace_every = 25;

    double lambda_C = 1.0e-8;
    double rcond = 1.0e-4;
    bool u_split_trotter = false;
    bool enforce_norm = true;
    bool autoinit = false;

    std::string basis_from;
    std::string out_root;

    std::string potential_label;
    std::string reference_units_label;
    std::string expanded_potential_label;
    double lambda_l_nm = 0.0;
    double lambda_s_nm = 0.0;
    ReferenceGap12 reference_gap12;
};

double recoil_energy_for_lattice(double lattice_a, double mass_code_units);

void set_lattice_depths_over_recoil(PumpConfig& cfg,
                                    double vs_depth_over_Er,
                                    double vl_depth_over_Er);

PumpConfig make_pump_config(const std::string& name);

std::vector<std::string> available_pump_configs();

}  // namespace pumpconfig
