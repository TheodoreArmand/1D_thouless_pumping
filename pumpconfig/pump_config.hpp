#pragma once

#include <string>
#include <utility>
#include <vector>

namespace pumpconfig {

constexpr double pi = 3.141592653589793238462643383279502884;

struct ManualFineDtWindow {
    double start_s = 0.0;
    double end_s = 0.0;
    double factor = 1.0;
};

struct PumpConfig {
    std::string config_name;
    std::string model;

    int N = 1;
    int K = 16;

    // Lattice and energy scale.
    double lattice_a = 8.0;  // long lattice period d_l; short spacing d_s=a/2
    double recoil_energy = 0.0;
    double Vs_depth_over_Er = 0.0;
    double Vl_depth_over_Er = 0.0;
    double Vs = 0.0;
    double Vl = 0.0;

    // Time grid.
    double pump_period = 0.0;
    double total_time = 0.0;
    double dt = 0.0;
    bool manual_fine_dt_enabled = false;
    double manual_fine_dt_factor = 1.0;
    std::vector<std::pair<double, double>> manual_fine_dt_windows_s;
    std::vector<ManualFineDtWindow> manual_fine_dt_schedule_s;

    // Phase schedule CSV with columns s=t/T and phi(s).
    std::string phase_schedule_csv;
    // Human-readable convention recorded in config.txt/summary.txt.
    std::string phase_schedule_convention;
    // Optional suffix that separates output produced by incompatible protocols.
    std::string output_protocol_tag;
    std::vector<double> phase_s;
    std::vector<double> phase_phi;

    // Output cadence.
    int snapshot_every = 25;
    int trace_every = 25;

    // TDVP linear solve.
    double lambda_C = 1.0e-8;
    double rcond = 1.0e-4;

    // Gaussian pair interaction settings for N >= 2 drivers.
    // Coupling is specified in the same short-lattice recoil unit as the
    // lattice depths and converted to code units by the driver.
    double g_gauss_over_Er = 0.0;
    double sigma_gauss = 1.0;

    // Input/output paths.
    std::string basis_from;
    std::string out_root;

    // Labels written to config.txt/summary.txt.
    std::string potential_label;
    std::string reference_units_label;
    std::string expanded_potential_label;
};

double recoil_energy_for_lattice(double lattice_a);

void set_lattice_depths_over_recoil(PumpConfig& cfg,
                                    double vs_depth_over_Er,
                                    double vl_depth_over_Er);

PumpConfig make_pump_config(const std::string& name);

std::vector<std::string> available_pump_configs();

}  // namespace pumpconfig
