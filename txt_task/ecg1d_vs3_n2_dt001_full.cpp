#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include "pumpconfig/pump_config.hpp"
#include "run_report.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

pumpconfig::PumpConfig make_vs3_n2_base() {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
    cfg.N = 2;
    cfg.K = 32;
    cfg.dt = 1.0e-2;
    cfg.manual_fine_dt_enabled = true;
    cfg.manual_fine_dt_factor = 0.1;
    cfg.manual_fine_dt_windows_s = {{0.2, 0.3}, {0.7, 0.8}};
    cfg.pump_period = 160.0 * pumpconfig::pi;
    cfg.total_time = 160.0 * pumpconfig::pi;
    cfg.trace_every = 25;
    cfg.snapshot_every = 25;
    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv";
    cfg.sigma_gauss = 1.0;
    return cfg;
}

std::string sigma_tag(double sigma) {
    return ecg1d::format_output_tag(sigma);
}

int run_free() {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base();
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K32_dt0p01_free";
    cfg.model = "legacy_prb_3_3_n2_pathpad_K32_dt0p01_free";
    cfg.out_root = "out/vs3_n2_dt0p01_T160pi_K32_sweep/free";
    cfg.potential_label = "legacy_prb_3_3_n2_free_cosine";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]";
    cfg.g_gauss_over_Er = 0.0;

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    opt.config_appendix =
        "driver=ecg1d_vs3_n2_dt001_full_free\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n"
        "dt_planning_estimate_seconds_per_step=0.35\n"
        "dt_planning_estimate_steps=140747\n"
        "dt_planning_estimate_wall_hours=13.7\n";
    return pump2::run_pump(cfg, opt);
}

int run_gauss(double sigma) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base();
    cfg.sigma_gauss = sigma;
    const std::string stag = sigma_tag(sigma);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K32_dt0p01_gauss_sigma" + stag;
    cfg.model = "legacy_prb_3_3_n2_pathpad_K32_dt0p01_gauss_sigma" + stag;
    cfg.out_root = "out/vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma" + stag;
    cfg.potential_label = "legacy_prb_3_3_n2_gaussian_pair";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    cfg.g_gauss_over_Er = 0.3;

    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    opt.config_appendix +=
        "driver=ecg1d_vs3_n2_dt001_full_gauss\n"
        "sweep_parameter=sigma_gauss\n"
        "sweep_sigma_gauss_code=" + std::to_string(sigma) + "\n"
        "dt_planning_estimate_seconds_per_step=0.80\n"
        "dt_planning_estimate_steps=140747\n"
        "dt_planning_estimate_wall_hours=31.3\n";
    return pump2::run_pump(cfg, opt);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) {
            std::cerr << "usage: " << argv[0] << " free|gauss [sigma_gauss]\n";
            return 2;
        }
        std::filesystem::create_directories("out/vs3_n2_dt0p01_T160pi_K32_sweep");
        const std::string mode = argv[1];
        if (mode == "free") return run_free();
        if (mode == "gauss") {
            const double sigma = (argc >= 3) ? std::stod(argv[2]) : 1.0;
            return run_gauss(sigma);
        }
        std::cerr << "unknown mode: " << mode << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "FATAL vs3 N2 dt=0.01 full run: " << e.what() << "\n";
        return 2;
    }
}
