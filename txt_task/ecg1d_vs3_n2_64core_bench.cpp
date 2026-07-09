#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include "pumpconfig/pump_config.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

namespace {

pumpconfig::PumpConfig make_vs3_n2_base(int steps) {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
    cfg.config_name = "legacy_prb_3_3_n2_pathpad";
    cfg.model = "legacy_prb_3_3_n2_pathpad";
    cfg.N = 2;
    cfg.K = 24;
    cfg.total_time = static_cast<double>(steps) * cfg.dt;
    cfg.trace_every = 1;
    cfg.snapshot_every = 0;
    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K24.csv";
    cfg.out_root = "out/bench_vs3_n2_64core_10step/free";
    cfg.potential_label = "legacy_prb_3_3_n2_free_cosine";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]";
    cfg.g_gauss_over_Er = 0.0;
    cfg.sigma_gauss = 1.0;
    return cfg;
}

int run_free(int steps) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(steps);
    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    opt.config_appendix =
        "driver=ecg1d_vs3_n2_64core_bench_free\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n";
    return pump2::run_pump(cfg, opt);
}

int run_gauss(int steps) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(steps);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_gauss";
    cfg.model = "legacy_prb_3_3_n2_pathpad_gauss";
    cfg.out_root = "out/bench_vs3_n2_64core_10step/gauss";
    cfg.potential_label = "legacy_prb_3_3_n2_gaussian_pair";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    cfg.g_gauss_over_Er = 0.3;

    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    opt.config_appendix += "driver=ecg1d_vs3_n2_64core_bench_gauss\n";
    return pump2::run_pump(cfg, opt);
}

}  // namespace

int main() {
    try {
        constexpr int steps = 10;
        std::filesystem::create_directories("out/bench_vs3_n2_64core_10step");
        const int free_code = run_free(steps);
        const int gauss_code = run_gauss(steps);
        return (free_code == 0 && gauss_code == 0) ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "FATAL vs3 N2 64-core bench: " << e.what() << "\n";
        return 2;
    }
}
