#include "pumpconfig/lohse_n2_free.hpp"

namespace pumpconfig {

PumpConfig make_lohse_n2_free_config() {
    PumpConfig cfg;
    cfg.config_name = "lohse_n2_free";
    cfg.model = "lohse_n2_free";

    cfg.N = 2;
    cfg.K = 24;

    cfg.lattice_a = 8.0;
    set_lattice_depths_over_recoil(cfg, 5.0, 2.5);

    cfg.pump_period = 500.0 * pi;
    cfg.total_time = 12.0;
    cfg.dt = 2.0e-3;

    cfg.phase_schedule_csv =
        "rice_mele_reference/lohes_experience/gap_adaptive_lohse_maincpp_schedule.csv";

    cfg.snapshot_every = 0;
    cfg.trace_every = 25;

    cfg.lambda_C = 1.0e-8;
    cfg.rcond = 1.0e-4;

    cfg.g_gauss_over_Er = 0.0;
    cfg.sigma_gauss = 1.0;

    cfg.basis_from = "initial_state/Vs10Ers_Vl5Ers/initial_lohse_N2_K24.csv";
    cfg.out_root = "out/n2_free_validation";

    cfg.potential_label = "lohse_n2_free_cosine";
    cfg.reference_units_label = "lohse_nphys2016_halfdepth_Ers";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]";

    return cfg;
}

}  // namespace pumpconfig
