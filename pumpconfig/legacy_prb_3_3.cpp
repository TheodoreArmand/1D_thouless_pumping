#include "pumpconfig/legacy_prb_3_3.hpp"

namespace pumpconfig {

PumpConfig make_legacy_prb_3_3_config() {
    PumpConfig cfg;
    cfg.config_name = "legacy_prb_3_3";
    cfg.model = "legacy_prb";

    cfg.N = 1;
    cfg.K = 16;

    // Lattice and energy scale.
    cfg.lattice_a = 8.0;
    set_lattice_depths_over_recoil(cfg, 3.0, 3.0);

    // Time grid.
    cfg.pump_period = 160.0 * pi;
    cfg.total_time = 160.0 * pi;
    cfg.dt = 2.0e-3;

    // Phase schedule.
    cfg.phase_schedule_csv =
        "rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_vs3vl3_maincpp_schedule.csv";

    // Output cadence.
    cfg.snapshot_every = 25;
    cfg.trace_every = 25;

    // TDVP linear solve.
    cfg.lambda_C = 1.0e-8;
    cfg.rcond = 1.0e-4;

    // Input/output paths.
    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N1_K16.csv";
    cfg.out_root = "out/pump_vs3pad_gapadaptive_T160pi";

    // Labels written to config.txt/summary.txt.
    cfg.potential_label = "legacy_prb_cosine";
    cfg.reference_units_label = "legacy_prb_hbar1_mass1";
    cfg.expanded_potential_label = "-Vs*cos(4*pi*x/a)-Vl*cos(2*pi*x/a+phi)";

    return cfg;
}

}  // namespace pumpconfig
