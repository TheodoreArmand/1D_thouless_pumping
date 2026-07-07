#include "pumpconfig/legacy_prb_3_3.hpp"

namespace pumpconfig {
namespace {

constexpr double legacy_prb_mass = 1.0;
constexpr double lattice_depth_over_Er = 3.0;

// Experimental wavelengths kept as run metadata. The simulation itself uses
// code units set by lattice_a, hbar, and mass_code_units.
constexpr double lambda_l_nm = 1534.0;
constexpr double lambda_s_nm = 767.0;

// Reference gap markers from rice_mele_reference/rice_mele_reference_data.npz.
// The adaptive clock uses these only to choose the numerical time grid; the
// physical pump phase remains phi(t)=2*pi*t/T.
constexpr double reference_g12_phi0 = 0.11129214004023469;
constexpr double reference_g12_phipi = 0.11129214004023469;
constexpr double reference_g12_max = 4.515385395187394;

}  // namespace

PumpConfig make_legacy_prb_3_3_config() {
    PumpConfig cfg;
    cfg.config_name = "legacy_prb_3_3";
    cfg.model = "legacy_prb";

    cfg.N = 1;
    cfg.K = 16;
    cfg.mass_code_units = legacy_prb_mass;
    cfg.lattice_a = 8.0;
    set_lattice_depths_over_recoil(cfg, lattice_depth_over_Er, lattice_depth_over_Er);
    cfg.omega = 0.0;

    cfg.pump_period = 160.0 * pi;
    cfg.total_time = 160.0 * pi;
    cfg.dt = 2.0e-3;
    cfg.dynamic_dt = false;
    cfg.dt_min = 1.0;
    cfg.dt_max = 50.0;
    cfg.dt_gap_power = 2.0;
    cfg.integrator = "rk4";
    cfg.adaptive_abs_tol = 1.0e-7;
    cfg.adaptive_rel_tol = 1.0e-4;
    cfg.adaptive_dt_min = 1.0e-6;
    cfg.adaptive_max_reject = 20;

    cfg.phase_schedule = "local_csv";
    cfg.phase_schedule_csv = "rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_vs3vl3_maincpp_schedule.csv";
    cfg.dphi_max = 0.0;

    cfg.snapshot_every = 25;
    cfg.trace_every = 25;
    cfg.lambda_C = 1.0e-8;
    cfg.rcond = 1.0e-4;
    cfg.u_split_trotter = false;
    cfg.enforce_norm = true;
    cfg.autoinit = false;

    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N1_K16.csv";
    cfg.out_root = "out/pump_vs3pad_gapadaptive_T160pi";

    cfg.potential_label = "legacy_prb_cosine";
    cfg.reference_units_label = "legacy_prb_hbar1_mass1";
    cfg.expanded_potential_label = "-Vs*cos(4*pi*x/a)-Vl*cos(2*pi*x/a+phi)";
    cfg.lambda_l_nm = lambda_l_nm;
    cfg.lambda_s_nm = lambda_s_nm;
    cfg.reference_gap12 = {
        reference_g12_phi0,
        reference_g12_phipi,
        reference_g12_max,
        true
    };

    return cfg;
}

}  // namespace pumpconfig
