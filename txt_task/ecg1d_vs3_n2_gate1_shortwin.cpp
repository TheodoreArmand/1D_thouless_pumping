// Gate-1 short-window driver for the R1 / R2 / R3a experiments of
// vs3_n2_k32_success_roadmap_report.html.
//
// One driver covers the cheapest discriminators for the N=2 K=32 failure:
//   free  : kinetic-only, A frozen (160-dim control).
//   gauss : Gaussian pair + A FROZEN (R1 / E1 — the never-run frozen-A gauss).
//   a01   : Gaussian pair + A_01-only (R2 — clean de-gauged 192-dim gauss).
// Each mode can be run with or without the R3a occupancy gate (4th arg):
//   gate   : occupancy gate ON  + occupancy_gate.csv diagnostics.
//   nogate : occupancy gate OFF + occupancy_gate.csv diagnostics (baseline H1
//            fingerprint — shows which block x tier owns the discarded RHS mass).
//
// Short window: total_time = 0.35 * pump_period with pump_period = 160*pi kept
// intact so phi(s) is the true full-cycle schedule; the run stops at s = 0.35
// (phi/pi ~ 0.523472 for the full-depth schedule), past the first pump
// transition. dt = 0.01 with a dt/10 fine window over s in [0.18, 0.32]
// (phi/pi ~ [0.484783, 0.515217]), centered on the minimum gap. The historical
// error starts near phi/pi ~ 0.42, but full20 already proved a coarse step is
// not the cause.
//
// Thresholds u_on / u_off default to 0.03 / 0.02 (roadmap suggests 0.02-0.05)
// and are overridable via ECG_OCC_U_ON / ECG_OCC_U_OFF.
//
// CLI: rcond1em4|rcond1em5  free|gauss|a01  [sigma_gauss]  [gate|nogate]

#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include "pumpconfig/pump_config.hpp"
#include "run_report.hpp"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

struct RcondCase {
    std::string tag;
    double value = 1.0e-4;
};

RcondCase parse_rcond_case(const std::string& arg) {
    if (arg == "rcond1em4" || arg == "1e-4" || arg == "1.0e-4") {
        return {"rcond1em4", 1.0e-4};
    }
    if (arg == "rcond1em5" || arg == "1e-5" || arg == "1.0e-5") {
        return {"rcond1em5", 1.0e-5};
    }
    throw std::runtime_error("unknown rcond case: " + arg);
}

double env_double(const char* name, double fallback) {
    const char* v = std::getenv(name);
    return v ? std::stod(v) : fallback;
}

pumpconfig::PumpConfig make_vs3_n2_base(const RcondCase& rc) {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
    cfg.N = 2;
    cfg.K = 32;
    cfg.dt = 1.0e-2;
    cfg.rcond = rc.value;
    cfg.manual_fine_dt_enabled = true;
    cfg.manual_fine_dt_factor = 0.1;
    cfg.manual_fine_dt_windows_s = {{0.18, 0.32}};
    cfg.pump_period = 160.0 * pumpconfig::pi;   // true full-cycle schedule
    cfg.total_time = 0.35 * 160.0 * pumpconfig::pi;  // stop at s = 0.35
    cfg.trace_every = 250;
    cfg.snapshot_every = 250;
    cfg.basis_from = "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv";
    cfg.sigma_gauss = 1.0;
    return cfg;
}

// Occupancy-gate switch shared by all modes. `gate` turns the R3a gate on;
// both settings still emit occupancy_gate.csv (diagnostics always on here).
struct GateCase {
    std::string tag;      // "gate" | "nogate"
    bool gate_on = false;
};

GateCase parse_gate_case(const std::string& arg) {
    if (arg == "gate")   return {"gate", true};
    if (arg == "nogate") return {"nogate", false};
    throw std::runtime_error("unknown gate case (expect gate|nogate): " + arg);
}

void apply_occupancy(pump2::RunOptions& opt, const GateCase& gc) {
    opt.occupancy_gate_enabled = gc.gate_on;
    opt.occupancy_diag_enabled = true;  // always emit the decomposition
    opt.occupancy_u_on = env_double("ECG_OCC_U_ON", 0.03);
    opt.occupancy_u_off = env_double("ECG_OCC_U_OFF", 0.02);
    opt.config_appendix +=
        "occupancy_gate_case=" + gc.tag + "\n"
        "occupancy_gate_enabled=" + std::to_string(gc.gate_on ? 1 : 0) + "\n"
        "occupancy_u_on=" + std::to_string(opt.occupancy_u_on) + "\n"
        "occupancy_u_off=" + std::to_string(opt.occupancy_u_off) + "\n";
}

std::string sigma_tag(double sigma) { return ecg1d::format_output_tag(sigma); }

std::string out_root_base(const RcondCase& rc) {
    return "out/vs3_n2_gate1_shortwin_" + rc.tag;
}

void apply_smoke_override(pumpconfig::PumpConfig& cfg) {
    const char* smoke_t = std::getenv("ECG_SMOKE_TOTAL_TIME");
    if (!smoke_t) return;
    cfg.total_time = std::stod(smoke_t);
    cfg.trace_every = 1;
    cfg.snapshot_every = 1;
    cfg.out_root += "_smoke";
}

const char* const kWindowProvenance =
    "gate1_window_s=0.0-0.35\n"
    "gate1_window_phi_over_pi=0.0-0.523472\n"
    "manual_fine_dt_factor=0.1\n"
    "manual_fine_dt_windows_s=0.18-0.32\n"
    "manual_fine_dt_windows_phi_over_pi=0.484783-0.515217\n"
    "dt_planning_estimate_seconds_per_step=0.55\n"
    "dt_planning_estimate_steps=80929\n"
    "dt_planning_estimate_wall_hours=12.4\n";

int run_free(const RcondCase& rc, const GateCase& gc) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(rc);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K32_gate1_" + rc.tag +
                      "_" + gc.tag + "_free";
    cfg.model = cfg.config_name;
    cfg.out_root = out_root_base(rc) + "/free_" + gc.tag;
    cfg.potential_label = "legacy_prb_3_3_n2_free_cosine";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]";
    cfg.g_gauss_over_Er = 0.0;
    apply_smoke_override(cfg);

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();  // A frozen (160-dim)
    opt.config_appendix =
        "driver=ecg1d_vs3_n2_gate1_shortwin_free\n"
        "gate1_mode=free_frozenA\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n"
        "tdvp_evolve_A=0\n"
        "basis_variant=pathpad_N2_K32\n"
        "rcond_case=" + rc.tag + "\n";
    opt.config_appendix += kWindowProvenance;
    apply_occupancy(opt, gc);
    return pump2::run_pump(cfg, opt);
}

int run_gauss_frozenA(const RcondCase& rc, double sigma, const GateCase& gc) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(rc);
    cfg.sigma_gauss = sigma;
    const std::string stag = sigma_tag(sigma);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K32_gate1_" + rc.tag +
                      "_" + gc.tag + "_gaussFrozenA_sigma" + stag;
    cfg.model = cfg.config_name;
    cfg.out_root = out_root_base(rc) + "/gaussFrozenA_" + gc.tag + "_sigma" + stag;
    cfg.potential_label = "legacy_prb_3_3_n2_gaussian_pair_frozenA";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    cfg.g_gauss_over_Er = 0.3;
    apply_smoke_override(cfg);

    // make_gaussian_options sets the g_gauss/sigma globals + gaussian term and
    // defaults A_01-only; R1/E1 wants A FROZEN, so override the evolve_A flags.
    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    opt.evolve_A = false;
    opt.evolve_A_offdiag_only = false;
    opt.config_appendix +=
        "driver=ecg1d_vs3_n2_gate1_shortwin_gaussFrozenA\n"
        "gate1_mode=gauss_frozenA\n"
        "gate1_tdvp_evolve_A_OVERRIDE=0\n"   // supersedes make_gaussian_options
        "basis_variant=pathpad_N2_K32\n"
        "sweep_parameter=sigma_gauss\n"
        "sweep_sigma_gauss_code=" + std::to_string(sigma) + "\n"
        "rcond_case=" + rc.tag + "\n";
    opt.config_appendix += kWindowProvenance;
    apply_occupancy(opt, gc);
    return pump2::run_pump(cfg, opt);
}

int run_gauss_a01(const RcondCase& rc, double sigma, const GateCase& gc) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(rc);
    cfg.sigma_gauss = sigma;
    const std::string stag = sigma_tag(sigma);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K32_gate1_" + rc.tag +
                      "_" + gc.tag + "_gaussA01_sigma" + stag;
    cfg.model = cfg.config_name;
    cfg.out_root = out_root_base(rc) + "/gaussA01_" + gc.tag + "_sigma" + stag;
    cfg.potential_label = "legacy_prb_3_3_n2_gaussian_pair_a01";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    cfg.g_gauss_over_Er = 0.3;
    apply_smoke_override(cfg);

    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);  // A_01-only
    opt.config_appendix +=
        "driver=ecg1d_vs3_n2_gate1_shortwin_gaussA01\n"
        "gate1_mode=gauss_a01_only\n"
        "basis_variant=pathpad_N2_K32\n"
        "sweep_parameter=sigma_gauss\n"
        "sweep_sigma_gauss_code=" + std::to_string(sigma) + "\n"
        "rcond_case=" + rc.tag + "\n";
    opt.config_appendix += kWindowProvenance;
    apply_occupancy(opt, gc);
    return pump2::run_pump(cfg, opt);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 5) {
            std::cerr << "usage: " << argv[0]
                      << " rcond1em4|rcond1em5 free|gauss|a01 [sigma_gauss] [gate|nogate]\n";
            return 2;
        }
        const RcondCase rc = parse_rcond_case(argv[1]);
        const std::string mode = argv[2];

        double sigma = 1.0;
        std::string gate_arg = "nogate";
        // Optional args are [sigma] then [gate|nogate], but tolerate the gate
        // token appearing in the sigma slot for the free mode's convenience.
        if (argc == 4) {
            const std::string a3 = argv[3];
            if (a3 == "gate" || a3 == "nogate") gate_arg = a3;
            else sigma = std::stod(a3);
        } else if (argc == 5) {
            sigma = std::stod(argv[3]);
            gate_arg = argv[4];
        }
        const GateCase gc = parse_gate_case(gate_arg);

        std::filesystem::create_directories(out_root_base(rc));
        if (mode == "free")  return run_free(rc, gc);
        if (mode == "gauss") return run_gauss_frozenA(rc, sigma, gc);
        if (mode == "a01")   return run_gauss_a01(rc, sigma, gc);
        std::cerr << "unknown mode: " << mode << " (expect free|gauss|a01)\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "FATAL vs3 N2 gate1 shortwin run: " << e.what() << "\n";
        return 2;
    }
}
