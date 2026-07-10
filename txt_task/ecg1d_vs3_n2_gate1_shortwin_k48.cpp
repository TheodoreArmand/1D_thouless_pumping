// K48 union-basis Gate-1 short-window driver (R5 of
// vs3_n2_k32_success_roadmap_report.html / txt_task/r5_basis_union_proposal.md).
//
// This is a copy of ecg1d_vs3_n2_gate1_shortwin.cpp with the ONLY knob changed
// being the basis: K = 48 with the width-enriched rich-union initial state
//   initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K48_richunion.csv
// (32 forced K32 path-pad terms + 16 greedy-selected reserves; A = 0 throughout;
//  static u-only span infidelity 0.97 -> 0.045 all-phi, 0.67 -> 0.032 Gate-1;
//  selection report: tmp/r5_basis_selection/R5_RESULT.md). Everything else --
// window (s in [0,0.35]), fine dt/10 in s in [0.18,0.32], dt=0.01, full-cycle
// schedule, occupancy diagnostics, mode set -- matches the K32 driver so the
// K48 basis is the single variable (roadmap P4).
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

// K48 rich-union basis. cfg.K MUST equal the CSV row count / 12 (loader checks).
constexpr int kK48 = 48;
const char* const kBasisK48 =
    "initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K48_richunion.csv";

pumpconfig::PumpConfig make_vs3_n2_base(const RcondCase& rc) {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("legacy_prb_3_3");
    cfg.N = 2;
    cfg.K = kK48;                                     // K48 (was 32)
    cfg.dt = 1.0e-2;
    cfg.rcond = rc.value;
    cfg.manual_fine_dt_enabled = true;
    cfg.manual_fine_dt_factor = 0.1;
    cfg.manual_fine_dt_windows_s = {{0.18, 0.32}};
    cfg.pump_period = 160.0 * pumpconfig::pi;         // true full-cycle schedule
    cfg.total_time = 0.35 * 160.0 * pumpconfig::pi;   // stop at s = 0.35
    cfg.trace_every = 250;
    cfg.snapshot_every = 250;
    cfg.basis_from = kBasisK48;                       // K48 rich-union (was K32)
    cfg.sigma_gauss = 1.0;
    return cfg;
}

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
    return "out/vs3_n2_gate1_shortwin_k48_" + rc.tag;
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
    "basis_variant=K48_richunion\n"
    "selection_report=tmp/r5_basis_selection/R5_RESULT.md\n"
    "pool_tiers=P-A+rich_width_enriched\n"
    "u_floor=1e-4max\n"
    // K48 ~2.25x pair table + param_dim 240 (frozen-A) vs K32 160; ~1.0-1.3 s/step.
    "dt_planning_estimate_seconds_per_step=1.2\n"
    "dt_planning_estimate_steps=80929\n"
    "dt_planning_estimate_wall_hours=27.0\n";

int run_free(const RcondCase& rc, const GateCase& gc) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(rc);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K48_gate1_" + rc.tag +
                      "_" + gc.tag + "_free";
    cfg.model = cfg.config_name;
    cfg.out_root = out_root_base(rc) + "/free_" + gc.tag;
    cfg.potential_label = "legacy_prb_3_3_n2_free_cosine";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]";
    cfg.g_gauss_over_Er = 0.0;
    apply_smoke_override(cfg);

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();  // A frozen (240-dim)
    opt.config_appendix =
        "driver=ecg1d_vs3_n2_gate1_shortwin_k48_free\n"
        "gate1_mode=free_frozenA\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n"
        "tdvp_evolve_A=0\n"
        "basis_variant=K48_richunion\n"
        "rcond_case=" + rc.tag + "\n";
    opt.config_appendix += kWindowProvenance;
    apply_occupancy(opt, gc);
    return pump2::run_pump(cfg, opt);
}

int run_gauss_frozenA(const RcondCase& rc, double sigma, const GateCase& gc) {
    pumpconfig::PumpConfig cfg = make_vs3_n2_base(rc);
    cfg.sigma_gauss = sigma;
    const std::string stag = sigma_tag(sigma);
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K48_gate1_" + rc.tag +
                      "_" + gc.tag + "_gaussFrozenA_sigma" + stag;
    cfg.model = cfg.config_name;
    cfg.out_root = out_root_base(rc) + "/gaussFrozenA_" + gc.tag + "_sigma" + stag;
    cfg.potential_label = "legacy_prb_3_3_n2_gaussian_pair_frozenA";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    cfg.g_gauss_over_Er = 0.3;
    apply_smoke_override(cfg);

    // make_gaussian_options defaults A_01-only; R1/E1 wants A FROZEN.
    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    opt.evolve_A = false;
    opt.evolve_A_offdiag_only = false;
    opt.config_appendix +=
        "driver=ecg1d_vs3_n2_gate1_shortwin_k48_gaussFrozenA\n"
        "gate1_mode=gauss_frozenA\n"
        "gate1_tdvp_evolve_A_OVERRIDE=0\n"
        "basis_variant=K48_richunion\n"
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
    cfg.config_name = "legacy_prb_3_3_n2_pathpad_K48_gate1_" + rc.tag +
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
        "driver=ecg1d_vs3_n2_gate1_shortwin_k48_gaussA01\n"
        "gate1_mode=gauss_a01_only\n"
        "basis_variant=K48_richunion\n"
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
        std::cerr << "FATAL vs3 N2 gate1 shortwin K48 run: " << e.what() << "\n";
        return 2;
    }
}
