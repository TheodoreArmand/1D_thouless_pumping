#include "pumpconfig/pump_config.hpp"

#include "physical_constants.hpp"
#include "pumpconfig/legacy_prb_3_3.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace pumpconfig {
namespace {

std::string lower_copy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::string normalized_config_name(std::string name) {
    name = lower_copy(name);
    if (name == "legacy" || name == "prb" || name == "legacy_prb") {
        return "legacy_prb_3_3";
    }
    if (name == "3_3" || name == "vs3_vl3" || name == "vs3er_vl3er") {
        return "legacy_prb_3_3";
    }
    return name;
}

}  // namespace

double recoil_energy_for_lattice(double lattice_a, double mass_code_units) {
    if (!(lattice_a > 0.0) || !std::isfinite(lattice_a)) {
        throw std::runtime_error("lattice_a must be positive and finite");
    }
    if (!(mass_code_units > 0.0) || !std::isfinite(mass_code_units)) {
        throw std::runtime_error("mass_code_units must be positive and finite");
    }
    const double k_s = 2.0 * pi / lattice_a;
    return ecg1d::hbar * ecg1d::hbar * k_s * k_s / (2.0 * mass_code_units);
}

void set_lattice_depths_over_recoil(PumpConfig& cfg,
                                    double vs_depth_over_Er,
                                    double vl_depth_over_Er) {
    if (!(vs_depth_over_Er > 0.0) || !std::isfinite(vs_depth_over_Er)) {
        throw std::runtime_error("Vs depth over Er must be positive and finite");
    }
    if (!(vl_depth_over_Er > 0.0) || !std::isfinite(vl_depth_over_Er)) {
        throw std::runtime_error("Vl depth over Er must be positive and finite");
    }
    cfg.recoil_energy = recoil_energy_for_lattice(cfg.lattice_a, cfg.mass_code_units);
    cfg.Vs_depth_over_Er = vs_depth_over_Er;
    cfg.Vl_depth_over_Er = vl_depth_over_Er;
    cfg.Vs = cfg.Vs_depth_over_Er * cfg.recoil_energy;
    cfg.Vl = cfg.Vl_depth_over_Er * cfg.recoil_energy;
}

PumpConfig make_pump_config(const std::string& name) {
    const std::string normalized = normalized_config_name(name);
    if (normalized == "legacy_prb_3_3") {
        return make_legacy_prb_3_3_config();
    }

    std::ostringstream oss;
    oss << "unknown pump config '" << name << "'. Available configs:";
    for (const auto& config_name : available_pump_configs()) {
        oss << " " << config_name;
    }
    throw std::runtime_error(oss.str());
}

std::vector<std::string> available_pump_configs() {
    return {"legacy_prb_3_3"};
}

}  // namespace pumpconfig
