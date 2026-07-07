#include "pumpconfig/pump_config.hpp"

#include "physical_constants.hpp"
#include "pumpconfig/legacy_prb_3_3.hpp"

#include <cmath>
#include <stdexcept>

namespace pumpconfig {
double recoil_energy_for_lattice(double lattice_a) {
    if (!(lattice_a > 0.0) || !std::isfinite(lattice_a)) {
        throw std::runtime_error("lattice_a must be positive and finite");
    }
    if (!(ecg1d::mass > 0.0) || !std::isfinite(ecg1d::mass)) {
        throw std::runtime_error("ecg1d::mass must be positive and finite");
    }
    const double k_s = 2.0 * pi / lattice_a;
    return ecg1d::hbar * ecg1d::hbar * k_s * k_s / (2.0 * ecg1d::mass);
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
    cfg.recoil_energy = recoil_energy_for_lattice(cfg.lattice_a);
    cfg.Vs_depth_over_Er = vs_depth_over_Er;
    cfg.Vl_depth_over_Er = vl_depth_over_Er;
    cfg.Vs = cfg.Vs_depth_over_Er * cfg.recoil_energy;
    cfg.Vl = cfg.Vl_depth_over_Er * cfg.recoil_energy;
}

PumpConfig make_pump_config(const std::string& name) {
    if (name == "legacy_prb_3_3") {
        return make_legacy_prb_3_3_config();
    }

    throw std::runtime_error(
        "unknown pump config '" + name + "'. Available configs: legacy_prb_3_3");
}

std::vector<std::string> available_pump_configs() {
    return {"legacy_prb_3_3"};
}

}  // namespace pumpconfig
