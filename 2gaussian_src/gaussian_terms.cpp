#include "2gaussian_src/gaussian_terms.hpp"

#include "physical_constants.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace pump2gaussian {

pump2::RunOptions make_gaussian_options(const pumpconfig::PumpConfig& cfg) {
    if (cfg.N < 2) {
        throw std::runtime_error("Gaussian pair driver requires N >= 2");
    }
    if (!std::isfinite(cfg.g_gauss_over_Er)) {
        throw std::runtime_error("g_gauss_over_Er must be finite");
    }
    if (!(cfg.sigma_gauss > 0.0) || !std::isfinite(cfg.sigma_gauss)) {
        throw std::runtime_error("sigma_gauss must be positive and finite");
    }

    // The engine integrates V_int(x_a-x_b) =
    // g_gauss * exp(-(x_a-x_b)^2 / sigma_gauss^2). For Gaussian basis pairs,
    // compute_H_Mijab gives sigma/sqrt(sigma^2+h) * exp(-p^2/(sigma^2+h)),
    // the analytic convolution of that unnormalized potential with the
    // relative-coordinate Gaussian of variance h/2.
    ecg1d::g_gauss = cfg.g_gauss_over_Er * cfg.recoil_energy;
    ecg1d::sigma_gauss = cfg.sigma_gauss;

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    opt.base_terms.gaussian = true;
    opt.evolve_A = true;
    opt.config_appendix =
        "driver=main_2gaussian\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=1\n"
        "tdvp_evolve_A=1\n";

    std::ostringstream oss;
    oss << "g_gauss_code=" << ecg1d::g_gauss << "\n"
        << "sigma_gauss_code_runtime=" << ecg1d::sigma_gauss << "\n"
        << "V_int_convention=g*exp(-(x0-x1)^2/sigma^2)\n";
    opt.config_appendix += oss.str();
    return opt;
}

}  // namespace pump2gaussian
