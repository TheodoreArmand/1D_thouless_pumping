#include "nointeraction_src/pump_common.hpp"

#include "pumpconfig/pump_config.hpp"

int main() {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("lohse_n2_free");

    pump2::RunOptions opt;
    opt.base_terms = ecg1d::HamiltonianTerms::kinetic_only();
    opt.config_appendix =
        "driver=main_2ninteraction\n"
        "hamiltonian_kinetic=1\n"
        "hamiltonian_delta=0\n"
        "hamiltonian_gaussian=0\n";

    return pump2::run_pump(cfg, opt);
}
