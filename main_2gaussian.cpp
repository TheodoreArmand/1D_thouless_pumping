#include "2gaussian_src/gaussian_terms.hpp"
#include "nointeraction_src/pump_common.hpp"

#include "pumpconfig/pump_config.hpp"

int main() {
    pumpconfig::PumpConfig cfg = pumpconfig::make_pump_config("lohse_n2_gauss");
    pump2::RunOptions opt = pump2gaussian::make_gaussian_options(cfg);
    return pump2::run_pump(cfg, opt);
}
