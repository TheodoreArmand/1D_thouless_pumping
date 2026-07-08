#include "pumpconfig/lohse_n2_gauss.hpp"
#include "pumpconfig/lohse_n2_free.hpp"

namespace pumpconfig {

PumpConfig make_lohse_n2_gauss_config() {
    PumpConfig cfg = make_lohse_n2_free_config();
    cfg.config_name = "lohse_n2_gauss";
    cfg.model = "lohse_n2_gauss";
    cfg.g_gauss_over_Er = 0.3;
    cfg.sigma_gauss = 1.0;
    cfg.out_root = "out/n2_gauss_validation";
    cfg.potential_label = "lohse_n2_gaussian_pair";
    cfg.expanded_potential_label =
        "sum_a[-Vs*cos(4*pi*x_a/a)-Vl*cos(2*pi*x_a/a+phi)]"
        "+g_gauss*exp(-(x_0-x_1)^2/sigma_gauss^2)";
    return cfg;
}

}  // namespace pumpconfig
