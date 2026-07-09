#define make_legacy_prb_3_3_config make_legacy_prb_3_3_config_tmp100
#define main ecg1d_n1_vs3_tmp100_main
#include "/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/main.cpp"
#undef main
#undef make_legacy_prb_3_3_config

namespace pumpconfig {

PumpConfig make_legacy_prb_3_3_config();

PumpConfig make_legacy_prb_3_3_config_tmp100() {
    PumpConfig cfg = make_legacy_prb_3_3_config();
    cfg.config_name = "legacy_prb_3_3_n1_tmp100";
    cfg.model = "legacy_prb_tmp100";
    cfg.total_time = 100.0 * cfg.dt;
    cfg.out_root =
        "/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/out";
    return cfg;
}

}  // namespace pumpconfig

int main() {
    return ecg1d_n1_vs3_tmp100_main();
}
