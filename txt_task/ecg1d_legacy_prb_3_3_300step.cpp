#define make_legacy_prb_3_3_config make_legacy_prb_3_3_config_300step
#define main ecg1d_legacy_prb_3_3_locked_main
#include "../main.cpp"
#undef main
#undef make_legacy_prb_3_3_config

namespace pumpconfig {

PumpConfig make_legacy_prb_3_3_config();

PumpConfig make_legacy_prb_3_3_config_300step() {
    PumpConfig cfg = make_legacy_prb_3_3_config();
    cfg.config_name = "legacy_prb_3_3_300step_8thread_bench";
    cfg.model = "legacy_prb_300step_bench";
    cfg.total_time = 300.0 * cfg.dt;
    cfg.out_root = "out/bench_legacy_prb_3_3_300step_8thread";
    return cfg;
}

}  // namespace pumpconfig

int main() {
    return ecg1d_legacy_prb_3_3_locked_main();
}
