#pragma once

#include "nointeraction_src/pump_common.hpp"
#include "pumpconfig/pump_config.hpp"

namespace pump2gaussian {

// Build the N=2 Gaussian-interaction run options and set the global Gaussian
// constants once before TDVP assembly starts.
pump2::RunOptions make_gaussian_options(const pumpconfig::PumpConfig& cfg);

}  // namespace pump2gaussian
