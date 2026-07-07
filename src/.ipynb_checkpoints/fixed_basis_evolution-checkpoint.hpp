#pragma once

#include "basis_params.hpp"
#include "types.hpp"

#include <vector>

namespace ecg1d {

// Free evolution with fixed Gaussian shapes: solve i S du/dt = H u exactly.
// Only the linear coefficients u change; A, B, and R stay fixed.
void free_evolve_fixed_basis(std::vector<BasisParams>& basis,
                             const MatrixXcd& H,
                             const MatrixXcd& S,
                             double duration);

} // namespace ecg1d
