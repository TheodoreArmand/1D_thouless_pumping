#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include <vector>

namespace ecg1d {

// Hamiltonian partial derivatives: d(<H>/S)/dz_alpha
// result = (1/S)*d<H>/dz - (<H>/S²)*dS/dz
// where d<H>/dz = d<H>_dMG + d<H>_dKernel

Cd calculate_Hamiltonian_kinetic_partial(int a1, int a2, int a3, int a4,
                                          bool Real,
                                          const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_delta_partial(int a1, int a2, int a3, int a4,
                                        bool Real,
                                        const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_gaussian_partial(int a1, int a2, int a3, int a4,
                                           bool Real,
                                           const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_general_cosine_partial(
    int a1, int a2, int a3, int a4,
    bool Real,
    const std::vector<BasisParams>& basis,
    const std::vector<CosineTerm>& terms);

} // namespace ecg1d
