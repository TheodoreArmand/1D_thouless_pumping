#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include <vector>

namespace ecg1d {

// Hamiltonian partial derivatives: d(<H>/S)/dz_alpha
// result = (1/S)*d<H>/dz - (<H>/S²)*dS/dz
// where d<H>/dz = d<H>_dMG + d<H>_dKernel
struct HamiltonianGradientContext {
    Cd S = Cd(1.0, 0.0);
    Cd kinetic_value = Cd(0.0, 0.0);
    Cd delta_value = Cd(0.0, 0.0);
    Cd gaussian_value = Cd(0.0, 0.0);
    Cd cosine_value = Cd(0.0, 0.0);
};

Cd calculate_Hamiltonian_kinetic_partial(int a1, int a2, int a3, int a4,
                                          bool Real,
                                          const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_kinetic_partial(int a1, int a2, int a3, int a4,
                                          bool Real,
                                          const std::vector<BasisParams>& basis,
                                          const HamiltonianGradientContext& ctx,
                                          Cd dS);

Cd calculate_Hamiltonian_delta_partial(int a1, int a2, int a3, int a4,
                                        bool Real,
                                        const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_delta_partial(int a1, int a2, int a3, int a4,
                                        bool Real,
                                        const std::vector<BasisParams>& basis,
                                        const HamiltonianGradientContext& ctx,
                                        Cd dS);

Cd calculate_Hamiltonian_gaussian_partial(int a1, int a2, int a3, int a4,
                                           bool Real,
                                           const std::vector<BasisParams>& basis);

Cd calculate_Hamiltonian_gaussian_partial(int a1, int a2, int a3, int a4,
                                           bool Real,
                                           const std::vector<BasisParams>& basis,
                                           const HamiltonianGradientContext& ctx,
                                           Cd dS);

Cd calculate_Hamiltonian_general_cosine_partial(
    int a1, int a2, int a3, int a4,
    bool Real,
    const std::vector<BasisParams>& basis,
    const std::vector<CosineTerm>& terms);

Cd calculate_Hamiltonian_general_cosine_partial(
    int a1, int a2, int a3, int a4,
    bool Real,
    const std::vector<BasisParams>& basis,
    const std::vector<CosineTerm>& terms,
    const HamiltonianGradientContext& ctx,
    Cd dS);

} // namespace ecg1d
