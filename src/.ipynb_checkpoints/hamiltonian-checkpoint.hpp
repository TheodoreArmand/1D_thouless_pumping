#pragma once
#include "basis_params.hpp"
#include "permutation.hpp"
#include "pair_cache.hpp"
#include <vector>

namespace ecg1d {

// Overlap: S = sum_{i,j} conj(u_i) u_j sum_p sign(p) M_G(i,j,p)
Cd overlap(const std::vector<BasisParams>& basis);

// Kinetic energy: <T> = (-hbar^2/2m) sum_{i,j} conj(u_i) u_j sum_p sign(p) M_G*P_Mij
Cd kinetic_energy_functional(const std::vector<BasisParams>& basis);

// Delta contact interaction: <V_delta> = g * sum_{i,j,p} ... M_G * G_Mij
Cd Delta_contact_functional(const std::vector<BasisParams>& basis);

// Gaussian interaction: <V_gauss> = g_gauss * sum_{i,j,p} ... M_G * H_Mij
Cd Gaussian_interaction_functional(const std::vector<BasisParams>& basis);

// General cosine potential: sum_r C_r sum_a cos(A_r x_a + b_r)
Cd general_cosine_functional(const std::vector<BasisParams>& basis,
                             const std::vector<CosineTerm>& terms);

} // namespace ecg1d
