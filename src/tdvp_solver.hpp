#pragma once

#include "basis_params.hpp"
#include "types.hpp"

#include <vector>

namespace ecg1d {

// TDVP parameter index: (a1, a2, a3, a4)
// a1=1: u, a1=2: B, a1=3: R, a1=4: A
struct AlphaIndex {
    int a1, a2, a3, a4;
};

// Configuration: which Hamiltonian terms to include.
struct HamiltonianTerms {
    bool kinetic  = true;
    bool delta    = true;
    bool gaussian = true;
    double one_body_constant = 0.0;  // constant c in c * sum_a 1
    std::vector<CosineTerm> cosine_terms;

    static HamiltonianTerms all() { return {}; }
    static HamiltonianTerms kinetic_only() {
        HamiltonianTerms t;
        t.delta = false;
        t.gaussian = false;
        return t;
    }
};

// Linear solve settings used by real-time TDVP.
struct SolverConfig {
    double lambda_C = 1e-8;  // Tikhonov regularization for C matrix
    double rcond    = 1e-4;  // Relative SVD truncation threshold
};

// Assemble the full C metric tensor matrix.
MatrixXcd assemble_C(const std::vector<AlphaIndex>& alpha_z_list,
                     const std::vector<BasisParams>& basis);

// Assemble the real-time Hamiltonian gradient vector.
VectorXcd assemble_grad(const std::vector<AlphaIndex>& alpha_z_list,
                        const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms = HamiltonianTerms::all());

// Update basis function parameters with a dz step.
void update_basis_function(std::vector<BasisParams>& basis,
                           const VectorXcd& dz,
                           double dt,
                           const std::vector<AlphaIndex>& alpha_z_list);

// Compute total energy <H>/<psi|psi>.
Cd compute_total_energy(const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms = HamiltonianTerms::all());

}  // namespace ecg1d
