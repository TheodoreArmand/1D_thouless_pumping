#include "tdvp_solver.hpp"

#include "derivatives.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"

namespace ecg1d {

MatrixXcd assemble_C(const std::vector<AlphaIndex>& alpha_z_list,
                     const std::vector<BasisParams>& basis) {
    const int d = static_cast<int>(alpha_z_list.size());
    MatrixXcd C_mat = MatrixXcd::Zero(d, d);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int a = 0; a < d; a++) {
        for (int b = 0; b < d; b++) {
            const auto& alpha = alpha_z_list[a];
            const auto& beta = alpha_z_list[b];
            C_mat(a, b) = calculate_C(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                       beta.a1, beta.a2, beta.a3, beta.a4,
                                       basis);
        }
    }
    return C_mat;
}

namespace {

Cd grad_H_for_alpha(const AlphaIndex& alpha,
                    const std::vector<BasisParams>& basis,
                    const HamiltonianTerms& terms) {
    const int N = basis[0].N();
    Cd result(0, 0);

    if (terms.kinetic) {
        result += calculate_Hamiltonian_kinetic_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                        /*Real=*/false, basis);
    }
    if (terms.delta && N >= 2) {
        result += calculate_Hamiltonian_delta_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                      /*Real=*/false, basis);
    }
    if (terms.gaussian && N >= 2) {
        result += calculate_Hamiltonian_gaussian_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                         /*Real=*/false, basis);
    }
    if (!terms.cosine_terms.empty()) {
        result += calculate_Hamiltonian_general_cosine_partial(
            alpha.a1, alpha.a2, alpha.a3, alpha.a4, /*Real=*/false, basis, terms.cosine_terms);
    }

    return result;
}

}  // namespace

VectorXcd assemble_grad(const std::vector<AlphaIndex>& alpha_z_list,
                        const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms) {
    const int d = static_cast<int>(alpha_z_list.size());
    VectorXcd g = VectorXcd::Zero(d);

    #pragma omp parallel for schedule(static)
    for (int a = 0; a < d; a++) {
        g(a) = grad_H_for_alpha(alpha_z_list[a], basis, terms);
    }
    return g;
}

void update_basis_function(std::vector<BasisParams>& basis,
                           const VectorXcd& dz,
                           double dt,
                           const std::vector<AlphaIndex>& alpha_z_list) {
    const int length = static_cast<int>(alpha_z_list.size());
    for (int i = 0; i < length; i++) {
        if (i >= dz.size()) continue;

        const auto& idx = alpha_z_list[i];
        if (idx.a1 == 1) {
            basis[idx.a2].u += dz(i) * dt;
        } else if (idx.a1 == 2) {
            basis[idx.a2].B(idx.a3, idx.a3) += dz(i) * dt;
        } else if (idx.a1 == 3) {
            basis[idx.a2].R(idx.a3) += dz(i) * dt;
        } else if (idx.a1 == 4) {
            const Cd val = basis[idx.a2].A(idx.a3, idx.a4) + dz(i) * dt;
            basis[idx.a2].A(idx.a3, idx.a4) = val;
            if (idx.a3 != idx.a4) {
                basis[idx.a2].A(idx.a4, idx.a3) = val;
            }
        }
    }
}

Cd compute_total_energy(const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms) {
    const int N = basis[0].N();
    const Cd S = overlap(basis);
    Cd E_num(0, 0);

    if (terms.kinetic) {
        E_num += kinetic_energy_functional(basis);
    }
    if (terms.delta && N >= 2) {
        E_num += Delta_contact_functional(basis);
    }
    if (terms.gaussian && N >= 2) {
        E_num += Gaussian_interaction_functional(basis);
    }
    if (!terms.cosine_terms.empty()) {
        E_num += general_cosine_functional(basis, terms.cosine_terms);
    }
    if (terms.one_body_constant != 0.0) {
        E_num += terms.one_body_constant * N * S;
    }

    return E_num / S;
}

}  // namespace ecg1d
