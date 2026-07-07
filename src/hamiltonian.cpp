#include "hamiltonian.hpp"
#include "interaction_kernels.hpp"
#include "physical_constants.hpp"

namespace ecg1d {

Cd overlap(const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    Cd result(0.0, 0.0);

    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            result += con_ui * basis[j].u * sum_p;
        }
    }
    return result;
}

Cd kinetic_energy_functional(const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    Cd result(0.0, 0.0);

    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G * compute_P_Mij(c);
            }
            result += con_ui * basis[j].u * sum_p;
        }
    }
    return result * (-hbar * hbar) / (2.0 * mass);
}

Cd Delta_contact_functional(const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    Cd result(0.0, 0.0);

    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G * compute_G_Mij(c);
            }
            result += con_ui * basis[j].u * sum_p;
        }
    }
    return result * g_contact;
}

Cd Gaussian_interaction_functional(const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    Cd result(0.0, 0.0);

    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G * compute_H_Mij(c);
            }
            result += con_ui * basis[j].u * sum_p;
        }
    }
    return result * g_gauss;
}

Cd general_cosine_functional(const std::vector<BasisParams>& basis,
                             const std::vector<CosineTerm>& terms) {
    if (terms.empty()) return Cd(0.0, 0.0);

    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    Cd result(0.0, 0.0);
    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                Cd kernel(0.0, 0.0);
                for (const auto& term : terms) {
                    kernel += term.C * compute_cosine_Mij(c, term.A, term.b0);
                }
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G * kernel;
            }
            result += con_ui * basis[j].u * sum_p;
        }
    }
    return result;
}

} // namespace ecg1d
