#include "2gaussian_src/gaussian_terms.hpp"
#include "derivatives.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "nointeraction_src/pump_common.hpp"
#include "physical_constants.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using ecg1d::BasisParams;
using ecg1d::Cd;

double gaussian_energy(const std::vector<BasisParams>& basis) {
    return (ecg1d::Gaussian_interaction_functional(basis) / ecg1d::overlap(basis)).real();
}

void normalize(std::vector<BasisParams>& basis) {
    const double nrm = ecg1d::overlap(basis).real();
    const double scale = std::sqrt(1.0 / nrm);
    for (auto& b : basis) b.u *= scale;
}

void perturb_real(std::vector<BasisParams>& basis, int a1, int a2, int a3, int a4, double eps) {
    if (a1 == 1) {
        basis[a2].u += eps;
    } else if (a1 == 2) {
        basis[a2].B(a3, a3) += eps;
    } else if (a1 == 3) {
        basis[a2].R(a3) += eps;
    } else if (a1 == 4) {
        basis[a2].A(a3, a4) += eps;
        if (a3 != a4) basis[a2].A(a4, a3) += eps;
    }
}

void check(const std::vector<BasisParams>& basis,
           const std::string& label,
           int a1, int a2, int a3, int a4,
           double eps) {
    auto plus = basis;
    auto minus = basis;
    perturb_real(plus, a1, a2, a3, a4, eps);
    perturb_real(minus, a1, a2, a3, a4, -eps);
    const double fd = (gaussian_energy(plus) - gaussian_energy(minus)) / (2.0 * eps);

    const Cd g_false =
        ecg1d::calculate_Hamiltonian_gaussian_partial(a1, a2, a3, a4, false, basis);
    const Cd g_true =
        ecg1d::calculate_Hamiltonian_gaussian_partial(a1, a2, a3, a4, true, basis);
    const double wirtinger_real = (g_false + g_true).real();

    std::cout << std::left << std::setw(10) << label
              << " fd=" << std::scientific << std::setprecision(12) << fd
              << " analytic_real=" << wirtinger_real
              << " abs_diff=" << std::abs(fd - wirtinger_real)
              << " g_false=" << g_false
              << " g_true=" << g_true << "\n";
}

}  // namespace

int main() {
    auto basis = pump2::load_basis_csv("initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv", 2);
    normalize(basis);
    ecg1d::g_gauss = 0.3 * (2.0 * M_PI * M_PI / (8.0 * 8.0));
    ecg1d::sigma_gauss = 1.0;

    std::cout << "E_gauss=" << std::setprecision(17) << gaussian_energy(basis) << "\n";
    std::cout << "Finite-difference check of real parameter perturbations.\n";
    check(basis, "u0", 1, 0, 0, 0, 1.0e-6);
    check(basis, "B00_0", 2, 0, 0, 0, 1.0e-6);
    check(basis, "B11_0", 2, 0, 1, 0, 1.0e-6);
    check(basis, "R0_0", 3, 0, 0, 0, 1.0e-6);
    check(basis, "R1_0", 3, 0, 1, 0, 1.0e-6);
    check(basis, "A00_0", 4, 0, 0, 0, 1.0e-6);
    check(basis, "A01_0", 4, 0, 0, 1, 1.0e-6);
    check(basis, "A11_0", 4, 0, 1, 1, 1.0e-6);
    return 0;
}
