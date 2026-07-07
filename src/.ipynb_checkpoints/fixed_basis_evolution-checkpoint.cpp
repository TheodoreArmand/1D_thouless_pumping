#include "fixed_basis_evolution.hpp"

#include <Eigen/Dense>
#include <iostream>

namespace ecg1d {

void free_evolve_fixed_basis(std::vector<BasisParams>& basis,
                             const MatrixXcd& H,
                             const MatrixXcd& S,
                             double duration) {
    const int K = static_cast<int>(basis.size());

    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_S(S);
    if (es_S.info() != Eigen::Success) {
        std::cerr << "free_evolve_fixed_basis: S eigendecomposition failed\n";
        return;
    }

    Eigen::VectorXd w = es_S.eigenvalues();
    const double w_max = w(K - 1);
    const double w_cutoff = w_max * 1.0e-10;

    int n_keep = 0;
    for (int i = 0; i < K; i++) {
        if (w(i) > w_cutoff) n_keep++;
    }
    if (n_keep == 0) {
        std::cerr << "free_evolve_fixed_basis: no positive S subspace kept\n";
        return;
    }

    MatrixXcd V_keep = es_S.eigenvectors().rightCols(n_keep);
    Eigen::VectorXd w_keep = w.tail(n_keep);
    MatrixXcd S_inv_half =
        V_keep * w_keep.array().rsqrt().matrix().asDiagonal() * V_keep.adjoint();

    MatrixXcd H_tilde = S_inv_half * H * S_inv_half;
    H_tilde = (0.5 * (H_tilde + H_tilde.adjoint())).eval();

    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_H(H_tilde);
    if (es_H.info() != Eigen::Success) {
        std::cerr << "free_evolve_fixed_basis: H_tilde eigendecomposition failed\n";
        return;
    }

    VectorXcd u(K);
    for (int i = 0; i < K; i++) u(i) = basis[i].u;

    MatrixXcd S_half =
        V_keep * w_keep.array().sqrt().matrix().asDiagonal() * V_keep.adjoint();
    VectorXcd u_tilde = S_half * u;
    VectorXcd coeffs = es_H.eigenvectors().adjoint() * u_tilde;

    for (int n = 0; n < n_keep; n++) {
        coeffs(n) *= std::exp(Cd(0.0, -1.0) * es_H.eigenvalues()(n) * duration);
    }

    VectorXcd u_new = S_inv_half * (es_H.eigenvectors() * coeffs);
    for (int i = 0; i < K; i++) {
        basis[i].u = u_new(i);
    }
}

} // namespace ecg1d
