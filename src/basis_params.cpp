#include "basis_params.hpp"

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

namespace ecg1d {

BasisParams BasisParams::conj_params() const {
    BasisParams result;
    result.u = std::conj(u);
    result.A = A.conjugate();
    result.B = B.conjugate();
    result.R = R.conjugate();
    result.name = name;
    return result;
}

BasisParams BasisParams::from_arrays(Cd u, const MatrixXcd& A, const MatrixXcd& B,
                                     const VectorXcd& R, int name) {
    BasisParams bp;
    bp.u = u;
    bp.A = A;
    bp.B = B;
    bp.R = R;
    bp.name = name;
    return bp;
}

double min_real_width_eigenvalue(const BasisParams& basis) {
    const int N = basis.N();
    if (N <= 0 || basis.A.rows() != N || basis.A.cols() != N ||
        basis.B.rows() != N || basis.B.cols() != N) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    Eigen::MatrixXd re_width(N, N);
    const MatrixXcd width = basis.A + basis.B;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const double rij = width(i, j).real();
            const double rji = width(j, i).real();
            if (!std::isfinite(rij) || !std::isfinite(rji)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            re_width(i, j) = 0.5 * (rij + rji);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(re_width, Eigen::EigenvaluesOnly);
    if (eig.info() != Eigen::Success || eig.eigenvalues().size() == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return eig.eigenvalues()(0);
}

} // namespace ecg1d
