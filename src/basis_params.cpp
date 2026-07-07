#include "basis_params.hpp"

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

} // namespace ecg1d
