#pragma once

#include "types.hpp"

namespace ecg1d {

struct BasisParams {
    Cd u;
    MatrixXcd A;  // (N,N)
    MatrixXcd B;  // (N,N) diagonal
    VectorXcd R;  // (N,)
    int name;

    int N() const { return static_cast<int>(R.size()); }

    BasisParams conj_params() const;

    static BasisParams from_arrays(Cd u, const MatrixXcd& A, const MatrixXcd& B,
                                   const VectorXcd& R, int name = 0);
};

// Normalizability monitor for one primitive. For N > 1 the width condition is
// positive definiteness of Re(A + B), not just positivity of its diagonal.
double min_real_width_eigenvalue(const BasisParams& basis);

} // namespace ecg1d
