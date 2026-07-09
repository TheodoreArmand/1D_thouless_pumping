#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include <cstddef>
#include <vector>

namespace ecg1d {

struct PairCache {
    int N;
    MatrixXcd K;       // (N,N) full K matrix
    MatrixXcd K_inv;   // K^{-1}
    VectorXcd b;       // (N,) linear term
    Cd C_val;          // scalar constant
    Cd det_K;          // det(K)
    Cd M_G;            // overlap matrix element
    VectorXcd mu;      // 0.5 * K_inv * b (mean position)
    MatrixXcd K_Mj;    // P^T A_j P + P^T B_j P (ket-side K)
    VectorXcd g_Mj;    // 2 * R_j^T B_j P (ket-side linear term, as column vector)

    static PairCache build(const BasisParams& pi, const BasisParams& pj,
                           const MatrixXi& perm_matrix);
};

class PairCacheTable {
public:
    explicit PairCacheTable(const std::vector<BasisParams>& basis);

    const PairCache& get(int i, int j, int p) const;
    int basis_size() const { return K_; }
    int particle_count() const { return perms_.N; }
    int permutation_count() const { return perms_.SN; }
    int permutation_index(const MatrixXi& perm_matrix) const;

private:
    std::size_t offset(int i, int j, int p) const;

    const std::vector<BasisParams>* basis_;
    PermutationSet perms_;
    int K_;
    std::vector<PairCache> entries_;
};

class PairCacheTableScope {
public:
    explicit PairCacheTableScope(const PairCacheTable& table);
    ~PairCacheTableScope();

    PairCacheTableScope(const PairCacheTableScope&) = delete;
    PairCacheTableScope& operator=(const PairCacheTableScope&) = delete;

private:
    const PairCacheTable* previous_;
};

} // namespace ecg1d
