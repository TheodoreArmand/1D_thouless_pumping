#include "pair_cache.hpp"
#include <cmath>
#include <stdexcept>

namespace ecg1d {
namespace {

const PairCacheTable* g_active_pair_cache_table = nullptr;

PairCache build_pair_cache_uncached(const BasisParams& pi, const BasisParams& pj,
                                    const MatrixXi& perm_matrix) {
    PairCache c;
    c.N = pi.N();
    int N = c.N;

    BasisParams pi_con = pi.conj_params();

    // Permutation matrix as double for arithmetic
    MatrixXcd P = perm_matrix.cast<double>().cast<Cd>();
    MatrixXcd PT = P.transpose();

    // K = conj(A_i) + conj(B_i) + P^T A_j P + P^T B_j P
    MatrixXcd AA = PT * pj.A * P;
    MatrixXcd BB = PT * pj.B * P;
    c.K = pi_con.A + pi_con.B + AA + BB;

    // K_Mj = P^T A_j P + P^T B_j P (ket side only)
    c.K_Mj = AA + BB;

    // b: bT = 2*conj(R_i)^T*conj(B_i) + 2*R_j^T*B_j*P
    // bT is a row vector; b = bT^T is column
    // conj(R_i)^T * conj(B_i) is a (1,N) row vector
    // R_j^T * B_j * P is a (1,N) row vector
    Eigen::RowVectorXcd ii = pi_con.R.transpose() * pi_con.B;  // (1,N)
    Eigen::RowVectorXcd jj = pj.R.transpose() * pj.B * P;      // (1,N)
    Eigen::RowVectorXcd bT = 2.0 * ii + 2.0 * jj;
    c.b = bT.transpose();  // column vector

    // g_Mj: g_Mj_T = 2 * R_j^T B_j P (row vector); g_Mj = column
    c.g_Mj = (2.0 * jj).transpose();

    // C = -conj(R_i)^T conj(B_i) conj(R_i) - R_j^T B_j R_j
    Cd C_ii = (pi_con.R.transpose() * pi_con.B * pi_con.R)(0);
    Cd C_jj = (pj.R.transpose() * pj.B * pj.R)(0);
    c.C_val = -C_ii - C_jj;

    // K_inv and det via LU
    Eigen::PartialPivLU<MatrixXcd> lu(c.K);
    c.K_inv = lu.inverse();
    c.det_K = lu.determinant();

    // mu = 0.5 * K_inv * b
    c.mu = 0.5 * c.K_inv * c.b;

    // M_G = pi^(N/2) * det(K)^(-1/2) * exp(C + 0.25 * bT * K_inv * b)
    double PI = std::pow(M_PI, N / 2.0);
    Cd sqrt_detK_inv = std::pow(c.det_K, -0.5);
    Cd exponent = c.C_val + 0.25 * (bT * c.K_inv * c.b)(0);
    c.M_G = PI * sqrt_detK_inv * std::exp(exponent);

    return c;
}

}  // namespace

PairCache PairCache::build(const BasisParams& pi, const BasisParams& pj,
                           const MatrixXi& perm_matrix) {
    const PairCacheTable* table = g_active_pair_cache_table;
    if (table != nullptr &&
        pi.N() == table->particle_count() &&
        pj.N() == table->particle_count() &&
        pi.name >= 0 && pj.name >= 0 &&
        pi.name < table->basis_size() &&
        pj.name < table->basis_size()) {
        const int p = table->permutation_index(perm_matrix);
        if (p >= 0) return table->get(pi.name, pj.name, p);
    }

    return build_pair_cache_uncached(pi, pj, perm_matrix);
}

PairCacheTable::PairCacheTable(const std::vector<BasisParams>& basis)
    : basis_(&basis),
      perms_(PermutationSet::generate(basis.empty() ? 0 : basis[0].N())),
      K_(static_cast<int>(basis.size())) {
    if (basis.empty()) {
        throw std::runtime_error("PairCacheTable requires a non-empty basis");
    }

    const std::size_t total =
        static_cast<std::size_t>(K_) * static_cast<std::size_t>(K_) *
        static_cast<std::size_t>(perms_.SN);
    entries_.resize(total);

    #pragma omp parallel for schedule(static)
    for (long long idx = 0; idx < static_cast<long long>(total); ++idx) {
        const int p = static_cast<int>(idx % perms_.SN);
        const int pair_idx = static_cast<int>(idx / perms_.SN);
        const int j = pair_idx % K_;
        const int i = pair_idx / K_;
        entries_[static_cast<std::size_t>(idx)] =
            build_pair_cache_uncached((*basis_)[i], (*basis_)[j], perms_.matrices[p]);
    }
}

const PairCache& PairCacheTable::get(int i, int j, int p) const {
    return entries_[offset(i, j, p)];
}

int PairCacheTable::permutation_index(const MatrixXi& perm_matrix) const {
    for (int p = 0; p < perms_.SN; ++p) {
        const MatrixXi& m = perms_.matrices[p];
        if (m.rows() == perm_matrix.rows() &&
            m.cols() == perm_matrix.cols() &&
            (m.array() == perm_matrix.array()).all()) {
            return p;
        }
    }
    return -1;
}

std::size_t PairCacheTable::offset(int i, int j, int p) const {
    return (static_cast<std::size_t>(i) * static_cast<std::size_t>(K_) +
            static_cast<std::size_t>(j)) *
           static_cast<std::size_t>(perms_.SN) +
           static_cast<std::size_t>(p);
}

PairCacheTableScope::PairCacheTableScope(const PairCacheTable& table)
    : previous_(g_active_pair_cache_table) {
    g_active_pair_cache_table = &table;
}

PairCacheTableScope::~PairCacheTableScope() {
    g_active_pair_cache_table = previous_;
}

} // namespace ecg1d
