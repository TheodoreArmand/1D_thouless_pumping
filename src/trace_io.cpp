#include "trace_io.hpp"

#include "hamiltonian.hpp"
#include "interaction_kernels.hpp"
#include "pair_cache.hpp"
#include "permutation.hpp"
#include "physical_constants.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace ecg1d {

WidthMonitors compute_width_monitors(const std::vector<BasisParams>& basis) {
    // Min over basis functions k AND particles a of the diagonal widths.
    // argmin encodes both as k*N + a (equals k for N = 1, so the N = 1
    // trace output is unchanged).
    WidthMonitors w;
    double best_b = std::numeric_limits<double>::infinity();
    double best_ab = std::numeric_limits<double>::infinity();
    for (int k = 0; k < static_cast<int>(basis.size()); ++k) {
        const int N = basis[k].N();
        for (int a = 0; a < N; ++a) {
            const double reB = basis[k].B(a, a).real();
            const double reAB = (basis[k].A(a, a) + basis[k].B(a, a)).real();
            if (reB < best_b) {
                best_b = reB;
                w.min_re_B = reB;
                w.argmin_re_B = k * N + a;
            }
            if (reAB < best_ab) {
                best_ab = reAB;
                w.min_re_AplusB = reAB;
                w.argmin_re_AplusB = k * N + a;
            }
        }
    }
    return w;
}

void sample_observables(const std::vector<BasisParams>& basis,
                        double t,
                        double phi,
                        const HamiltonianTerms& terms,
                        double lattice_a,
                        Trace& trace) {
    const Cd S = overlap(basis);
    const Cd T_raw = kinetic_energy_functional(basis);
    const Cd V_raw = general_cosine_functional(basis, terms.cosine_terms);
    const Cd E = compute_total_energy(basis, terms);
    const double V_const = terms.one_body_constant * basis[0].N();
    const int N = basis[0].N();
    const bool include_v_gauss = (N >= 2);
    const double V_gauss = (include_v_gauss && terms.gaussian)
        ? (Gaussian_interaction_functional(basis) / S).real()
        : 0.0;

    double x_mean_v = 0.0;
    double p_mean_v = 0.0;
    const int K = static_cast<int>(basis.size());
    const PermutationSet perms = PermutationSet::generate(N);
    {
        Cd amp_x(0.0, 0.0);
        Cd amp_p(0.0, 0.0);
        for (int i = 0; i < K; ++i) {
            const Cd conj_ui = std::conj(basis[i].u);
            for (int j = 0; j < K; ++j) {
                Cd term_x(0.0, 0.0);
                Cd term_p(0.0, 0.0);
                for (int p = 0; p < perms.SN; ++p) {
                    PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                    const double sign = static_cast<double>(perms.signs[p]);
                    const Cd mu_sum = c.mu.sum();
                    const Cd Km_mu = (c.K_Mj * c.mu).sum();
                    const Cd g_sum = c.g_Mj.sum();
                    term_x += sign * c.M_G * mu_sum;
                    term_p += sign * Cd(0.0, -1.0) * c.M_G * (-2.0 * Km_mu + g_sum);
                }
                amp_x += conj_ui * basis[j].u * term_x;
                amp_p += conj_ui * basis[j].u * term_p;
            }
        }
        x_mean_v = (amp_x / S).real();
        p_mean_v = (amp_p / S).real();
    }

    const double p2 = 2.0 * mass * (T_raw / S).real();
    double x2 = 0.0;
    {
        Cd amp_x2(0.0, 0.0);
        for (int i = 0; i < K; ++i) {
            const Cd conj_ui = std::conj(basis[i].u);
            for (int j = 0; j < K; ++j) {
                Cd term(0.0, 0.0);
                for (int p = 0; p < perms.SN; ++p) {
                    PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                    // sum_a <x_a^2> kernel; for N = 1 this equals the old
                    // expression 0.5*K_inv(0,0) + mu(0)*mu(0) bit-for-bit.
                    const Cd second = compute_rTr_Mij(c);
                    term += static_cast<double>(perms.signs[p]) * c.M_G * second;
                }
                amp_x2 += conj_ui * basis[j].u * term;
            }
        }
        x2 = (amp_x2 / S).real();
    }

    trace.t.push_back(t);
    trace.tau.push_back(phi);
    trace.norm.push_back(S.real());
    trace.E_total.push_back(E.real());
    trace.T_kin.push_back((T_raw / S).real());
    trace.V_cos.push_back((V_raw / S).real());
    trace.V_const.push_back(V_const);
    trace.V_lattice.push_back((V_raw / S).real() + V_const);
    if (include_v_gauss) trace.V_gauss.push_back(V_gauss);
    trace.x_mean.push_back(x_mean_v);
    trace.p_mean.push_back(p_mean_v);
    trace.x2.push_back(x2);
    trace.p2.push_back(p2);
    trace.polarization_cell.push_back(x_mean_v / lattice_a);
    trace.basis_size.push_back(K);
}

void append_trace_diagnostics(Trace& trace,
                              const WidthMonitors& width,
                              const RealtimeStepResult* step_result) {
    const double dnan = std::numeric_limits<double>::quiet_NaN();
    trace.raw_cond.push_back(step_result ? step_result->raw_cond : dnan);
    trace.actual_solve_cond.push_back(step_result ? step_result->actual_solve_cond : dnan);
    trace.actual_solve_rank.push_back(step_result ? step_result->effective_rank : -1);
    trace.sv_min.push_back(step_result ? step_result->sv_small[0] : dnan);
    trace.relative_raw_residual.push_back(step_result ? step_result->relative_raw_residual : dnan);
    trace.discarded_rhs_fraction.push_back(step_result ? step_result->discarded_rhs_fraction : dnan);
    trace.dz_norm.push_back(step_result ? step_result->dz_norm : dnan);
    trace.metric_norm.push_back(step_result ? step_result->metric_norm : dnan);
    trace.min_re_B.push_back(width.min_re_B);
    trace.argmin_re_B.push_back(width.argmin_re_B);
    trace.min_re_AplusB.push_back(width.min_re_AplusB);
    trace.argmin_re_AplusB.push_back(width.argmin_re_AplusB);
}

void write_trace_csv(const std::string& path, const Trace& tr) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    const bool include_v_gauss = !tr.V_gauss.empty();
    if (include_v_gauss && tr.V_gauss.size() != tr.t.size()) {
        throw std::runtime_error("trace V_gauss length does not match sample count");
    }
    f << "t,phi,norm,E_total,T_kin,V_cos,V_const,V_lattice,x_mean,p_mean,x2,p2,polarization_cell,basis_size,"
         "raw_cond,actual_solve_cond,actual_solve_rank,sv_min,"
         "relative_raw_residual,discarded_rhs_fraction,"
         "dz_norm,metric_norm,min_re_B,argmin_re_B,min_re_AplusB,argmin_re_AplusB";
    if (include_v_gauss) f << ",V_gauss";
    f << "\n";
    f << std::setprecision(17);
    for (size_t i = 0; i < tr.t.size(); ++i) {
        f << tr.t[i] << ","
          << tr.tau[i] << ","
          << tr.norm[i] << ","
          << tr.E_total[i] << ","
          << tr.T_kin[i] << ","
          << tr.V_cos[i] << ","
          << tr.V_const[i] << ","
          << tr.V_lattice[i] << ","
          << tr.x_mean[i] << ","
          << tr.p_mean[i] << ","
          << tr.x2[i] << ","
          << tr.p2[i] << ","
          << tr.polarization_cell[i] << ","
          << tr.basis_size[i] << ","
          << tr.raw_cond[i] << ","
          << tr.actual_solve_cond[i] << ","
          << tr.actual_solve_rank[i] << ","
          << tr.sv_min[i] << ","
          << tr.relative_raw_residual[i] << ","
          << tr.discarded_rhs_fraction[i] << ","
          << tr.dz_norm[i] << ","
          << tr.metric_norm[i] << ","
          << tr.min_re_B[i] << ","
          << tr.argmin_re_B[i] << ","
          << tr.min_re_AplusB[i] << ","
          << tr.argmin_re_AplusB[i];
        if (include_v_gauss) f << "," << tr.V_gauss[i];
        f << "\n";
    }
}

}  // namespace ecg1d
