#pragma once
#include "pair_cache.hpp"

namespace ecg1d {

// Kinetic energy kernel: P_Mij
// = -2*tr(K_Mj) + g^T g - 4*g^T K_Mj mu + 4*mu^T K_Mj^2 mu + 2*tr(K_Mj^2 K_inv)
Cd compute_P_Mij(const PairCache& c);

// Second-position-moment kernel: sum_a <x_a^2> / M_G
// = 0.5*tr(K_inv) + mu^T mu
Cd compute_rTr_Mij(const PairCache& c);

// Delta contact interaction kernel for particle pair (a,b)
// G_Mijab = 1/sqrt(pi*h) * exp(-p^2/h)
// where h = K_inv[a,a] + K_inv[b,b] - 2*K_inv[a,b], p = mu[a] - mu[b]
Cd compute_G_Mijab(const PairCache& c, int a, int b);

// Sum of G_Mijab over all pairs a < b
Cd compute_G_Mij(const PairCache& c);

// Gaussian interaction kernel for particle pair (a,b)
// H_Mijab = sigma / sqrt(sigma^2 + h) * exp(-p^2 / (sigma^2 + h))
Cd compute_H_Mijab(const PairCache& c, int a, int b);

// Sum of H_Mijab over all pairs a < b
Cd compute_H_Mij(const PairCache& c);

// General one-body cosine kernel for particle a:
// exp[-A^2 K_inv[a,a]/4] cos(b0 + A mu[a]).
Cd compute_cosine_Mija(const PairCache& c, int a, double A, double b0);

// Sum over all particles for C*cos(A*x + b0), without the strength C.
Cd compute_cosine_Mij(const PairCache& c, double A, double b0);

} // namespace ecg1d
