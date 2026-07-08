#pragma once

namespace ecg1d {

constexpr double hbar   = 1.0;
// Reference-unit convention used by rice_mele_reference:
//   d_s = 1, E_{r,s} = 1, k_s = pi, so hbar^2/(2m) = 1/pi^2.
inline double mass   = 1.0;   // mass in code units
constexpr double g_contact = 1.0;   // delta contact coupling
// Gaussian pair interaction V_int(x_a - x_b) = g_gauss * exp(-(x_a-x_b)^2 / sigma_gauss^2)
// (this is exactly what compute_H_Mijab integrates; no 1/(2 sigma^2), no norm prefactor).
// Mutable like `mass`: a driver may set these ONCE in main() before the first
// OpenMP parallel region; they are read-only during assembly (no race).
inline double g_gauss = 1.0;        // Gaussian interaction coupling
inline double sigma_gauss = 1.0;    // Gaussian interaction width

} // namespace ecg1d
