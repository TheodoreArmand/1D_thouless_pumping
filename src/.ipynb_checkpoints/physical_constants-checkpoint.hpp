#pragma once

namespace ecg1d {

constexpr double hbar   = 1.0;
// Reference-unit convention used by rice_mele_reference:
//   d_s = 1, E_{r,s} = 1, k_s = pi, so hbar^2/(2m) = 1/pi^2.
inline double mass   = 4.934802200544679;  // pi^2 / 2 with hbar = 1
inline double omega  = 1.0;   // harmonic-trap frequency; set at runtime from ECG_TP_OMEGA
constexpr double g_contact = 1.0;   // delta contact coupling
constexpr double k_L    = 0.5;      // kicking lattice wavevector
constexpr double kappa  = 0.83;     // kicking amplitude
constexpr double g_gauss = 1.0;     // Gaussian interaction coupling
constexpr double sigma_gauss = 1.0; // Gaussian interaction width

} // namespace ecg1d
