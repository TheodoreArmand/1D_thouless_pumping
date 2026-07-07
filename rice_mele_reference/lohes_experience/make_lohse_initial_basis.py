#!/usr/bin/env python3
"""
make_lohse_initial_basis.py
===========================

Generate the ECG initial basis CSV for the Lohse pump config
(pumpconfig/lohse_10_5.cpp), in the SAME path-pad style as
initial_state/Vs3Er_Vl3Er/initial_pathpad_N1_K16.csv.

Target state: the band-0 Wannier function of the main.cpp-convention
potential at the schedule start phase phi^m(s=0) = 3*pi/2:

    V(x) = -Vs^m cos(4 pi x/a) - Vl^m cos(2 pi x/a + 3 pi/2)
         = -Vs^m cos(pi x/2)   - Vl^m sin(pi x/4)          (a = 8)

Code units: hbar = m = 1, a = 8, E_rs = pi^2/32.
Vs^m = 5 E_rs, Vl^m = 2.5 E_rs  (paper depths Vs=10 E_rs, Vl=5 E_rs, halved
because the paper's sin^2 depth equals twice the cosine amplitude).

At phi^m = 3*pi/2 (paper phi = 0) the double wells are SYMMETRIC
(Delta = 0, J1 >> J2 across the lowered barrier at x = 2 mod 8), so the
Wannier state is a two-peak bonding orbital over the sites x = 0 and x = 4,
centered at x = 2.

PUMP DIRECTION: the schedule has phi^m DECREASING from 3*pi/2 (that is the
paper's forward pump phi: 0 -> 2pi mapped through phi^m = -phi - pi/2), and
the frame map between conventions is a pure translation, so the packet moves
toward +x: one +a per cycle (0/4 -> 4 -> 4/8 -> 8). Microscopic check:
Delta(phi^m) = eps(x=0) - eps(x=4) = -2*Vl~*cos(phi^m), so just below
3*pi/2 the x=4 site drops and the atom tunnels 0 -> 4, then bonds across
the x=6 barrier toward 8. NOTE this is OPPOSITE to the Vs3Vl3 legacy run,
which drives phi^m increasing and pumps -a.

Basis layout (K = 16, N = 1, A = 0, B real > 0, R real):
  - 11 "fit" Gaussians: least-squares fit of the exact Wannier function
    (centers on/around sites 0 and 4 and the bond center 2).
  - 5 "path-pad" Gaussians with u = 1e-3 seeded along the pump path
    (x = +6, +8, +10, +12) so TDVP has tangent directions to move into.

Checks printed: fit fidelity vs exact Wannier, <H> of the written state vs
exact, and the overlap-matrix condition number compared with the legacy
Vs3Vl3 initial basis.

Output: ../../initial_state/Vs10Ers_Vl5Ers/initial_lohse_N1_K16.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUTDIR = REPO / "initial_state" / "Vs10Ers_Vl5Ers"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTCSV = OUTDIR / "initial_lohse_N1_K16.csv"
LEGACY = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N1_K16.csv"

# ---- code units ----
A_LAT = 8.0
ER = 2 * np.pi**2 / A_LAT**2          # = pi^2/32: short-lattice recoil E_rs
VS, VL = 5.0 * ER, 2.5 * ER           # main.cpp cosine amplitudes
PHI0 = 3 * np.pi / 2                  # schedule start phase phi^m(s=0)


def Vpot(x):
    return -VS * np.cos(4 * np.pi * x / A_LAT) - VL * np.cos(2 * np.pi * x / A_LAT + PHI0)


# ---- exact band-0 Wannier on a ring of NC cells ----
NC, PPC = 16, 128
L = NC * A_LAT
M = NC * PPC
DX = L / M
X = np.arange(M) * DX
PG = 2 * np.pi * np.fft.fftfreq(M, d=DX)

print(f"[1] exact Wannier: ring L={L:g}, M={M} points, dx={DX:g}")
F = np.fft.fft(np.eye(M), axis=0)
TK = np.fft.ifft((0.5 * PG**2)[:, None] * F, axis=0).real
H = TK + np.diag(Vpot(X))
E, PSI = np.linalg.eigh(H)
B0 = PSI[:, :NC]                       # band-0 subspace
Xp = B0.T @ (X[:, None] * B0)
Xp = 0.5 * (Xp + Xp.T)
cc, CC = np.linalg.eigh(Xp)
W = B0 @ CC
j = int(np.argmin(np.abs(cc - L / 2)))
w0 = W[:, j].copy()
w0 /= np.sqrt(np.sum(w0**2) * DX)
if w0[np.argmax(np.abs(w0))] < 0:
    w0 = -w0
xc = float(np.sum(w0**2 * X) * DX)
print(f"    Wannier center = {xc:.4f} (expect near L/2 + 2 = {L/2 + 2:g}, "
      f"bond center of a symmetric double well)")

# shift coordinates so the double-well bond center sits at x = 2
# (packet peaks at x = 0 and x = 4, like the fit dictionary assumes).
# The shift must be an exact multiple of the lattice period a so that the
# ring potential expressed in Xrel equals main.cpp's absolute potential.
shift = A_LAT * round((xc - 2.0) / A_LAT)
Xrel = X - shift
print(f"    ring->line shift = {shift:g} (exact multiple of a; "
      f"Wannier center at Xrel = {xc - shift:.6f})")


def gauss(b, r):
    """exp(-b (x-r)^2) sampled on the (shifted) grid, wrapped on the ring."""
    d = (Xrel - r + L / 2) % L - L / 2
    return np.exp(-b * d**2)


# ---- basis dictionary ----
BB, BN = 0.9, 2.6      # broad / narrow widths [1/x^2]; harmonic on-site B ~ 0.98
fit_spec = [
    (BB, 0.0), (BN, 0.0),          # site x = 0
    (BB, 4.0), (BN, 4.0),          # site x = 4
    (0.35, 2.0), (BB, 2.0),        # bond center (two-peak envelope + barrier)
    (BN, -1.0), (BN, 5.0),         # inner shoulders
    (BB, -4.0), (BB, 8.0),         # neighbour sites (tails)
    (0.35, 6.0),                   # broad tail toward the pump path (+x)
]
pad_spec = [
    (BB, 6.0), (BN, 8.0), (BB, 10.0), (BN, 12.0), (BB, 12.0),
]
K = len(fit_spec) + len(pad_spec)
assert K == 16, K

print(f"[2] least-squares fit with {len(fit_spec)} Gaussians "
      f"+ {len(pad_spec)} path pads (u=1e-3)")
G_fit = np.stack([gauss(b, r) for b, r in fit_spec], axis=1)
u_fit, *_ = np.linalg.lstsq(G_fit, w0, rcond=None)
psi_fit = G_fit @ u_fit
fid_fit = (np.sum(psi_fit * w0) * DX) ** 2 / (np.sum(psi_fit**2) * DX)
print(f"    fit-only fidelity |<fit|W>|^2 = {fid_fit:.6f}")

# assemble full K=16 state: fit coefficients + 1e-3 pads, then normalize
spec = fit_spec + pad_spec
u = np.concatenate([u_fit, 1e-3 * np.ones(len(pad_spec))])
G = np.stack([gauss(b, r) for b, r in spec], axis=1)
psi = G @ u
nrm = np.sqrt(np.sum(psi**2) * DX)
u /= nrm
psi = G @ u

fid = (np.sum(psi * w0) * DX) ** 2 / (np.sum(psi**2) * DX)
E_w0 = float(w0 @ (H @ w0) * DX)
E_fit = float(psi @ (H @ psi) * DX) / (np.sum(psi**2) * DX)
print(f"    full-basis fidelity |<psi|W>|^2 = {fid:.6f}")
print(f"    <H>: exact Wannier = {E_w0:.6f}, fitted state = {E_fit:.6f} "
      f"({abs(E_fit-E_w0)/ER:.4f} E_rs apart)")

# overlap-matrix condition number (TDVP health at t=0), vs legacy basis
S = (G.T @ G) * DX
cond_new = float(np.linalg.cond(S))


def load_legacy_cond():
    rows = [l.strip() for l in open(LEGACY) if l.strip()]
    Bs, Rs = [], []
    for k in range(16):
        _u, _a, b, r, _n = rows[5 * k:5 * k + 5]
        Bs.append(float(b.split(",")[0]))
        Rs.append(float(r.split(",")[0]))
    Gl = np.stack([np.exp(-b * ((X - L / 2 - r + L / 2) % L - L / 2) ** 2)
                   for b, r in zip(Bs, Rs)], axis=1)
    return float(np.linalg.cond((Gl.T @ Gl) * DX))


cond_legacy = load_legacy_cond()
print(f"    overlap cond: new = {cond_new:.3e}, legacy Vs3Vl3 = {cond_legacy:.3e}")

# ---- write CSV in main.cpp load_basis_csv format ----
# per basis function: u / A (N^2) / B (N^2) / R (N) as "re,im" lines + name
with open(OUTCSV, "w") as f:
    for k, ((b, r), uk) in enumerate(zip(spec, u)):
        f.write(f"{uk:.17g},0\n")          # u
        f.write("0,0\n")                   # A = 0 (frozen width channel)
        f.write(f"{b:.17g},0\n")           # B
        f.write(f"{r:.17g},0\n")           # R (absolute main.cpp coordinates)
        f.write(f"{k}\n")
print(f"[write] {OUTCSV}")

# reload-and-verify: reconstruct from the CSV in the Xrel frame
rows = [l.strip() for l in open(OUTCSV) if l.strip()]
assert len(rows) == 5 * K
psi_re = np.zeros(M)
for k in range(K):
    uk = float(rows[5 * k].split(",")[0])
    bk = float(rows[5 * k + 2].split(",")[0])
    rk = float(rows[5 * k + 3].split(",")[0])
    d = (Xrel - rk + L / 2) % L - L / 2
    psi_re += uk * np.exp(-bk * d**2)
fid_reload = (np.sum(psi_re * w0) * DX) ** 2 / (np.sum(psi_re**2) * DX)
print(f"    reload check fidelity = {fid_reload:.6f} (must equal full-basis fidelity)")
assert abs(fid_reload - fid) < 1e-9
print("done.")
