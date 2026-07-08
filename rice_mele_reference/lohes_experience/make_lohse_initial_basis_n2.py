#!/usr/bin/env python3
"""
Generate the N=2 adjacent-cell Lohse initial ECG basis.

The state is the bosonic product of two band-0 Wannier packets centered at
x ~= 2 and x ~= 10. The ECG basis uses the engine's permutation symmetrisation:
one product Gaussian (left_i, right_j) represents the normalized symmetric
combination of that product and its particle-swapped partner.

Output:
  initial_state/Vs10Ers_Vl5Ers/initial_lohse_N2_K24.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUTDIR = REPO / "initial_state" / "Vs10Ers_Vl5Ers"
OUTDIR.mkdir(parents=True, exist_ok=True)
N1CSV = OUTDIR / "initial_lohse_N1_K16.csv"
OUTCSV = OUTDIR / "initial_lohse_N2_K24.csv"

A_LAT = 8.0
ER = 2 * np.pi**2 / A_LAT**2
VS, VL = 5.0 * ER, 2.5 * ER
PHI0 = 3 * np.pi / 2
K_OUT = 24
N_KEEP_1D = 6


def vpot(x: np.ndarray) -> np.ndarray:
    return -VS * np.cos(4 * np.pi * x / A_LAT) - VL * np.cos(2 * np.pi * x / A_LAT + PHI0)


def load_n1_basis(path: Path):
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(rows) % 5 != 0:
        raise RuntimeError(f"unexpected N=1 CSV length: {path}")
    basis = []
    for k in range(len(rows) // 5):
        u = complex(*map(float, rows[5 * k].split(",")))
        b = complex(*map(float, rows[5 * k + 2].split(",")))
        r = complex(*map(float, rows[5 * k + 3].split(",")))
        basis.append((k, u.real, b.real, r.real))
    return basis


def exact_wannier():
    nc, ppc = 16, 64
    length = nc * A_LAT
    m = nc * ppc
    dx = length / m
    x = np.arange(m) * dx
    p = 2 * np.pi * np.fft.fftfreq(m, d=dx)

    print(f"[1] exact N=1 Wannier grid: L={length:g}, M={m}, dx={dx:g}")
    fmat = np.fft.fft(np.eye(m), axis=0)
    tmat = np.fft.ifft((0.5 * p**2)[:, None] * fmat, axis=0).real
    hmat = tmat + np.diag(vpot(x))
    eigvals, eigvecs = np.linalg.eigh(hmat)
    band0 = eigvecs[:, :nc]
    xpos = band0.T @ (x[:, None] * band0)
    xpos = 0.5 * (xpos + xpos.T)
    centers, coeffs = np.linalg.eigh(xpos)
    wanniers = band0 @ coeffs
    j = int(np.argmin(np.abs(centers - length / 2)))
    w = wanniers[:, j].copy()
    w /= np.sqrt(np.sum(w**2) * dx)
    if w[np.argmax(np.abs(w))] < 0:
        w = -w
    xc = float(np.sum(w**2 * x) * dx)
    shift = A_LAT * round((xc - 2.0) / A_LAT)
    xrel = x - shift
    print(f"    center after shift: {xc - shift:.6f}")
    return xrel, w, dx, length


def gauss(xrel: np.ndarray, length: float, b: float, r: float) -> np.ndarray:
    d = (xrel - r + length / 2) % length - length / 2
    return np.exp(-b * d**2)


def sym_product(g_left_a, g_right_b, g_right_a, g_left_b):
    return (g_left_a[:, None] * g_right_b[None, :]
            + g_right_a[:, None] * g_left_b[None, :]) / np.sqrt(2.0)


def kinetic_plus_lattice_energy(psi: np.ndarray, xrel: np.ndarray, dx: float) -> float:
    m = psi.shape[0]
    p = 2 * np.pi * np.fft.fftfreq(m, d=dx)
    p2 = p[:, None] ** 2 + p[None, :] ** 2
    hpsi = np.fft.ifft2(0.5 * p2 * np.fft.fft2(psi)).real
    v = vpot(xrel)
    hpsi += (v[:, None] + v[None, :]) * psi
    norm = np.sum(psi * psi) * dx * dx
    return float(np.sum(psi * hpsi) * dx * dx / norm)


def gaussian_interaction_energy(psi: np.ndarray, xrel: np.ndarray, dx: float,
                                g_code: float, sigma: float) -> float:
    d = xrel[:, None] - xrel[None, :]
    vint = g_code * np.exp(-(d * d) / (sigma * sigma))
    norm = np.sum(psi * psi) * dx * dx
    return float(np.sum(psi * psi * vint) * dx * dx / norm)


def main() -> None:
    if not N1CSV.exists():
        raise RuntimeError(f"missing {N1CSV}; run make_lohse_initial_basis.py first")

    n1 = load_n1_basis(N1CSV)
    chosen = sorted(n1, key=lambda row: abs(row[1]), reverse=True)[:N_KEEP_1D]
    print("[2] selected one-particle Gaussians:")
    for k, u, b, r in chosen:
        print(f"    k={k:2d} |u|={abs(u):.6e} B={b:.6g} R={r:.6g}")

    products = []
    for left in chosen:
        for right in chosen:
            i, ui, bi, ri = left
            j, uj, bj, rj = right
            products.append((abs(ui * uj), i, j, ui * uj, bi, bj, ri, rj + A_LAT))
    products.sort(reverse=True, key=lambda row: row[0])
    products = products[:K_OUT]

    xrel, w_left, dx, length = exact_wannier()
    shift_pts = int(round(A_LAT / dx))
    w_right = np.roll(w_left, shift_pts)
    psi_exact = (w_left[:, None] * w_right[None, :]
                 + w_right[:, None] * w_left[None, :])
    psi_exact /= np.sqrt(np.sum(psi_exact * psi_exact) * dx * dx)

    g_cache = {}

    def get_g(b: float, r: float):
        key = (float(b), float(r))
        if key not in g_cache:
            g_cache[key] = gauss(xrel, length, b, r)
        return g_cache[key]

    psi = np.zeros((xrel.size, xrel.size), dtype=float)
    for _, _i, _j, coeff, bi, bj, ri, rj in products:
        gl = get_g(bi, ri)
        gr = get_g(bj, rj)
        gl_swap = get_g(bi, ri)
        gr_swap = get_g(bj, rj)
        psi += coeff * sym_product(gl, gr, gr_swap, gl_swap)

    norm = np.sqrt(np.sum(psi * psi) * dx * dx)
    products = [(score, i, j, coeff / norm, bi, bj, ri, rj)
                for score, i, j, coeff, bi, bj, ri, rj in products]
    psi /= norm

    overlap = float(np.sum(psi_exact * psi) * dx * dx)
    fidelity = overlap * overlap
    e_free = kinetic_plus_lattice_energy(psi, xrel, dx)
    v_int = gaussian_interaction_energy(psi, xrel, dx, 0.3 * ER, 1.0)
    print(f"[3] N=2 grid checks on {xrel.size}^2 grid")
    print(f"    fidelity vs exact sym product = {fidelity:.8f}")
    print(f"    <H_free> = {e_free:.10f}")
    print(f"    <V_gauss> for g=0.3 E_rs, sigma=1 = {v_int:.10e}")
    print(f"    <H_gauss> = {e_free + v_int:.10f}")

    with OUTCSV.open("w") as f:
        for k, (_score, i, j, coeff, bi, bj, ri, rj) in enumerate(products):
            f.write(f"{coeff:.17g},0\n")
            f.write("0,0\n")   # A00
            f.write("0,0\n")   # A01
            f.write("0,0\n")   # A10
            f.write("0,0\n")   # A11
            f.write(f"{bi:.17g},0\n")
            f.write("0,0\n")
            f.write("0,0\n")
            f.write(f"{bj:.17g},0\n")
            f.write(f"{ri:.17g},0\n")
            f.write(f"{rj:.17g},0\n")
            # BasisParams::name must equal the C++ basis index; the derivative
            # code uses it to match AlphaIndex::a2 to bra/ket basis functions.
            f.write(f"{k}\n")
            print(f"    basis {k:2d}: left={i:2d} right={j:2d} u={coeff:+.6e}")
    print(f"[write] {OUTCSV}")


if __name__ == "__main__":
    main()
