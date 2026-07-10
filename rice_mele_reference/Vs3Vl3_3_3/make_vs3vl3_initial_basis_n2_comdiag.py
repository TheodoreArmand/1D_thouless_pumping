#!/usr/bin/env python3
"""
Generate a COM-diagonal reserve N=2 K32 basis for the Vs3/Vl3 3/3 pump.

This keeps the same 16 dominant core product terms as the path-pad K32 basis,
but replaces the 16 one-particle-anchor path pads with four two-particle
COM-diagonal centers,

  (-2, 6), (-4, 4), (-6, 2), (-8, 0),

each with the four broad/narrow joint width combinations. Coefficients are
least-squares fitted to the same initial symmetrized product target used by the
path-pad generator, then normalized in that grid convention.

Output:
  initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32_COMdiag.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUTDIR = REPO / "initial_state" / "Vs3Er_Vl3Er"
N1CSV = OUTDIR / "initial_pathpad_N1_K16.csv"
OUTCSV = OUTDIR / "initial_pathpad_N2_K32_COMdiag.csv"

A_LAT = 8.0
LENGTH = 16 * A_LAT
M_GRID = 2048
K_OUT = 32
CORE_COUNT = 16


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


def gauss_on_ring(xrel: np.ndarray, length: float, b: float, r: float) -> np.ndarray:
    d = (xrel - r + length / 2) % length - length / 2
    return np.exp(-b * d * d)


def sym_product(gl: np.ndarray, gr: np.ndarray) -> np.ndarray:
    return (gl[:, None] * gr[None, :] + gr[:, None] * gl[None, :]) / np.sqrt(2.0)


def write_n2_csv(path: Path, terms, coeffs):
    with path.open("w") as f:
        for k, (term, coeff) in enumerate(zip(terms, coeffs)):
            _tag, bi, ri, bj, rj = term
            f.write(f"{coeff:.17g},0\n")
            f.write("0,0\n")  # A00
            f.write("0,0\n")  # A01
            f.write("0,0\n")  # A10
            f.write("0,0\n")  # A11
            f.write(f"{bi:.17g},0\n")
            f.write("0,0\n")
            f.write("0,0\n")
            f.write(f"{bj:.17g},0\n")
            f.write(f"{ri:.17g},0\n")
            f.write(f"{rj:.17g},0\n")
            f.write(f"{k}\n")


def main() -> None:
    if not N1CSV.exists():
        raise RuntimeError(f"missing {N1CSV}")

    n1 = load_n1_basis(N1CSV)
    if len(n1) != 16:
        raise RuntimeError(f"expected K=16 one-particle basis, got {len(n1)}")

    fit_terms = list(range(8))
    core_products = []
    for i in fit_terms:
        for j in fit_terms:
            core_products.append((abs(n1[i][1] * n1[j][1]), i, j))
    core_products.sort(reverse=True, key=lambda row: row[0])

    terms = []
    for _score, i, j in core_products[:CORE_COUNT]:
        _ki, _ui, bi, ri = n1[i]
        _kj, _uj, bj, rj = n1[j]
        terms.append((f"core:{i},{j}", bi, ri, bj, rj + A_LAT))

    broad = n1[8][2]
    narrow = n1[9][2]
    width_pairs = [(broad, broad), (broad, narrow), (narrow, broad), (narrow, narrow)]
    for shift in (-2.0, -4.0, -6.0, -8.0):
        for bi, bj in width_pairs:
            terms.append((f"comdiag:{shift:g}", bi, shift, bj, shift + A_LAT))

    if len(terms) != K_OUT:
        raise RuntimeError(f"internal selection error: {len(terms)} terms")

    dx = LENGTH / M_GRID
    x = np.arange(M_GRID) * dx
    xrel = x - LENGTH / 2

    g_cache = {}

    def get_n1_g(i: int, shift: float):
        key = ("n1", i, shift)
        if key not in g_cache:
            _k, _u, b, r = n1[i]
            g_cache[key] = gauss_on_ring(xrel, LENGTH, b, r + shift)
        return g_cache[key]

    def get_g(b: float, r: float):
        key = ("custom", b, r)
        if key not in g_cache:
            g_cache[key] = gauss_on_ring(xrel, LENGTH, b, r)
        return g_cache[key]

    psi_left = np.zeros(M_GRID)
    psi_right = np.zeros(M_GRID)
    for i, u, _b, _r in n1:
        psi_left += u * get_n1_g(i, 0.0)
        psi_right += u * get_n1_g(i, A_LAT)

    psi_left /= np.sqrt(np.sum(psi_left * psi_left) * dx)
    psi_right /= np.sqrt(np.sum(psi_right * psi_right) * dx)
    psi_target = (psi_left[:, None] * psi_right[None, :]
                  + psi_right[:, None] * psi_left[None, :])
    psi_target /= np.sqrt(np.sum(psi_target * psi_target) * dx * dx)

    columns = []
    for _tag, bi, ri, bj, rj in terms:
        columns.append(sym_product(get_g(bi, ri), get_g(bj, rj)).reshape(-1))
    design = np.stack(columns, axis=1)
    target = psi_target.reshape(-1)
    coeffs, *_ = np.linalg.lstsq(design, target, rcond=None)

    psi_fit = (design @ coeffs).reshape(M_GRID, M_GRID)
    fit_norm = np.sqrt(np.sum(psi_fit * psi_fit) * dx * dx)
    coeffs /= fit_norm
    psi_fit /= fit_norm

    overlap = float(np.sum(psi_fit * psi_target) * dx * dx)
    fidelity = overlap * overlap

    print(f"[1] selected K={K_OUT} products:")
    for k, (term, coeff) in enumerate(zip(terms, coeffs)):
        tag, bi, ri, bj, rj = term
        kind = "core" if k < CORE_COUNT else "comdiag"
        print(
            f"    {k:2d} {kind:7s}: {tag:14s} "
            f"B=({bi:.8g},{bj:.8g}) R=({ri:+.8g},{rj:+.8g}) "
            f"coeff={coeff:+.8e}"
        )
    print(f"[2] fidelity vs full N1 path-pad sym product = {fidelity:.10f}")
    print(f"    grid norm check = {np.sum(psi_fit * psi_fit) * dx * dx:.12f}")

    write_n2_csv(OUTCSV, terms, coeffs)
    print(f"[write] {OUTCSV}")

    rows = [line.strip() for line in OUTCSV.read_text().splitlines() if line.strip()]
    if len(rows) != 12 * K_OUT:
        raise RuntimeError(f"bad output row count: {len(rows)}")
    print("done.")


if __name__ == "__main__":
    main()
