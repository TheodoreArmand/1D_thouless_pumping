#!/usr/bin/env python3
"""
Generate a path-pad N=2 initial ECG basis for the Vs3/Vl3 3/3 pump.

The input one-particle state is

  initial_state/Vs3Er_Vl3Er/initial_pathpad_N1_K16.csv

which already contains the N=1 fit plus path-pad Gaussians seeded along the
legacy 3/3 pump direction (-x). This script builds an adjacent-cell bosonic
product state:

  psi(x0, x1) ~ phi_left(x0) phi_right(x1) + phi_right(x0) phi_left(x1),

where phi_right is the same N=1 packet shifted by one long-lattice period
a = 8. The output keeps K=32:

  - 16 core product Gaussians selected from the dominant N=1 fit terms;
  - 16 path-pad reserve products, eight for the left packet and eight for the
    right packet. Each packet gets the same four downstream centers as the N=1
    path-pad state, with both broad and narrow widths.

The coefficients are least-squares fitted to the full K16 x K16 symmetrized
product on a periodic grid, then normalized in that same grid convention.

Output:
  initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUTDIR = REPO / "initial_state" / "Vs3Er_Vl3Er"
N1CSV = OUTDIR / "initial_pathpad_N1_K16.csv"
OUTCSV = OUTDIR / "initial_pathpad_N2_K32.csv"

A_LAT = 8.0
LENGTH = 16 * A_LAT
M_GRID = 2048
K_OUT = 32


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


def write_n2_csv(path: Path, selected, coeffs, n1):
    with path.open("w") as f:
        for k, ((i, j), coeff) in enumerate(zip(selected, coeffs)):
            _ki, _ui, bi, ri = n1[i]
            _kj, _uj, bj, rj = n1[j]
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
            f.write(f"{rj + A_LAT:.17g},0\n")
            f.write(f"{k}\n")


def main() -> None:
    if not N1CSV.exists():
        raise RuntimeError(f"missing {N1CSV}")

    n1 = load_n1_basis(N1CSV)
    if len(n1) != 16:
        raise RuntimeError(f"expected K=16 one-particle basis, got {len(n1)}")

    fit_terms = list(range(8))
    pad_terms = list(range(8, 16))

    core_products = []
    for i in fit_terms:
        for j in fit_terms:
            core_products.append((abs(n1[i][1] * n1[j][1]), i, j))
    core_products.sort(reverse=True, key=lambda row: row[0])
    selected = [(i, j) for _score, i, j in core_products[:16]]

    # N=1 path pads at R = -2, -4, -6, -8, with broad+narrow widths.
    # Include the same downstream pads once for each packet. The right packet
    # is shifted by +a, so these become x = 6, 4, 2, 0 rather than +x pads.
    path_pads = list(range(8, 16))
    anchor = 0
    selected += [(p, anchor) for p in path_pads]
    selected += [(anchor, p) for p in path_pads]
    if len(selected) != K_OUT:
        raise RuntimeError(f"internal selection error: {len(selected)} terms")

    dx = LENGTH / M_GRID
    x = np.arange(M_GRID) * dx
    xrel = x - LENGTH / 2

    g_cache = {}

    def get_g(i: int, shift: float):
        key = (i, shift)
        if key not in g_cache:
            _k, _u, b, r = n1[i]
            g_cache[key] = gauss_on_ring(xrel, LENGTH, b, r + shift)
        return g_cache[key]

    psi_left = np.zeros(M_GRID)
    psi_right = np.zeros(M_GRID)
    for i, u, _b, _r in n1:
        psi_left += u * get_g(i, 0.0)
        psi_right += u * get_g(i, A_LAT)

    psi_left /= np.sqrt(np.sum(psi_left * psi_left) * dx)
    psi_right /= np.sqrt(np.sum(psi_right * psi_right) * dx)
    psi_target = (psi_left[:, None] * psi_right[None, :]
                  + psi_right[:, None] * psi_left[None, :])
    psi_target /= np.sqrt(np.sum(psi_target * psi_target) * dx * dx)

    columns = []
    for i, j in selected:
        columns.append(sym_product(get_g(i, 0.0), get_g(j, A_LAT)).reshape(-1))
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
    for k, (i, j) in enumerate(selected):
        tag = "core" if k < 16 else "path-pad"
        print(f"    {k:2d} {tag:8s}: left={i:2d} right={j:2d} coeff={coeffs[k]:+.8e}")
    print(f"[2] fidelity vs full N1 path-pad sym product = {fidelity:.10f}")
    print(f"    grid norm check = {np.sum(psi_fit * psi_fit) * dx * dx:.12f}")

    write_n2_csv(OUTCSV, selected, coeffs, n1)
    print(f"[write] {OUTCSV}")

    # Basic CSV shape check: N=2 rows are u + A(4) + B(4) + R(2) + name = 12.
    rows = [line.strip() for line in OUTCSV.read_text().splitlines() if line.strip()]
    if len(rows) != 12 * K_OUT:
        raise RuntimeError(f"bad output row count: {len(rows)}")
    print("done.")


if __name__ == "__main__":
    main()
