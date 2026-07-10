#!/usr/bin/env python3
"""
Create a conservative A-seeded variant of the Vs3/Vl3 N=2 K32 path-pad basis.

The original K32 basis keeps A=0 for every term.  This variant preserves the
first 16 core product terms exactly and adds a small relative-coordinate
correlation only to the 16 path-pad reserve terms:

    A = alpha [[ 1, -1],
               [-1,  1]]

R is adjusted so that the single-Gaussian mean remains at the original target
center despite the nonzero A:

    R_input = B^{-1} (A + B) m

Output:
  initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32_Aseed0p05.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
INDIR = REPO / "initial_state" / "Vs3Er_Vl3Er"
INCSV = INDIR / "initial_pathpad_N2_K32.csv"
OUTCSV = INDIR / "initial_pathpad_N2_K32_Aseed0p05.csv"

ALPHA = 0.05
K_CORE = 16
N = 2
ROWS_PER_TERM = 12


def parse_cd(text: str) -> complex:
    re, im = text.split(",", 1)
    return complex(float(re), float(im))


def fmt_cd(z: complex) -> str:
    return f"{z.real:.17g},{z.imag:.17g}"


def main() -> None:
    rows = [line.strip() for line in INCSV.read_text().splitlines() if line.strip()]
    if len(rows) % ROWS_PER_TERM != 0:
        raise RuntimeError(f"bad row count in {INCSV}: {len(rows)}")
    k_total = len(rows) // ROWS_PER_TERM
    if k_total != 32:
        raise RuntimeError(f"expected K=32, got {k_total}")

    out: list[str] = []
    seeded = 0
    for k in range(k_total):
        off = ROWS_PER_TERM * k
        u = parse_cd(rows[off])
        A = np.array(
            [[parse_cd(rows[off + 1]), parse_cd(rows[off + 2])],
             [parse_cd(rows[off + 3]), parse_cd(rows[off + 4])]],
            dtype=np.complex128,
        )
        B = np.array(
            [[parse_cd(rows[off + 5]), parse_cd(rows[off + 6])],
             [parse_cd(rows[off + 7]), parse_cd(rows[off + 8])]],
            dtype=np.complex128,
        )
        R = np.array([parse_cd(rows[off + 9]), parse_cd(rows[off + 10])], dtype=np.complex128)
        name = rows[off + 11]

        if k >= K_CORE:
            target_center = R.copy()
            A = ALPHA * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=np.complex128)
            bdiag = np.array([B[0, 0], B[1, 1]], dtype=np.complex128)
            R = ((A + B) @ target_center) / bdiag
            seeded += 1

        out.append(fmt_cd(u))
        out.extend(fmt_cd(A[i, j]) for i in range(N) for j in range(N))
        out.extend(fmt_cd(B[i, j]) for i in range(N) for j in range(N))
        out.extend(fmt_cd(R[i]) for i in range(N))
        out.append(name)

    OUTCSV.write_text("\n".join(out) + "\n")
    print(f"wrote {OUTCSV}")
    print(f"alpha={ALPHA:g}")
    print(f"seeded_path_pad_terms={seeded}")


if __name__ == "__main__":
    main()
