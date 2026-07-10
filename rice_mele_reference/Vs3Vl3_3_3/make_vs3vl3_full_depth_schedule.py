#!/usr/bin/env python3
"""Generate the canonical full-depth Vs=Vl=3 E_R pump schedule.

The evolved one-particle lattice is

    V(x, phi) = -Vs cos(4 pi x / a) - Vl cos(2 pi x / a + phi),

with ``a=8`` and ``Vs=Vl=3 E_R`` by default.  The schedule obeys
``dphi/dt proportional to gap(phi)^2`` and is written as ``s=t/T,phi`` for
the C++ and grid drivers.  The band gap is computed with a plane-wave Bloch
Hamiltonian in the same full-depth convention as those drivers.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "gap_adaptive_vs3vl3_full_depth_schedule.csv"


def band_gap_curve(
    *,
    lattice_a: float,
    vs: float,
    vl: float,
    phase_points: int,
    k_points: int,
    plane_wave_cutoff: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the minimum first-band gap over quasimomentum for each phase."""
    reciprocal_g = 2.0 * np.pi / lattice_a
    ks = np.linspace(-np.pi / lattice_a, np.pi / lattice_a, k_points)
    phases = np.linspace(0.0, 2.0 * np.pi, phase_points)
    modes = np.arange(-plane_wave_cutoff, plane_wave_cutoff + 1)
    n_modes = modes.size
    adjacent = np.arange(n_modes - 1)
    next_adjacent = np.arange(n_modes - 2)
    gaps = np.empty(phase_points)

    for phase_index, phase in enumerate(phases):
        gap_at_k = np.empty(k_points)
        long_forward = -0.5 * vl * np.exp(1j * phase)
        long_backward = np.conj(long_forward)
        for k_index, k in enumerate(ks):
            hamiltonian = np.diag(
                0.5 * (k + modes * reciprocal_g) ** 2
            ).astype(complex)
            # -Vl cos(G x + phi): Fourier couplings m <-> m+1.
            hamiltonian[adjacent + 1, adjacent] += long_forward
            hamiltonian[adjacent, adjacent + 1] += long_backward
            # -Vs cos(2 G x): Fourier couplings m <-> m+2.
            hamiltonian[next_adjacent + 2, next_adjacent] += -0.5 * vs
            hamiltonian[next_adjacent, next_adjacent + 2] += -0.5 * vs
            eigenvalues = np.linalg.eigvalsh(hamiltonian)
            gap_at_k[k_index] = eigenvalues[1] - eigenvalues[0]
        gaps[phase_index] = gap_at_k.min()

    return phases, gaps


def build_schedule(
    phases: np.ndarray,
    gaps: np.ndarray,
    schedule_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dense_phase = np.linspace(0.0, 2.0 * np.pi, schedule_points)
    dense_gap = np.interp(dense_phase, phases, gaps)
    inverse_gap_squared = 1.0 / dense_gap**2
    cumulative = np.concatenate(
        [
            [0.0],
            np.cumsum(
                0.5
                * (inverse_gap_squared[1:] + inverse_gap_squared[:-1])
                * np.diff(dense_phase)
            ),
        ]
    )
    normalized_time = cumulative / cumulative[-1]
    normalized_time[0] = 0.0
    normalized_time[-1] = 1.0
    return normalized_time, dense_phase, dense_gap, inverse_gap_squared


def write_schedule(
    path: Path,
    *,
    normalized_time: np.ndarray,
    phases: np.ndarray,
    gaps: np.ndarray,
    inverse_gap_squared: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as output:
        # Keep the header on the first line for csv.DictReader and
        # numpy.genfromtxt(names=True) compatibility. Provenance lives in this
        # generator, the convention note, and each run's config.txt.
        output.write("s,phi,g,inv_g2\n")
        for s, phase, gap, inv_gap2 in zip(
            normalized_time, phases, gaps, inverse_gap_squared
        ):
            output.write(
                f"{s:.17g},{phase:.17g},{gap:.17g},{inv_gap2:.17g}\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--lattice-a", type=float, default=8.0)
    parser.add_argument("--vs-over-Er", type=float, default=3.0)
    parser.add_argument("--vl-over-Er", type=float, default=3.0)
    # Match the full-depth grid reference's native phase grid exactly.
    parser.add_argument("--phase-points", type=int, default=241)
    parser.add_argument("--k-points", type=int, default=121)
    parser.add_argument("--plane-wave-cutoff", type=int, default=20)
    parser.add_argument("--schedule-points", type=int, default=6001)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if min(
        args.lattice_a,
        args.vs_over_Er,
        args.vl_over_Er,
        args.phase_points,
        args.k_points,
        args.plane_wave_cutoff,
        args.schedule_points,
    ) <= 0:
        raise ValueError("all depths, sizes, and point counts must be positive")
    if args.phase_points < 2 or args.k_points < 2 or args.schedule_points < 2:
        raise ValueError("phase, k, and schedule grids need at least two points")

    recoil_energy = 2.0 * np.pi**2 / args.lattice_a**2
    phases, gaps = band_gap_curve(
        lattice_a=args.lattice_a,
        vs=args.vs_over_Er * recoil_energy,
        vl=args.vl_over_Er * recoil_energy,
        phase_points=args.phase_points,
        k_points=args.k_points,
        plane_wave_cutoff=args.plane_wave_cutoff,
    )
    normalized_time, dense_phase, dense_gap, inverse_gap_squared = build_schedule(
        phases, gaps, args.schedule_points
    )

    if not np.all(np.diff(normalized_time) > 0.0):
        raise RuntimeError("generated normalized time is not strictly increasing")
    if not np.all(np.isfinite(dense_gap)) or dense_gap.min() <= 0.0:
        raise RuntimeError("generated band gap is nonpositive or nonfinite")

    write_schedule(
        args.output,
        normalized_time=normalized_time,
        phases=dense_phase,
        gaps=dense_gap,
        inverse_gap_squared=inverse_gap_squared,
    )
    gap_min_index = int(np.argmin(dense_gap))
    print(f"[write] {args.output} ({args.schedule_points} schedule points)")
    print(
        f"[full-depth] g_min/E_R={dense_gap[gap_min_index] / recoil_energy:.12g} "
        f"at phi/pi={dense_phase[gap_min_index] / np.pi:.12g}"
    )


if __name__ == "__main__":
    main()
