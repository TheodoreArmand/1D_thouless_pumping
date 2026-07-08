#!/usr/bin/env python3
"""
Split-step 2D grid reference for the adjacent-cell N=2 Lohse validation.

This script is intentionally standalone and conservative. Defaults match the
short C++ validation window; reduce --grid for quick local checks and increase
to --grid 1024 for the reference comparison.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DEFAULT_BASIS = REPO / "initial_state" / "Vs10Ers_Vl5Ers" / "initial_lohse_N2_K24.csv"
DEFAULT_SCHEDULE = HERE / "gap_adaptive_lohse_maincpp_schedule.csv"

A_LAT = 8.0
ER = 2 * np.pi**2 / A_LAT**2
VS, VL = 5.0 * ER, 2.5 * ER


def load_schedule(path: Path):
    ss, phis = [], []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("s,"):
            continue
        cols = line.split(",")
        ss.append(float(cols[0]))
        phis.append(float(cols[1]))
    return np.asarray(ss), np.asarray(phis)


def phi_at(t: float, period: float, ss: np.ndarray, phis: np.ndarray) -> float:
    return float(np.interp(np.clip(t / period, 0.0, 1.0), ss, phis))


def v_lattice(x: np.ndarray, phi: float) -> np.ndarray:
    return -VS * np.cos(4 * np.pi * x / A_LAT) - VL * np.cos(2 * np.pi * x / A_LAT + phi)


def load_n2_basis(path: Path):
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    stride = 1 + 4 + 4 + 2 + 1
    if len(rows) % stride != 0:
        raise RuntimeError(f"unexpected N=2 CSV length: {path}")
    basis = []
    for k in range(len(rows) // stride):
        off = stride * k
        u = complex(*map(float, rows[off].split(","))).real
        b00 = complex(*map(float, rows[off + 5].split(","))).real
        b11 = complex(*map(float, rows[off + 8].split(","))).real
        r0 = complex(*map(float, rows[off + 9].split(","))).real
        r1 = complex(*map(float, rows[off + 10].split(","))).real
        basis.append((u, b00, b11, r0, r1))
    return basis


def gaussian(x: np.ndarray, length: float, b: float, r: float) -> np.ndarray:
    d = (x - r + length / 2) % length - length / 2
    return np.exp(-b * d * d)


def initial_wavefunction(path: Path, x: np.ndarray, length: float):
    psi = np.zeros((x.size, x.size), dtype=np.complex128)
    cache = {}

    def get(b, r):
        key = (float(b), float(r))
        if key not in cache:
            cache[key] = gaussian(x, length, b, r)
        return cache[key]

    for u, b0, b1, r0, r1 in load_n2_basis(path):
        g0 = get(b0, r0)
        g1 = get(b1, r1)
        psi += u * (g0[:, None] * g1[None, :] + g1[:, None] * g0[None, :]) / np.sqrt(2.0)
    return psi


def observables(psi, x, dx, g_code, sigma):
    prob = np.abs(psi) ** 2
    norm = float(prob.sum() * dx * dx)
    x_mean = float(((x[:, None] + x[None, :]) * prob).sum() * dx * dx / norm)
    d = x[:, None] - x[None, :]
    v_gauss = float((g_code * np.exp(-(d * d) / (sigma * sigma)) * prob).sum() * dx * dx / norm)
    return norm, x_mean, v_gauss


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", type=Path, default=DEFAULT_BASIS)
    ap.add_argument("--schedule", type=Path, default=DEFAULT_SCHEDULE)
    ap.add_argument("--grid", type=int, default=512)
    ap.add_argument("--cells", type=int, default=16)
    ap.add_argument("--dt", type=float, default=2.0e-3)
    ap.add_argument("--tmax", type=float, default=12.0)
    ap.add_argument("--period", type=float, default=500.0 * np.pi)
    ap.add_argument("--sample-every", type=int, default=25)
    ap.add_argument("--g-over-Er", type=float, default=0.0)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--out", type=Path, default=HERE / "n2_grid_reference.npz")
    args = ap.parse_args()

    length = args.cells * A_LAT
    dx = length / args.grid
    x = np.arange(args.grid) * dx - length / 2
    p = 2 * np.pi * np.fft.fftfreq(args.grid, d=dx)
    p2 = p[:, None] ** 2 + p[None, :] ** 2
    g_code = args.g_over_Er * ER

    ss, phis = load_schedule(args.schedule)
    psi = initial_wavefunction(args.basis, x, length)
    psi /= np.sqrt((np.abs(psi) ** 2).sum() * dx * dx)

    t_values, norm_values, x_values, vint_values = [], [], [], []

    def sample(t):
        n, xm, vg = observables(psi, x, dx, g_code, args.sigma)
        t_values.append(t)
        norm_values.append(n)
        x_values.append(xm)
        vint_values.append(vg)

    sample(0.0)
    n_steps = int(np.ceil(args.tmax / args.dt))
    t = 0.0
    for step in range(1, n_steps + 1):
        dt = min(args.dt, args.tmax - t)
        if dt <= 0:
            break
        phi_mid = phi_at(t + 0.5 * dt, args.period, ss, phis)
        v = v_lattice(x, phi_mid)
        d = x[:, None] - x[None, :]
        v2 = v[:, None] + v[None, :] + g_code * np.exp(-(d * d) / (args.sigma * args.sigma))
        psi *= np.exp(-0.5j * dt * v2)
        psi = np.fft.ifft2(np.exp(-1j * dt * 0.5 * p2) * np.fft.fft2(psi))
        psi *= np.exp(-0.5j * dt * v2)
        t += dt
        if step % args.sample_every == 0 or t >= args.tmax - 1e-15:
            sample(t)

    np.savez(args.out,
             t=np.asarray(t_values),
             norm=np.asarray(norm_values),
             x_mean=np.asarray(x_values),
             V_gauss=np.asarray(vint_values),
             grid=args.grid,
             dt=args.dt,
             tmax=args.tmax,
             g_over_Er=args.g_over_Er,
             sigma=args.sigma)
    print(f"[write] {args.out} ({len(t_values)} samples, grid={args.grid})")


if __name__ == "__main__":
    main()
