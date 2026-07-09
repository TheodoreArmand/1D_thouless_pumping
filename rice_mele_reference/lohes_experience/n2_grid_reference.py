#!/usr/bin/env python3
"""
Split-step 2D grid reference for adjacent-cell N=2 pump validation.

This script is intentionally standalone and conservative. Defaults match the
short C++ validation window; reduce --grid for quick local checks and increase
to --grid 1024 for the reference comparison.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

try:
    import scipy.fft as scipy_fft
except ImportError:
    scipy_fft = None

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DEFAULT_BASIS = REPO / "initial_state" / "Vs10Ers_Vl5Ers" / "initial_lohse_N2_K24.csv"
DEFAULT_SCHEDULE = HERE / "gap_adaptive_lohse_maincpp_schedule.csv"

A_LAT = 8.0
ER = 2 * np.pi**2 / A_LAT**2
DEFAULT_VS_OVER_ER = 5.0
DEFAULT_VL_OVER_ER = 2.5


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


def recoil_energy(a_lat: float) -> float:
    return 2 * np.pi**2 / (a_lat * a_lat)


def v_lattice(x: np.ndarray, phi: float, a_lat: float, vs: float, vl: float) -> np.ndarray:
    return -vs * np.cos(4 * np.pi * x / a_lat) - vl * np.cos(2 * np.pi * x / a_lat + phi)


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


def observables(psi, x_sum, x_diff, dx, a_lat, g_code, gauss_kernel):
    prob = np.abs(psi) ** 2
    norm = float(prob.sum() * dx * dx)
    x_mean = float((x_sum * prob).sum() * dx * dx / norm)
    r12_sq = float((x_diff * x_diff * prob).sum() * dx * dx / norm)
    v_gauss = float((g_code * gauss_kernel * prob).sum() * dx * dx / norm)
    return norm, x_mean, x_mean / a_lat, r12_sq, v_gauss


def make_fft_pair(workers: int):
    workers = max(1, int(workers))
    if scipy_fft is None:
        if workers != 1:
            print("[warn] scipy.fft unavailable; falling back to numpy.fft with one worker", flush=True)
        return np.fft.fft2, np.fft.ifft2, "numpy.fft", 1

    def fft2(a):
        return scipy_fft.fft2(a, workers=workers)

    def ifft2(a):
        return scipy_fft.ifft2(a, workers=workers)

    return fft2, ifft2, "scipy.fft", workers


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", type=Path, default=DEFAULT_BASIS)
    ap.add_argument("--schedule", type=Path, default=DEFAULT_SCHEDULE)
    ap.add_argument("--grid", type=int, default=512)
    ap.add_argument("--cells", type=int, default=16)
    ap.add_argument("--a-lat", type=float, default=A_LAT)
    ap.add_argument("--vs-over-Er", type=float, default=DEFAULT_VS_OVER_ER)
    ap.add_argument("--vl-over-Er", type=float, default=DEFAULT_VL_OVER_ER)
    ap.add_argument("--dt", type=float, default=2.0e-3)
    ap.add_argument("--tmax", type=float, default=12.0)
    ap.add_argument("--period", type=float, default=500.0 * np.pi)
    ap.add_argument("--sample-every", type=int, default=25)
    ap.add_argument("--progress-every", type=int, default=250)
    ap.add_argument("--fft-workers", type=int, default=1)
    ap.add_argument("--g-over-Er", type=float, default=0.0)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--out", type=Path, default=HERE / "n2_grid_reference.npz")
    ap.add_argument("--csv", type=Path, default=None)
    args = ap.parse_args()

    er = recoil_energy(args.a_lat)
    vs_code = args.vs_over_Er * er
    vl_code = args.vl_over_Er * er
    g_code = args.g_over_Er * er

    length = args.cells * args.a_lat
    dx = length / args.grid
    x = np.arange(args.grid) * dx - length / 2
    p = 2 * np.pi * np.fft.fftfreq(args.grid, d=dx)
    p2 = p[:, None] ** 2 + p[None, :] ** 2
    x_sum = x[:, None] + x[None, :]
    x_diff = x[:, None] - x[None, :]
    gauss_kernel = np.exp(-(x_diff * x_diff) / (args.sigma * args.sigma))
    kinetic_phase_full = np.exp(-1j * args.dt * 0.5 * p2)
    fft2, ifft2, fft_backend, fft_workers = make_fft_pair(args.fft_workers)
    print(
        f"[config] grid={args.grid} dt={args.dt} tmax={args.tmax} "
        f"Vs/Er={args.vs_over_Er} Vl/Er={args.vl_over_Er} "
        f"g/Er={args.g_over_Er} fft_backend={fft_backend} "
        f"fft_workers={fft_workers}",
        flush=True,
    )

    ss, phis = load_schedule(args.schedule)
    psi = initial_wavefunction(args.basis, x, length)
    psi /= np.sqrt((np.abs(psi) ** 2).sum() * dx * dx)

    t_values, phi_values = [], []
    norm_values, x_values, pcell_values = [], [], []
    dp_values, r12_sq_values, vint_values = [], [], []

    def sample(t):
        phi = phi_at(t, args.period, ss, phis)
        n, xm, pcell, r12_sq, vg = observables(
            psi, x_sum, x_diff, dx, args.a_lat, g_code, gauss_kernel)
        t_values.append(t)
        phi_values.append(phi)
        norm_values.append(n)
        x_values.append(xm)
        pcell_values.append(pcell)
        dp_values.append(pcell - pcell_values[0])
        r12_sq_values.append(r12_sq)
        vint_values.append(vg)

    sample(0.0)
    n_steps = int(np.ceil(args.tmax / args.dt))
    t = 0.0
    start = time.perf_counter()
    for step in range(1, n_steps + 1):
        dt = min(args.dt, args.tmax - t)
        if dt <= 0:
            break
        phi_mid = phi_at(t + 0.5 * dt, args.period, ss, phis)
        v = v_lattice(x, phi_mid, args.a_lat, vs_code, vl_code)
        v2 = v[:, None] + v[None, :] + g_code * gauss_kernel
        psi *= np.exp(-0.5j * dt * v2)
        kinetic_phase = kinetic_phase_full if dt == args.dt else np.exp(-1j * dt * 0.5 * p2)
        psi = ifft2(kinetic_phase * fft2(psi))
        psi *= np.exp(-0.5j * dt * v2)
        t += dt
        if step % args.sample_every == 0 or t >= args.tmax - 1e-15:
            sample(t)
        if args.progress_every > 0 and (
                step % args.progress_every == 0 or t >= args.tmax - 1e-15):
            wall = time.perf_counter() - start
            s_per_step = wall / step
            eta = s_per_step * max(n_steps - step, 0)
            print(
                f"progress grid N=2 step {step}/{n_steps} "
                f"{100.0 * step / n_steps:.1f}% t={t:.5f} eta_s={eta:.1f} "
                f"P={pcell_values[-1]:.8f} dP={dp_values[-1]:.8f} "
                f"r12={np.sqrt(max(r12_sq_values[-1], 0.0)):.8f} "
                f"Vg={vint_values[-1]:.8g} norm={norm_values[-1]:.8f}",
                flush=True,
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out,
             t=np.asarray(t_values),
             phi=np.asarray(phi_values),
             norm=np.asarray(norm_values),
             x_mean=np.asarray(x_values),
             polarization_cell=np.asarray(pcell_values),
             delta_polarization=np.asarray(dp_values),
             r12_sq=np.asarray(r12_sq_values),
             V_gauss=np.asarray(vint_values),
             grid=args.grid,
             cells=args.cells,
             a_lat=args.a_lat,
             dt=args.dt,
             tmax=args.tmax,
             period=args.period,
             vs_over_Er=args.vs_over_Er,
             vl_over_Er=args.vl_over_Er,
             g_over_Er=args.g_over_Er,
             sigma=args.sigma,
             fft_backend=fft_backend,
             fft_workers=fft_workers,
             basis=str(args.basis),
             schedule=str(args.schedule))
    csv_path = args.csv if args.csv is not None else args.out.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_data = np.column_stack([
        np.asarray(t_values),
        np.asarray(phi_values),
        np.asarray(norm_values),
        np.asarray(x_values),
        np.asarray(pcell_values),
        np.asarray(dp_values),
        np.asarray(r12_sq_values),
        np.sqrt(np.maximum(np.asarray(r12_sq_values), 0.0)),
        np.asarray(vint_values),
    ])
    np.savetxt(
        csv_path,
        csv_data,
        delimiter=",",
        header="t,phi,norm,x_mean,polarization_cell,delta_polarization,r12_sq,r12_rms,V_gauss",
        comments="",
    )
    print(f"[write] {args.out} ({len(t_values)} samples, grid={args.grid})")
    print(f"[write] {csv_path}")


if __name__ == "__main__":
    main()
