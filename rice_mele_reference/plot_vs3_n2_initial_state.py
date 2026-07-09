#!/usr/bin/env python3
"""Plot the Vs3/Vl3 N=2 initial ECG state on the lattice potential."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
REF_ROOT = REPO / "rice_mele_reference"
BASIS = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N2_K24.csv"

A_LAT = 8.0
CELLS = 16
LENGTH = CELLS * A_LAT
GRID = 1024
VS_OVER_ER = 3.0
VL_OVER_ER = 3.0
SIGMA_GAUSS = 1.0
WINDOW = (-16.0, 24.0)


def load_n2_basis(path: Path):
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    stride = 12
    if len(rows) % stride != 0:
        raise RuntimeError(f"unexpected N=2 CSV length: {path}")
    terms = []
    for k in range(len(rows) // stride):
        off = stride * k
        u = complex(*map(float, rows[off].split(","))).real
        b00 = complex(*map(float, rows[off + 5].split(","))).real
        b11 = complex(*map(float, rows[off + 8].split(","))).real
        r0 = complex(*map(float, rows[off + 9].split(","))).real
        r1 = complex(*map(float, rows[off + 10].split(","))).real
        terms.append((u, b00, b11, r0, r1))
    return terms


def gaussian_on_ring(x: np.ndarray, length: float, b: float, r: float) -> np.ndarray:
    d = (x - r + length / 2.0) % length - length / 2.0
    return np.exp(-b * d * d)


def initial_wavefunction(terms, x: np.ndarray) -> np.ndarray:
    psi = np.zeros((x.size, x.size), dtype=np.complex128)
    cache = {}

    def get(b: float, r: float) -> np.ndarray:
        key = (float(b), float(r))
        if key not in cache:
            cache[key] = gaussian_on_ring(x, LENGTH, b, r)
        return cache[key]

    for u, b0, b1, r0, r1 in terms:
        g0 = get(b0, r0)
        g1 = get(b1, r1)
        psi += u * (g0[:, None] * g1[None, :] + g1[:, None] * g0[None, :]) / np.sqrt(2.0)
    return psi


def lattice_potential_over_er(x: np.ndarray, phi: float = 0.0) -> np.ndarray:
    return (-VS_OVER_ER * np.cos(4.0 * np.pi * x / A_LAT)
            -VL_OVER_ER * np.cos(2.0 * np.pi * x / A_LAT + phi))


def write_summary(path: Path, *, mode: str, g_over_er: float, norm: float,
                  x_mean: float, polarization: float, r12_rms: float,
                  v_gauss_over_er: float) -> None:
    path.write_text(
        "\n".join([
            f"mode={mode}",
            f"basis={BASIS.relative_to(REPO)}",
            f"grid={GRID}",
            f"cells={CELLS}",
            f"a_lat={A_LAT}",
            f"Vs_over_Er={VS_OVER_ER}",
            f"Vl_over_Er={VL_OVER_ER}",
            f"g_over_Er={g_over_er}",
            "initial_state_same_as_ecg_basis=yes",
            "ecg_basis_form=real coefficients, A=0, diagonal real B",
            "grid_reconstruction=sum_k u_k [g0(x0) g1(x1) + g1(x0) g0(x1)]/sqrt(2)",
            f"norm={norm:.16g}",
            f"x_mean_sum={x_mean:.16g}",
            f"polarization_cell={polarization:.16g}",
            f"r12_rms={r12_rms:.16g}",
            f"V_gauss_over_Er={v_gauss_over_er:.16g}",
            "",
        ])
    )


def make_plot(out_dir: Path, *, mode: str, g_over_er: float,
              x: np.ndarray, dx: float, psi: np.ndarray) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    prob = np.abs(psi) ** 2
    norm = float(prob.sum() * dx * dx)
    prob /= norm

    x_sum = x[:, None] + x[None, :]
    x_diff = x[:, None] - x[None, :]
    x_mean = float((x_sum * prob).sum() * dx * dx)
    polarization = x_mean / A_LAT
    r12_sq = float((x_diff * x_diff * prob).sum() * dx * dx)
    r12_rms = float(np.sqrt(max(r12_sq, 0.0)))
    pair_kernel = np.exp(-(x_diff * x_diff) / (SIGMA_GAUSS * SIGMA_GAUSS))
    v_gauss_over_er = float((g_over_er * pair_kernel * prob).sum() * dx * dx)

    rho0 = prob.sum(axis=1) * dx
    rho1 = prob.sum(axis=0) * dx
    rho_total = rho0 + rho1

    mask = (x >= WINDOW[0]) & (x <= WINDOW[1])
    xw = x[mask]
    rho_w = rho_total[mask]
    v_w = lattice_potential_over_er(xw)

    prob_crop = prob[np.ix_(mask, mask)]
    prob_rel = prob_crop / prob_crop.max()
    v1 = lattice_potential_over_er(xw)
    pair_crop = g_over_er * np.exp(-((xw[:, None] - xw[None, :]) ** 2)
                                    / (SIGMA_GAUSS * SIGMA_GAUSS))
    v2_over_er = v1[:, None] + v1[None, :] + pair_crop

    fig = plt.figure(figsize=(9.6, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25])

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(xw, v_w, color="#222222", lw=1.8, label="lattice V(x)/Er")
    for well in np.arange(-64.0, 65.0, A_LAT):
        if WINDOW[0] <= well <= WINDOW[1]:
            ax.axvline(well, color="#9a9a9a", lw=0.8, ls=":", zorder=0)
    ax.set_xlim(WINDOW)
    ax.set_ylabel("V(x) / Er")
    ax.set_xlabel("x")
    ax.grid(True, color="#d6d6d6", lw=0.6, alpha=0.8)

    ax_rho = ax.twinx()
    ax_rho.fill_between(xw, rho_w, color="#2f80ed", alpha=0.28,
                        label="one-body density")
    ax_rho.plot(xw, rho_w, color="#2f80ed", lw=1.2)
    ax_rho.set_ylabel("rho1(x), integral = 2")

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax_rho.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc="upper right", frameon=False)
    ax.set_title(f"Vs=Vl=3Er N=2 initial state on lattice wells ({mode})")

    ax2 = fig.add_subplot(gs[1, 0])
    extent = [xw[0], xw[-1], xw[0], xw[-1]]
    im = ax2.imshow(
        np.log10(prob_rel.T + 1.0e-8),
        origin="lower",
        extent=extent,
        cmap="magma",
        vmin=-8.0,
        vmax=0.0,
        interpolation="nearest",
        aspect="equal",
    )
    levels = np.linspace(np.nanmin(v2_over_er), min(2.0, np.nanmax(v2_over_er)), 9)
    ax2.contour(xw, xw, v2_over_er.T, levels=levels, colors="#f2f2f2",
                linewidths=0.55, alpha=0.65)
    for well in np.arange(-64.0, 65.0, A_LAT):
        if WINDOW[0] <= well <= WINDOW[1]:
            ax2.axvline(well, color="#ffffff", lw=0.45, ls=":", alpha=0.45)
            ax2.axhline(well, color="#ffffff", lw=0.45, ls=":", alpha=0.45)
    ax2.set_xlabel("x0")
    ax2.set_ylabel("x1")
    ax2.set_title("log10 normalized pair density, with two-particle potential contours")
    cbar = fig.colorbar(im, ax=ax2, shrink=0.92)
    cbar.set_label("log10(|psi|^2 / max)")

    fig.savefig(out_dir / "initial_state_on_lattice.png", dpi=180)
    plt.close(fig)

    write_summary(
        out_dir / "initial_state_summary.txt",
        mode=mode,
        g_over_er=g_over_er,
        norm=norm,
        x_mean=x_mean,
        polarization=polarization,
        r12_rms=r12_rms,
        v_gauss_over_er=v_gauss_over_er,
    )


def main() -> None:
    terms = load_n2_basis(BASIS)
    dx = LENGTH / GRID
    x = np.arange(GRID) * dx - LENGTH / 2.0
    psi = initial_wavefunction(terms, x)
    make_plot(REF_ROOT / "Vs3Vl3_3_3_2noin", mode="free/no interaction",
              g_over_er=0.0, x=x, dx=dx, psi=psi)
    make_plot(REF_ROOT / "Vs3Vl3_3_3_2gauss", mode="gaussian pair",
              g_over_er=0.3, x=x, dx=dx, psi=psi)
    print("wrote:")
    print(REF_ROOT / "Vs3Vl3_3_3_2noin" / "initial_state_on_lattice.png")
    print(REF_ROOT / "Vs3Vl3_3_3_2gauss" / "initial_state_on_lattice.png")


if __name__ == "__main__":
    main()
