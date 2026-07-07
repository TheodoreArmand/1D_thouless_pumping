#!/usr/bin/env python3
"""
make_lohse_reference.py
=======================

Self-contained reference computation for the Lohse et al. (2016) Thouless pump
[Nature Physics 12, 350; ../../reference/lohes_nature_physics.pdf], at the
paper's experimental depths.

Continuum Hamiltonian (paper convention)
----------------------------------------
    H(phi) = -(hbar^2/2m) d^2/dx^2 + V(x, phi)
    V(x, phi) = Vs sin^2(pi x/ds + pi/2) + Vl sin^2(pi x/dl - phi/2),  dl = 2 ds
              = Vs cos^2(pi x/ds)       + Vl sin^2(pi x/dl - phi/2)

Units: ds = 1, E_rs = hbar^2 k_s^2/2m = 1 (k_s = pi), so hbar^2/2m = 1/pi^2.
Depths: Vs = 10 E_rs, Vl = 20 E_rl = 5 E_rs (E_rl = E_rs/4).

Outputs (this directory):
  lohse_reference_data.npz        all curves (bands, gap, J1/J2/Delta, scans)
  lohse_reference_summary.json    headline numbers for the HTML report
  gap_adaptive_lohse_schedule.csv           s,phi  (paper convention)
  gap_adaptive_lohse_maincpp_schedule.csv   s,phi  (main.cpp convention phi^m = 3pi/2 - phi)
  figs/fig1_potential.png         V(x,phi) snapshots
  figs/fig2_bands.png             continuum bands at phi = 0, pi/2
  figs/fig3_rice_mele_params.png  J1, J2, Delta, band gaps vs phi
  figs/fig4_pump_path.png         (J1-J2, Delta) pump loop
  figs/fig5_gap_adaptive_pump.png 3-panel pump demo at T = 400 pi (like the
                                  Vs3Vl3 reference figure)
  figs/fig6_adiabatic_time_scan.png  final P0 and transport vs T, both schedules
  figs/fig7_uniform_T2000pi.png   uniform pump at T = 2000 pi (adiabatic)

Cross-checks against ../../../periodic/rice_mele_reference (old project):
  J1(0) = 0.0607, J2(0) = 0.0051, Delta_max = 4.519, g12_pump = 0.1113 E_rs,
  eta = max |<1|dH/dphi|0>|/g^2 = 179.5 / E_rs.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ecg_rice_mele_mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figs"
FIGDIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# model constants (ds = 1, E_rs = 1)
# ---------------------------------------------------------------------------
VS = 10.0                  # short-lattice depth [E_rs]
VL = 5.0                   # long-lattice depth  [E_rs] (= 20 E_rl)
DS, DL = 1.0, 2.0
HB2 = 1.0 / np.pi**2       # hbar^2/2m

# physical conversion (87Rb, lambda_l = 1534 nm => ds = 383.5 nm)
H_PLANCK = 6.62607015e-34
HBAR = H_PLANCK / (2 * np.pi)
M_RB87 = 86.909180527 * 1.66053906892e-27
LAMBDA_S = 767.0e-9
E_RS_J = H_PLANCK**2 / (2 * M_RB87 * LAMBDA_S**2)
E_RS_HZ = E_RS_J / H_PLANCK                # E_rs/h in Hz
T_UNIT_S = HBAR / E_RS_J                   # one time unit (hbar/E_rs) in s
# ECG main.cpp units: hbar=m=1, a=8 -> E_rs^code = (pi/2)^2/2 = pi^2/8
T_CODE_PER_UNIT = 1.0 / (np.pi**2 / 8.0)   # inverse: T_code = T_ours / (pi^2/8)?  see note
# NOTE: E_rs expressed in code units is pi^2/8 ~ 1.2337, so one of our time
# units (hbar/E_rs) equals 1/1.2337 = 0.8106 code time units.
T_CODE_FACTOR = 8.0 / np.pi**2             # multiply T[hbar/E_rs] by this -> code units

# ---------------------------------------------------------------------------
# house style (matches the Vs3Vl3 reference figures)
# ---------------------------------------------------------------------------
INK, MUTED, GRIDC, AXISC = "#0b0b0b", "#898781", "#e1e0d9", "#c3c2b7"
C_J1, C_J2, C_DELTA = "#2a78d6", "#eda100", "#1baf7a"
C_GAP, C_GAP23 = "#4a3aa7", "#eb6834"
C_ADA, C_UNI, C_SCH = "#1baf7a", "#c0392b", "#2a78d6"
plt.rcParams.update({
    "font.size": 11.0, "font.family": "DejaVu Sans", "text.color": INK,
    "axes.edgecolor": AXISC, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def potential(x, phi):
    return VS * np.cos(np.pi * x / DS) ** 2 + VL * np.sin(np.pi * x / DL - phi / 2.0) ** 2


# ===========================================================================
# 1. plane-wave (Bloch) bands: basis e^{i(kappa+m) pi x}, kappa in [-1/2, 1/2]
# ===========================================================================
CUT = 25
_M = np.arange(-CUT, CUT + 1)
_NB = _M.size


def bloch_h_batch(kappas, phi):
    """Batched Bloch Hamiltonians, shape (nk, n, n), energies in E_rs."""
    n = _NB
    H = np.zeros((len(kappas), n, n), dtype=complex)
    i1 = np.arange(n - 1)
    i2 = np.arange(n - 2)
    H[:, i1 + 1, i1] = -(VL / 4.0) * np.exp(-1j * phi)
    H[:, i1, i1 + 1] = -(VL / 4.0) * np.exp(+1j * phi)
    H[:, i2 + 2, i2] = VS / 4.0
    H[:, i2, i2 + 2] = VS / 4.0
    idx = np.arange(n)
    H[:, idx, idx] = (np.asarray(kappas)[:, None] + _M) ** 2
    return H


def dh_dphi(phi):
    n = _NB
    D = np.zeros((n, n), dtype=complex)
    i1 = np.arange(n - 1)
    D[i1 + 1, i1] = +1j * (VL / 4.0) * np.exp(-1j * phi)
    D[i1, i1 + 1] = -1j * (VL / 4.0) * np.exp(+1j * phi)
    return D


def band_fit_AB(kap, E01):
    """Fit exact (gap/2)^2 to A + B cos(k dl): A=(Delta/2)^2+J1^2+J2^2, B=2 J1 J2."""
    y = ((E01[:, 1] - E01[:, 0]) / 2.0) ** 2
    c = np.cos(kap * np.pi * DL)
    coef, *_ = np.linalg.lstsq(np.vstack([np.ones_like(c), c]).T, y, rcond=None)
    return float(coef[0]), float(coef[1])


def compute_band_curves(nphi=241, nk=401):
    phis = np.linspace(0.0, 2 * np.pi, nphi)
    kap = np.linspace(-0.5, 0.5, nk)
    E12min = np.empty(nphi)
    E23min = np.empty(nphi)
    E0_phi = np.empty(nphi)
    A_phi = np.empty(nphi)
    B_phi = np.empty(nphi)
    rm_err = np.empty(nphi)
    eta_phi = np.empty(nphi)      # max_k |<1|dH/dphi|0>| / g^2
    for i, phi in enumerate(phis):
        Hb = bloch_h_batch(kap, phi)
        w, v = np.linalg.eigh(Hb)                     # batched
        E = w[:, :3]
        g12 = E[:, 1] - E[:, 0]
        E12min[i] = g12.min()
        E23min[i] = (E[:, 2] - E[:, 1]).min()
        E0_phi[i] = 0.5 * (E[:, 0] + E[:, 1]).mean()
        A, B = band_fit_AB(kap, E[:, :2])
        A_phi[i], B_phi[i] = A, B
        half = np.sqrt(np.maximum(A + B * np.cos(kap * np.pi * DL), 0.0))
        rm_err[i] = max(np.max(np.abs(E[:, 0] - (E0_phi[i] - half))),
                        np.max(np.abs(E[:, 1] - (E0_phi[i] + half))))
        D = dh_dphi(phi)
        cpl = np.abs(np.einsum("ki,ij,kj->k", v[:, :, 1].conj(), D, v[:, :, 0]))
        eta_phi[i] = np.max(cpl / g12**2)
    return dict(phis=phis, E12min=E12min, E23min=E23min, E0_phi=E0_phi,
                A_phi=A_phi, B_phi=B_phi, rm_err=rm_err, eta_phi=eta_phi)


# ===========================================================================
# 2. Wannier projection: J1, J2, Delta as <w|H|w'> matrix elements
# ===========================================================================
def _build_kinetic(M, L):
    p = 2 * np.pi * np.fft.fftfreq(M, d=L / M)
    F = np.fft.fft(np.eye(M), axis=0)
    return np.fft.ifft((HB2 * p**2)[:, None] * F, axis=0).real


class WannierExtractor:
    def __init__(self, N=16, ppc=48):
        self.N, self.ppc = N, ppc
        self.L, self.M = DL * N, ppc * N
        self.x = np.arange(self.M) * (self.L / self.M)
        self.T = _build_kinetic(self.M, self.L)

    def params(self, phi):
        """J1 (even-integer barrier, strong at phi=0), J2 (odd barrier),
        Delta = eps_A - eps_B with A the (x-0.5)-even sublattice."""
        H = self.T + np.diag(potential(self.x, phi))
        E, psi = np.linalg.eigh(H)
        occ = 2 * self.N
        Eo, P = E[:occ], psi[:, :occ]
        Xp = P.conj().T @ (self.x[:, None] * P)
        Xp = 0.5 * (Xp + Xp.conj().T)
        centers, C = np.linalg.eigh(Xp)
        HW = (C.conj().T @ np.diag(Eo) @ C).real
        order = np.argsort(centers)
        centers = centers[order]
        HW = HW[np.ix_(order, order)]
        lo, hi = 6, 2 * self.N - 6
        onsite = np.diag(HW)
        j_even, j_odd, epsA, epsB = [], [], [], []
        for a in range(lo, hi):
            mid = 0.5 * (centers[a] + centers[a + 1])
            J = -HW[a, a + 1]
            (j_even if round(mid) % 2 == 0 else j_odd).append(abs(J))
        for a in range(lo, hi + 1):
            (epsA if round(centers[a] - 0.5) % 2 == 0 else epsB).append(onsite[a])
        return (float(np.mean(j_even)), float(np.mean(j_odd)),
                float(np.mean(epsA) - np.mean(epsB)))


def compute_wannier_curves(nphi_w=121):
    wex = WannierExtractor(N=16, ppc=48)
    phis_w = np.linspace(0.0, 2 * np.pi, nphi_w)
    J1 = np.empty(nphi_w)
    J2 = np.empty(nphi_w)
    Delta = np.empty(nphi_w)
    for i, phi in enumerate(phis_w):
        J1[i], J2[i], Delta[i] = wex.params(phi)
        if i % 20 == 0:
            print(f"    wannier phi={phi/np.pi:4.2f}pi  J1={J1[i]:.5f} "
                  f"J2={J2[i]:.5f} Delta={Delta[i]:+.4f}", flush=True)
    idx_half = int(np.argmin(np.abs(phis_w - np.pi / 2)))
    if Delta[idx_half] < 0:
        Delta = -Delta
    return dict(phis_w=phis_w, J1=J1, J2=J2, Delta=Delta)


# ===========================================================================
# 3. real-time pump dynamics (split-step Fourier on a ring of N cells)
# ===========================================================================
NRING, PPC = 20, 48
LRING = DL * NRING
MPTS = PPC * NRING
DX = LRING / MPTS
XG = np.arange(MPTS) * DX
PG = 2 * np.pi * np.fft.fftfreq(MPTS, d=DX)

NS = 61
SPH = np.linspace(0, 2 * np.pi, NS)

_TK = _build_kinetic(MPTS, LRING)
QB = None          # instantaneous band-0 subspaces at SPH
W0 = None          # initial band-0 Wannier packet at phi=0
X0 = None


def init_dynamics():
    global QB, W0, X0
    print("    building band-0 subspaces / Wannier packet ...", flush=True)
    QB = np.array([np.linalg.eigh(_TK + np.diag(potential(XG, ph)))[1][:, :NRING]
                   for ph in SPH])
    B0 = QB[0]
    Xp = B0.conj().T @ (XG[:, None] * B0)
    Xp = 0.5 * (Xp + Xp.conj().T)
    cc, CC = np.linalg.eigh(Xp)
    Wf = B0 @ CC
    j = int(np.argmin(np.abs(cc - LRING / 2)))
    W0 = Wf[:, j].astype(complex)
    W0 /= np.sqrt(np.sum(np.abs(W0) ** 2) * DX)
    X0 = float(np.sum(np.abs(W0) ** 2 * XG) * DX)
    assert abs(DX * np.sum(np.abs(B0.conj().T @ W0) ** 2) - 1.0) < 1e-6


# gap-adaptive schedule dphi/dt = C g(phi)^2 built from the computed gap curve
PHIG = None
TCUM_UNIT = None   # cumulative time for T=1 (normalized), on PHIG


def build_schedule(phis, gap):
    global PHIG, TCUM_UNIT
    PHIG = np.linspace(0, 2 * np.pi, 6001)
    ginv2 = 1.0 / np.interp(PHIG, phis, gap) ** 2
    tc = np.concatenate([[0], np.cumsum(0.5 * (ginv2[1:] + ginv2[:-1]) * np.diff(PHIG))])
    TCUM_UNIT = tc / tc[-1]     # s(phi) in [0,1]


def phi_of_t_factory(Ttot, mode):
    if mode == "adaptive":
        return lambda t: np.interp(t / Ttot, TCUM_UNIT, PHIG)
    return lambda t: 2 * np.pi * t / Ttot


def run_pump(Ttot, mode, dt=0.02, traces=True):
    """One pump cycle phi: 0 -> 2pi. Returns final transport fraction f,
    final band-0 population P0, and (optionally) traces."""
    phit = phi_of_t_factory(Ttot, mode)
    Kin = np.exp(-1j * HB2 * PG**2 * dt)
    ns = int(round(Ttot / dt))
    psi = W0.copy()
    de = max(1, ns // 4000)
    tr_phi, tr_x = [], []
    p0 = np.full(NS, np.nan)
    ti = 0
    for s in range(ns + 1):
        ph_now = phit(s * dt)
        if traces and s % de == 0:
            tr_phi.append(ph_now)
            tr_x.append((np.sum(np.abs(psi) ** 2 * XG) * DX - X0) / DL)
        while ti < NS and ph_now >= SPH[ti] - 1e-9:
            p0[ti] = DX * np.sum(np.abs(QB[ti].conj().T @ psi) ** 2)
            ti += 1
        if s == ns:
            break
        pm = phit((s + 0.5) * dt)
        Vh = np.exp(-1j * potential(XG, pm) * 0.5 * dt)
        psi = Vh * psi
        psi = np.fft.ifft(Kin * np.fft.fft(psi))
        psi = Vh * psi
    f = (np.sum(np.abs(psi) ** 2 * XG) * DX - X0) / DL
    P0f = DX * np.sum(np.abs(QB[0].conj().T @ psi) ** 2)
    return dict(f=float(f), P0=float(P0f), tr_phi=np.array(tr_phi),
                tr_x=np.array(tr_x), p0=p0)


# ===========================================================================
# figures
# ===========================================================================
def fmt_pi_axis(ax, xmax=2 * np.pi):
    ax.set_xlim(0, xmax)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])


def fig_potential():
    x = np.linspace(0, 4, 1200)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 5.6), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (phi, lab) in zip(axes.flat,
                              [(0.0, "0"), (np.pi / 2, r"\pi/2"),
                               (np.pi, r"\pi"), (3 * np.pi / 2, r"3\pi/2")]):
        ax.plot(x, potential(x, phi), color=C_J1, lw=2.0)
        ax.set_title(rf"$\varphi = {lab}$")
        ax.axhline(0, color=AXISC, lw=0.8)
        ax.set_ylim(-0.6, VS + VL + 0.8)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$x / d_s$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$V(x,\varphi)\ [E_{r,s}]$")
    fig.suptitle(r"Lohse superlattice  $V(x,\varphi)=V_s\cos^2(\pi x/d_s)+V_l\sin^2(\pi x/d_l-\varphi/2)$"
                 "\n" r"$V_s=10\,E_{r,s},\ V_l=20\,E_{r,l}=5\,E_{r,s}$", fontsize=12)
    return fig


def fig_bands():
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4), sharey=True,
                             constrained_layout=True)
    cols = [C_J1, C_DELTA, C_GAP23, MUTED]
    labs = ["band 1", "band 2", "band 3", "band 4"]
    kap = np.linspace(-0.5, 0.5, 241)
    for ax, (phi, lab) in zip(axes, [(0.0, "0"), (np.pi / 2, r"\pi/2")]):
        E = np.linalg.eigvalsh(bloch_h_batch(kap, phi))[:, :4]
        for b in range(4):
            ax.plot(kap, E[:, b], color=cols[b], lw=2.2,
                    label=labs[b] if ax is axes[1] else None)
        ax.fill_between(kap, E[:, 0], E[:, 1], color=C_J1, alpha=0.08)
        ax.set_title(rf"$\varphi = {lab}$")
        ax.set_xlabel(r"$k\,d_l/2\pi$")
        ax.set_xlim(-0.5, 0.5)
        ax.set_xticks([-0.5, 0, 0.5])
    axes[0].set_ylabel(r"$E_n(k)\ [E_{r,s}]$")
    axes[1].legend(loc="center right", fontsize=9)
    axes[0].annotate("Rice–Mele\ntwo-band manifold", xy=(0.0, -1.7),
                     xytext=(0.0, 0.6), ha="center", fontsize=9, color=C_J1,
                     arrowprops=dict(arrowstyle="->", color=C_J1, lw=1.1))
    fig.suptitle("Exact continuum bands (plane-wave diagonalization)", fontsize=13)
    return fig


def fig_rm_params(d):
    fig, axes = plt.subplots(3, 1, figsize=(8.6, 9.0), sharex=True,
                             constrained_layout=True)
    ax = axes[0]
    ax.plot(d["phis_w"], d["J1"], color=C_J1, lw=2.4, label=r"$J_1$ (intracell)")
    ax.plot(d["phis_w"], d["J2"], color=C_J2, lw=2.4, label=r"$J_2$ (intercell)")
    ax.plot(d["phis_w"], d["J1"] - d["J2"], color=MUTED, lw=1.6, ls="--",
            label=r"$J_1-J_2$")
    ax.axhline(0, color=AXISC, lw=0.8)
    ax.set_ylabel(r"hopping $[E_{r,s}]$")
    ax.legend(ncol=3, loc="upper center", fontsize=9.5)
    ax.set_title(r"(a) hoppings — cross at $\varphi=\pi/2,\,3\pi/2$ "
                 r"($J_1{=}J_2$, tilt $\Delta$ maximal)")
    ax = axes[1]
    ax.plot(d["phis_w"], d["Delta"], color=C_DELTA, lw=2.4)
    ax.axhline(0, color=AXISC, lw=0.8)
    ax.set_ylabel(r"$\Delta(\varphi)\ [E_{r,s}]$")
    ax.set_title(r"(b) sublattice offset $\Delta=\epsilon_A-\epsilon_B$ "
                 rf"— amplitude ${d['Delta_max']:.2f}\,E_{{r,s}}$")
    ax = axes[2]
    ax.plot(d["phis"], d["E12min"], color=C_GAP, lw=2.4, label=r"1–2 gap $\min_k(E_2-E_1)$")
    ax.plot(d["phis"], d["E23min"], color=C_GAP23, lw=2.0, label=r"2–3 gap $\min_k(E_3-E_2)$")
    ax.axhline(d["g12_pump"], color=C_GAP, lw=0.9, ls=":")
    ax.text(0.15, d["g12_pump"] * 1.12, rf"$g_{{\rm pump}}={d['g12_pump']:.3f}\,E_{{r,s}}$",
            color=C_GAP, fontsize=9)
    ax.set_ylabel(r"band gap $[E_{r,s}]$")
    ax.set_yscale("log")
    ax.legend(loc="center right", fontsize=9.5)
    ax.set_title(r"(c) minimum band gaps — 1–2 gap dips at $\varphi=0,\pi$ "
                 "(the pump bottleneck)")
    ax.set_xlabel(r"pump phase $\varphi$")
    fmt_pi_axis(ax)
    return fig


def fig_pump_path(d):
    x = d["J1"] - d["J2"]
    y = d["Delta"]
    fig, ax = plt.subplots(figsize=(6.4, 6.0), constrained_layout=True)
    pts = np.array([x, y]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap="viridis", linewidth=3.0)
    lc.set_array(d["phis_w"][:-1] / np.pi)
    ax.add_collection(lc)
    ax.plot([0], [0], "x", color=C_UNI, ms=13, mew=3,
            label=r"degeneracy $J_1{=}J_2,\ \Delta{=}0$")
    for f in (0.12, 0.62):
        i = int(f * len(x))
        ax.annotate("", xy=(x[i + 1], y[i + 1]), xytext=(x[i], y[i]),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0))
    ax.set_xlabel(r"$J_1 - J_2\ [E_{r,s}]$")
    ax.set_ylabel(r"$\Delta\ [E_{r,s}]$")
    ax.axhline(0, color=AXISC, lw=0.8)
    ax.axvline(0, color=AXISC, lw=0.8)
    ax.set_title("Pump path winds once around the degeneracy\n"
                 r"$\Rightarrow$ Chern number $\nu_1=+1$ (band 1), $\nu_2=-1$ (band 2)")
    ax.legend(loc="lower right", fontsize=9.5)
    cb = fig.colorbar(lc, ax=ax, pad=0.02)
    cb.set_label(r"$\varphi/\pi$")
    ax.margins(0.12)
    return fig


def fig_gap_demo(d, ru, ra, Tpi, name, title_extra=""):
    Cc_rate = np.interp(PHIG, d["phis"], d["E12min"]) ** 2
    fig, ax = plt.subplots(3, 1, figsize=(8.6, 9.8), sharex=True,
                           constrained_layout=True)
    a0 = ax[0]
    a0.plot(ru["tr_phi"] / np.pi, ru["tr_x"], color=C_UNI, lw=1.0,
            label=rf"uniform $T={Tpi:g}\pi$ ($P_0$={ru['P0']:.3f}, $f$={ru['f']:+.3f})")
    a0.plot(ra["tr_phi"] / np.pi, ra["tr_x"], color=C_ADA, lw=1.8,
            label=rf"gap-adaptive $T={Tpi:g}\pi$ ($P_0$={ra['P0']:.3f}, $f$={ra['f']:+.3f})")
    a0.plot(SPH / np.pi, SPH / (2 * np.pi), color=MUTED, lw=1.2, ls=":",
            label=r"ideal $+d_l$/cycle")
    a0.axhline(1, color=AXISC, lw=0.8)
    a0.text(0.03, 1 + 0.03, r"$+d_l$", color=MUTED)
    a0.set_ylabel(r"$(\langle x\rangle-\langle x\rangle_0)/d_l$")
    a0.set_title(rf"(a) Lohse $V_s{{=}}10E_{{r,s}},V_l{{=}}5E_{{r,s}}$: "
                 rf"COM transport over one cycle{title_extra}")
    a0.legend(loc="upper left", fontsize=9.2)
    a1 = ax[1]
    a1.plot(SPH / np.pi, ru["p0"], "-o", color=C_UNI, ms=3.5,
            label=rf"uniform (final $P_0$={ru['P0']:.3f})")
    a1.plot(SPH / np.pi, ra["p0"], "-o", color=C_ADA, ms=3.5,
            label=rf"gap-adaptive (final $P_0$={ra['P0']:.3f})")
    a1.set_ylabel(r"band-0 population $P_0(\varphi)$")
    lo = min(0.55, min(np.nanmin(ru["p0"]), np.nanmin(ra["p0"])) - 0.05)
    a1.set_ylim(lo, 1.02)
    a1.set_title("(b) gap-adaptive stays in band 0 — no excitation, hence smooth")
    a1.legend(loc="lower left", fontsize=9.5)
    a2 = ax[2]
    a2.plot(PHIG / np.pi, Cc_rate / Cc_rate.max(), color=C_SCH, lw=1.8,
            label=r"$\dot\varphi\propto g(\varphi)^2$")
    a2.plot(d["phis"] / np.pi, d["E12min"] / d["E12min"].max(), color=MUTED,
            lw=1.2, ls="--", label=r"gap $g(\varphi)$ (norm.)")
    for xc in (0, 1, 2):
        a2.axvline(xc, color="#b08", lw=0.9, ls="--", alpha=0.35)
    a2.set_ylabel("normalized rate")
    a2.set_xlabel(r"pump phase $\varphi/\pi$")
    a2.set_xlim(0, 2)
    a2.set_title(r"(c) schedule crawls where the gap is smallest ($\varphi\approx 0,\pi,2\pi$)")
    a2.legend(loc="center right", fontsize=9.5)
    fig.suptitle(rf"Gap-adaptive pump — Lohse lattice, paper convention "
                 rf"($g_{{\min}}={d['g12_pump']:.3f}E_{{r,s}}$)", fontsize=12)
    fig.savefig(FIGDIR / name, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figs/{name}", flush=True)


def fig_scan(scan):
    fig, ax = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    Tu, Pu, fu = scan["T_uni"], scan["P0_uni"], scan["f_uni"]
    Ta, Pa, fa = scan["T_ada"], scan["P0_ada"], scan["f_ada"]
    ax.plot(Tu, Pu, "-o", color=C_UNI, lw=2.0, ms=5, label=r"uniform: $P_0(T)$")
    ax.plot(Tu, fu, "--o", color=C_UNI, lw=1.1, ms=3.5, alpha=0.6,
            label=r"uniform: transport $f(T)$")
    ax.plot(Ta, Pa, "-s", color=C_ADA, lw=2.0, ms=5, label=r"gap-adaptive: $P_0(T)$")
    ax.plot(Ta, fa, "--s", color=C_ADA, lw=1.1, ms=3.5, alpha=0.6,
            label=r"gap-adaptive: transport $f(T)$")
    ax.axhline(1.0, color=AXISC, lw=1.0, ls="--")
    ax.axhline(0.99, color=GRIDC, lw=0.9)
    ax.text(52, 0.992, r"$P_0=0.99$", color=MUTED, fontsize=8.5)
    Tbare = 2 * np.pi * scan["eta"] / np.pi
    ax.axvline(Tbare, color="#b08", lw=1.3, ls=":")
    ax.text(Tbare * 1.05, 0.18, rf"bare LZ scale $2\pi\eta\approx{Tbare:.0f}\pi$",
            color="#b08", rotation=90, va="bottom", fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel(r"pump cycle time $T\ [\pi\,\hbar/E_{r,s}]$")
    ax.set_ylabel("final-state metric")
    ax.set_ylim(0, 1.08)
    ax.set_title(r"How long to be adiabatic? uniform needs $T\gtrsim 2000\pi$; "
                 r"gap-adaptive $\sim 300\pi$")
    ax.legend(loc="lower right", fontsize=9.5, ncol=1)
    return fig


def save(fig, name):
    fig.savefig(FIGDIR / f"{name}.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figs/{name}.png", flush=True)


# ===========================================================================
# main
# ===========================================================================
def main():
    t_start = time.time()
    print("[1] continuum bands / gap / eta over the (k, phi) torus ...", flush=True)
    bd = compute_band_curves()
    print("[2] Wannier projection J1(phi), J2(phi), Delta(phi) ...", flush=True)
    wd = compute_wannier_curves()

    d = {**bd, **wd}
    d["g12_pump"] = float(bd["E12min"].min())
    d["g23_pump"] = float(bd["E23min"].min())
    d["eta"] = float(bd["eta_phi"].max())
    d["Delta_max"] = float(np.max(np.abs(wd["Delta"])))
    # phi = 0 exact band-edge combination
    E0c = np.linalg.eigvalsh(bloch_h_batch([0.0], 0.0))[0, :2]
    E0e = np.linalg.eigvalsh(bloch_h_batch([0.5], 0.0))[0, :2]
    g0, gpi = float(E0c[1] - E0c[0]), float(E0e[1] - E0e[0])
    d["g0"], d["gpi"] = g0, gpi
    d["J1_0"], d["J2_0"] = (g0 + gpi) / 4.0, (g0 - gpi) / 4.0
    d["g12_from_wannier"] = 2.0 * np.sqrt((wd["Delta"] / 2) ** 2
                                          + (wd["J1"] - wd["J2"]) ** 2)

    print("\n================ VALIDATION (old project values in parens) ============")
    print(f"  J1(0)      = {d['J1_0']:.5f} E_rs   (0.06074)")
    print(f"  J2(0)      = {d['J2_0']:.5f} E_rs   (0.00510)")
    print(f"  J1/J2 (0)  = {d['J1_0']/d['J2_0']:.2f}       (11.92)")
    print(f"  Delta_max  = {d['Delta_max']:.4f} E_rs  (4.5194)")
    print(f"  g12_pump   = {d['g12_pump']:.5f} E_rs  (0.11129)")
    print(f"  g23_pump   = {d['g23_pump']:.5f} E_rs  (1.04997)")
    print(f"  eta        = {d['eta']:.1f} /E_rs      (179.5)")
    print(f"  RM 2-band max residual = {float(d['rm_err'].max()):.4f} E_rs")
    print("========================================================================\n", flush=True)

    # schedule + CSV export (main.cpp two-column format)
    build_schedule(d["phis"], d["E12min"])
    s_grid = TCUM_UNIT
    with open(HERE / "gap_adaptive_lohse_schedule.csv", "w") as f:
        f.write("# gap-adaptive schedule for Lohse lattice Vs=10 E_rs, Vl=5 E_rs\n")
        f.write("# paper convention V = Vs cos^2(pi x/ds) + Vl sin^2(pi x/dl - phi/2)\n")
        f.write("# dphi/dt ~ g(phi)^2; s = t/T\n")
        f.write("s,phi\n")
        for si, pi_ in zip(s_grid, PHIG):
            f.write(f"{si:.12f},{pi_:.12f}\n")
    with open(HERE / "gap_adaptive_lohse_maincpp_schedule.csv", "w") as f:
        f.write("# same schedule mapped to the ECG main.cpp convention\n")
        f.write("# V^m = -Vs^m cos(4 pi x/a) - Vl^m cos(2 pi x/a + phi^m),\n")
        f.write("# Vs^m = Vs/2, Vl^m = Vl/2, phi^m(s) = 3 pi/2 - phi(s)  (pumps -a/cycle)\n")
        f.write("s,phi\n")
        for si, pi_ in zip(s_grid, PHIG):
            f.write(f"{si:.12f},{3*np.pi/2 - pi_:.12f}\n")
    print("[write] gap_adaptive_lohse_schedule.csv (+ maincpp variant)", flush=True)

    print("[3] pump dynamics: demo runs at T = 400 pi ...", flush=True)
    init_dynamics()
    Tpi_demo = 400.0
    ra = run_pump(Tpi_demo * np.pi, "adaptive")
    print(f"    adaptive T=400pi: f={ra['f']:+.4f}  P0={ra['P0']:.4f}", flush=True)
    ru = run_pump(Tpi_demo * np.pi, "uniform")
    print(f"    uniform  T=400pi: f={ru['f']:+.4f}  P0={ru['P0']:.4f}", flush=True)
    fig_gap_demo(d, ru, ra, Tpi_demo, "fig5_gap_adaptive_pump.png")

    print("[4] uniform pump at T = 2000 pi (adiabatic check) ...", flush=True)
    ru2 = run_pump(2000.0 * np.pi, "uniform", dt=0.025)
    ra2 = run_pump(2000.0 * np.pi, "adaptive", dt=0.025)
    print(f"    uniform  T=2000pi: f={ru2['f']:+.4f}  P0={ru2['P0']:.4f}", flush=True)
    print(f"    adaptive T=2000pi: f={ra2['f']:+.4f}  P0={ra2['P0']:.4f}", flush=True)
    fig_gap_demo(d, ru2, ra2, 2000.0, "fig7_uniform_T2000pi.png",
                 title_extra=" (uniform now nearly adiabatic)")

    print("[5] adiabatic time scan (uniform + gap-adaptive) ...", flush=True)
    T_uni = np.array([100, 200, 400, 600, 800, 1000, 1400, 1800, 2200, 2600, 3000, 3600])
    T_ada = np.array([50, 100, 150, 200, 300, 400, 600, 800])
    scan = dict(T_uni=T_uni, T_ada=T_ada, eta=d["eta"])
    for key, Ts, mode in (("uni", T_uni, "uniform"), ("ada", T_ada, "adaptive")):
        fs, ps = [], []
        for Tpi in Ts:
            r = run_pump(Tpi * np.pi, mode, dt=0.025, traces=False)
            fs.append(r["f"])
            ps.append(r["P0"])
            print(f"    {mode:8s} T={Tpi:5.0f}pi  f={r['f']:+.4f}  P0={r['P0']:.4f}"
                  f"  leak={100*(1-r['P0']):.2f}%", flush=True)
        scan[f"f_{key}"] = np.array(fs)
        scan[f"P0_{key}"] = np.array(ps)

    def threshold_T(Ts, Ps, thr):
        ok = Ps >= thr
        for i in range(len(Ts)):
            if ok[i:].all():
                return float(Ts[i])
        return float("nan")

    T99_uni = threshold_T(T_uni, scan["P0_uni"], 0.99)
    T99_ada = threshold_T(T_ada, scan["P0_ada"], 0.99)
    T999_ada = threshold_T(T_ada, scan["P0_ada"], 0.999)

    print("[6] figures ...", flush=True)
    save(fig_potential(), "fig1_potential")
    save(fig_bands(), "fig2_bands")
    save(fig_rm_params(d), "fig3_rice_mele_params")
    save(fig_pump_path(d), "fig4_pump_path")
    save(fig_scan(scan), "fig6_adiabatic_time_scan")

    np.savez(HERE / "lohse_reference_data.npz", **d, **scan,
             demo_Tpi=Tpi_demo,
             demo_uni_phi=ru["tr_phi"], demo_uni_x=ru["tr_x"], demo_uni_p0=ru["p0"],
             demo_ada_phi=ra["tr_phi"], demo_ada_x=ra["tr_x"], demo_ada_p0=ra["p0"],
             SPH=SPH)

    ms = 1e3 * np.pi * T_UNIT_S     # one "pi hbar/E_rs" in ms
    summary = dict(
        VS_Ers=VS, VL_Ers=VL, VL_Erl=4 * VL,
        E_rs_Hz=E_RS_HZ, hbar_over_Ers_us=1e6 * T_UNIT_S,
        T_code_factor=T_CODE_FACTOR,
        J1_0=d["J1_0"], J2_0=d["J2_0"], J1J2_ratio_0=d["J1_0"] / d["J2_0"],
        Delta_max=d["Delta_max"], g12_pump=d["g12_pump"], g23_pump=d["g23_pump"],
        g0=d["g0"], gpi=d["gpi"], eta=d["eta"], rm_err_max=float(d["rm_err"].max()),
        J1_0_Hz=d["J1_0"] * E_RS_HZ, Delta_max_Hz=d["Delta_max"] * E_RS_HZ,
        g12_pump_Hz=d["g12_pump"] * E_RS_HZ,
        demo=dict(Tpi=Tpi_demo,
                  uniform=dict(f=ru["f"], P0=ru["P0"]),
                  adaptive=dict(f=ra["f"], P0=ra["P0"]),
                  Tpi2000_uniform=dict(f=ru2["f"], P0=ru2["P0"]),
                  Tpi2000_adaptive=dict(f=ra2["f"], P0=ra2["P0"])),
        scan=dict(T_uni=T_uni.tolist(), P0_uni=scan["P0_uni"].tolist(),
                  f_uni=scan["f_uni"].tolist(),
                  T_ada=T_ada.tolist(), P0_ada=scan["P0_ada"].tolist(),
                  f_ada=scan["f_ada"].tolist()),
        T99_uni_pi=T99_uni, T99_ada_pi=T99_ada, T999_ada_pi=T999_ada,
        T_bare_LZ_pi=2 * d["eta"],
        ms_per_pi_unit=ms,
        T99_uni_ms=T99_uni * ms, T99_ada_ms=T99_ada * ms,
    )
    with open(HERE / "lohse_reference_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[write] lohse_reference_summary.json")
    print(f"\nkey answer:  uniform needs T >~ {T99_uni:.0f} pi hbar/E_rs "
          f"({T99_uni*ms:.0f} ms for 87Rb) for P0 >= 0.99;")
    print(f"             gap-adaptive reaches P0 >= 0.99 at T ~ {T99_ada:.0f} pi "
          f"({T99_ada*ms:.0f} ms), >= 0.999 at T ~ {T999_ada:.0f} pi.")
    print(f"total wall time: {time.time()-t_start:.0f} s")


if __name__ == "__main__":
    main()
