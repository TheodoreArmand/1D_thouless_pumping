#!/usr/bin/env python3
"""Rule out numerical artifacts behind the plateau infidelity floor:
 (a) dt convergence of the free evolution (dt=0.01 vs 0.002): norm drift and
     the single-particle infidelity of a broad-gaussian fit at the plateau;
 (b) single-particle free-gaussian matching-pursuit floor at phi=0.60pi
     (N=1..6 free gaussians) -- the intrinsic representability bound."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "rice_mele_reference" / "lohes_experience"))
import n2_grid_reference as gr

A_LAT = 8.0; CELLS = 16; LENGTH = CELLS*A_LAT; NG = 2048
DX = LENGTH/NG; X = np.arange(NG)*DX - LENGTH/2
P = 2*np.pi*np.fft.fftfreq(NG, d=DX); P2 = P*P
ER = gr.recoil_energy(A_LAT); VS = 3.0*ER; VL = 3.0*ER
PERIOD = 160.0*np.pi
SCHEDULE = REPO/"rice_mele_reference"/"Vs3Vl3_3_3"/"gap_adaptive_vs3vl3_full_depth_schedule.csv"
N1CSV = REPO/"initial_state"/"Vs3Er_Vl3Er"/"initial_pathpad_N1_K16.csv"


def load_n1():
    rows = [l.strip() for l in N1CSV.read_text().splitlines() if l.strip()]
    return [(complex(*map(float, rows[5*k].split(","))).real,
             complex(*map(float, rows[5*k+2].split(","))).real,
             complex(*map(float, rows[5*k+3].split(","))).real)
            for k in range(len(rows)//5)]

def gauss(b, r):
    return gr.gaussian(X, LENGTH, b, r)

def ov(f, h):
    return np.sum(np.conj(f)*h)*DX


def evolve_to_phi(dt, phi_pi_target):
    ss, phis = gr.load_schedule(SCHEDULE)
    n1 = load_n1()
    pL = np.zeros(NG, complex)
    for u, b, r in n1:
        pL += u*gauss(b, r)
    kin = np.exp(-1j*dt*0.5*P2); t = 0.0
    phin = gr.phi_at(t, PERIOD, ss, phis)
    norm0 = np.real(ov(pL, pL))
    while phin/np.pi < phi_pi_target - 1e-12:
        vm = gr.v_lattice(X, gr.phi_at(t+0.5*dt, PERIOD, ss, phis), A_LAT, VS, VL)
        pL = pL*np.exp(-0.5j*dt*vm); pL = np.fft.ifft(kin*np.fft.fft(pL))
        pL = pL*np.exp(-0.5j*dt*vm); t += dt
        phin = gr.phi_at(t, PERIOD, ss, phis)
    drift = np.real(ov(pL, pL))/norm0 - 1.0
    return pL, drift


def fit_free(psi, n, ntries=8):
    psi = psi/np.sqrt(np.real(ov(psi, psi)))
    prob = np.abs(psi)**2; xm = np.sum(X*prob)*DX
    rms = np.sqrt(np.sum((X-xm)**2*prob)*DX)
    def resid(pp):
        m = np.zeros(NG, complex)
        for k in range(n):
            a, bb, w, c = pp[4*k:4*k+4]; w = abs(w)+1e-3
            m += (a+1j*bb)*gauss(w, c)
        return np.concatenate([(m-psi).real, (m-psi).imag])
    rng = np.random.default_rng(0); best = 1.0
    for _ in range(ntries):
        p0 = []
        for k in range(n):
            p0 += [rng.uniform(-1, 1), rng.uniform(-1, 1),
                   float(1.0/(2*max(rms, 0.3)**2))*rng.uniform(0.4, 2.0),
                   xm + rng.uniform(-3, 3)]
        try:
            sol = least_squares(resid, np.array(p0), method="lm", max_nfev=6000)
            m = np.zeros(NG, complex)
            for k in range(n):
                a, bb, w, c = sol.x[4*k:4*k+4]; w = abs(w)+1e-3
                m += (a+1j*bb)*gauss(w, c)
            infid = 1 - np.abs(ov(m, psi))**2/np.real(ov(m, m))
            best = min(best, infid)
        except Exception:
            pass
    return best


def main():
    print("(a) dt convergence of free evolution to phi=0.60pi:")
    pL_coarse, d_c = evolve_to_phi(0.01, 0.60)
    pL_fine, d_f = evolve_to_phi(0.002, 0.60)
    a = pL_coarse/np.sqrt(np.real(ov(pL_coarse, pL_coarse)))
    b = pL_fine/np.sqrt(np.real(ov(pL_fine, pL_fine)))
    mism = 1 - np.abs(ov(a, b))**2
    print(f"    norm drift: dt0.01={d_c:.2e}  dt0.002={d_f:.2e}")
    print(f"    1-|<psi_dt0.01|psi_dt0.002>|^2 = {mism:.2e}  (numerical, not physical)")

    print("(b) single-particle free-gaussian matching-pursuit floor at phi=0.60pi:")
    for n in (1, 2, 3, 4, 5, 6):
        infid = fit_free(pL_fine, n)
        print(f"    N={n} free gaussians -> single-particle infid = {infid:.3e}")
    print("    (2-particle plateau floor ~ this^(product); needs N^2/2 sym terms)")


if __name__ == "__main__":
    main()
