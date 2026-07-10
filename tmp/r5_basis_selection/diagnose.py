#!/usr/bin/env python3
"""Diagnose the free single-particle representability: for each snapshot fit
phiL(t) (and phiR) with n free complex gaussians (center, width, phase) and
report infidelity. This localizes the K32 span gap: centers vs widths vs
genuine non-gaussian spreading."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "rice_mele_reference" / "lohes_experience"))
import n2_grid_reference as gr

A_LAT = 8.0; CELLS = 16; LENGTH = CELLS * A_LAT; NG = 2048
DX = LENGTH / NG; X = np.arange(NG) * DX - LENGTH / 2
P = 2 * np.pi * np.fft.fftfreq(NG, d=DX); P2 = P * P
ER = gr.recoil_energy(A_LAT); VS = 3.0 * ER; VL = 3.0 * ER
DT = 0.01; PERIOD = 160.0 * np.pi
SCHEDULE = REPO / "rice_mele_reference" / "Vs3Vl3_3_3" / "gap_adaptive_vs3vl3_full_depth_schedule.csv"
N1CSV = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N1_K16.csv"
PHI_PI = np.round(np.arange(0.30, 0.6001, 0.02), 4)


def load_n1():
    rows = [l.strip() for l in N1CSV.read_text().splitlines() if l.strip()]
    return [(complex(*map(float, rows[5*k].split(","))).real,
             complex(*map(float, rows[5*k+2].split(","))).real,
             complex(*map(float, rows[5*k+3].split(","))).real)
            for k in range(len(rows)//5)]


def gauss(b, r):
    d = (X - r + LENGTH/2) % LENGTH - LENGTH/2
    return np.exp(-b*d*d)


def evolve():
    ss, phis = gr.load_schedule(SCHEDULE)
    n1 = load_n1()
    pL = np.zeros(NG, complex); pR = np.zeros(NG, complex)
    for u, b, r in n1:
        pL += u*gauss(b, r); pR += u*gauss(b, r+A_LAT)
    kin = np.exp(-1j*DT*0.5*P2)
    def step(psi, t):
        vm = gr.v_lattice(X, gr.phi_at(t+0.5*DT, PERIOD, ss, phis), A_LAT, VS, VL)
        psi = psi*np.exp(-0.5j*DT*vm); psi = np.fft.ifft(kin*np.fft.fft(psi))
        return psi*np.exp(-0.5j*DT*vm)
    snaps = []; ti = 0; t = 0.0
    targets = list(PHI_PI)
    phin = gr.phi_at(t, PERIOD, ss, phis)
    while ti < len(targets):
        if phin/np.pi >= targets[ti]-1e-12:
            snaps.append((targets[ti], pL.copy(), pR.copy())); ti += 1; continue
        pL = step(pL, t); pR = step(pR, t); t += DT
        phin = gr.phi_at(t, PERIOD, ss, phis)
    return snaps


def ov(f, h):
    return np.sum(np.conj(f)*h)*DX


def fit_n_gauss(psi, n, starts):
    """Fit psi ~ sum_k (ak+i bk) exp(-wk (x-ck)^2); return infidelity."""
    psi = psi/np.sqrt(np.real(ov(psi, psi)))
    def resid(params):
        model = np.zeros(NG, complex)
        for k in range(n):
            a, bb, w, c = params[4*k:4*k+4]
            w = abs(w)+1e-3
            model += (a+1j*bb)*gauss(w, c)
        r = (model-psi)
        return np.concatenate([r.real, r.imag])
    best = None
    for st in starts:
        try:
            sol = least_squares(resid, st, method="lm", max_nfev=4000)
            model = np.zeros(NG, complex)
            for k in range(n):
                a, bb, w, c = sol.x[4*k:4*k+4]; w = abs(w)+1e-3
                model += (a+1j*bb)*gauss(w, c)
            fid = np.abs(ov(model, psi))**2/np.real(ov(model, model))
            infid = 1-fid
            if best is None or infid < best[0]:
                best = (infid, sol.x)
        except Exception:
            pass
    return best[0] if best else 1.0


def main():
    snaps = evolve()
    print("free single-particle representability (phiL packet):")
    print(f"{'phi/pi':>7} {'<x>':>7} {'rms':>6} {'infid_1g':>9} {'infid_2g':>9} {'infid_3g':>9}")
    for phipi, pL, pR in snaps:
        psi = pL/np.sqrt(np.real(ov(pL, pL)))
        prob = np.abs(psi)**2
        xm = np.sum(X*prob)*DX
        rms = np.sqrt(np.sum((X-xm)**2*prob)*DX)
        # starting guesses
        s1 = [[1.0, 0.0, 0.57, xm]]
        s2 = [[0.7, 0, 0.57, xm-2], [0.7, 0, 0.57, xm+2],
              [0.7, 0, 0.57, 0.0], [0.7, 0, 0.57, -4.0]]
        s2 = [np.array(s2[0]+s2[1]), np.array(s2[2]+s2[3])]
        s3 = [np.array([0.6,0,0.57,xm-4, 0.6,0,0.57,xm, 0.6,0,0.57,xm+4]),
              np.array([0.6,0,0.57,0.0, 0.6,0,0.57,-4.0, 0.6,0,0.57,4.0])]
        i1 = fit_n_gauss(pL, 1, s1)
        i2 = fit_n_gauss(pL, 2, s2)
        i3 = fit_n_gauss(pL, 3, s3)
        print(f"{phipi:7.2f} {xm:7.3f} {rms:6.3f} {i1:9.2e} {i2:9.2e} {i3:9.2e}")


if __name__ == "__main__":
    main()
