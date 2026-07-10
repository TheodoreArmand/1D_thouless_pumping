#!/usr/bin/env python3
"""
R5 union basis selection (Gate-1, first transition) for the Vs3/Vl3 N=2 K32 pump.

Self-contained, lives in tmp/, touches no production code. Implements the R5
proposal (txt_task/r5_basis_union_proposal.md) end-to-end:

  1. Generate the FREE grid reference wavefunction snapshots psi_grid(x1,x2) at
     phi/pi in {0.30, 0.32, ..., 0.60} using the NEW full-depth schedule.
     For g=0 the split-step is exactly separable, so
        psi(x1,x2,t) = ( phiL(x1,t) phiR(x2,t) + phiR(x1,t) phiL(x2,t) ) / sqrt2
     with phiL, phiR two independently-evolved 1D packets. Every overlap the
     selection needs then factorizes into exact 1D overlaps -- no 2D arrays.
  2. Build the structured candidate pool P-A (32 forced) + P-B + P-C1 + P-C3.
  3. Greedy, conditioning-guarded selection up to K_target (novelty veto).
  4. Report: per-snapshot infidelity K32 vs selection, S-spectrum, wake order.
  5. Refit initial-state u over the selected set (+ u-floor), verify fidelity,
     write the union basis CSV into tmp/.

Reference imports come from the canonical grid script so the physics
(potential, schedule interpolation, recoil energy, ring gaussians) is identical.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "rice_mele_reference" / "lohes_experience"))
import n2_grid_reference as gr  # noqa: E402  (canonical primitives)

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
A_LAT = 8.0
CELLS = 16
LENGTH = CELLS * A_LAT               # 128
NG = 2048                            # 1D grid; free evolution is cheap
DX = LENGTH / NG
X = np.arange(NG) * DX - LENGTH / 2
P = 2 * np.pi * np.fft.fftfreq(NG, d=DX)
P2 = P * P

ER = gr.recoil_energy(A_LAT)         # 2 pi^2 / a^2 = pi^2/32
VS = 3.0 * ER
VL = 3.0 * ER
DT = 0.01
PERIOD = 160.0 * np.pi               # T = 160 pi

SCHEDULE = REPO / "rice_mele_reference" / "Vs3Vl3_3_3" / "gap_adaptive_vs3vl3_full_depth_schedule.csv"
N1CSV = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N1_K16.csv"
K32CSV = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N2_K32.csv"

PHI_PI_TARGETS = np.round(np.arange(0.30, 0.6001, 0.02), 4)   # 16 snapshots
K_TARGET = 48
INFID_STOP = 1e-3          # early-stop when worst snapshot <= this
NOVELTY_VETO = 0.005       # skip candidate if 1-||P_I ghat_j||^2 < this
REG = 1e-10                # Tikhonov on normalized Gram solves
U_FLOOR_FRAC = 1e-4        # reserve u floor = frac * max|u|

BROAD = 0.57436677743869635
NARROW = 2.1961082666773684

# ----------------------------------------------------------------------------
# 1D helpers (real gaussians + exact overlaps)
# ----------------------------------------------------------------------------

def gauss1d(b, r):
    return gr.gaussian(X, LENGTH, b, r)

def ov(f, h):
    """<f|h> = sum conj(f) h dx on the 1D grid."""
    return np.sum(np.conj(f) * h) * DX


def load_n1():
    rows = [l.strip() for l in N1CSV.read_text().splitlines() if l.strip()]
    out = []
    for k in range(len(rows) // 5):
        u = complex(*map(float, rows[5 * k].split(","))).real
        b = complex(*map(float, rows[5 * k + 2].split(","))).real
        r = complex(*map(float, rows[5 * k + 3].split(","))).real
        out.append((u, b, r))
    return out


def load_k32():
    rows = [l.strip() for l in K32CSV.read_text().splitlines() if l.strip()]
    stride = 12
    out = []
    for k in range(len(rows) // stride):
        o = stride * k
        u = complex(*map(float, rows[o].split(","))).real
        b0 = complex(*map(float, rows[o + 5].split(","))).real
        b1 = complex(*map(float, rows[o + 8].split(","))).real
        r0 = complex(*map(float, rows[o + 9].split(","))).real
        r1 = complex(*map(float, rows[o + 10].split(","))).real
        out.append((u, b0, r0, b1, r1))
    return out


# ----------------------------------------------------------------------------
# Free single-particle evolution -> snapshots by phi
# ----------------------------------------------------------------------------

def evolve_snapshots():
    ss, phis = gr.load_schedule(SCHEDULE)
    n1 = load_n1()

    # exact initial adjacent-cell packets (the true symmetrized-product target)
    phiL = np.zeros(NG, dtype=np.complex128)
    phiR = np.zeros(NG, dtype=np.complex128)
    for u, b, r in n1:
        phiL += u * gauss1d(b, r)
        phiR += u * gauss1d(b, r + A_LAT)
    phiL = phiL.astype(np.complex128)
    phiR = phiR.astype(np.complex128)

    kin = np.exp(-1j * DT * 0.5 * P2)

    def step(psi, t):
        phi_mid = gr.phi_at(t + 0.5 * DT, PERIOD, ss, phis)
        v = gr.v_lattice(X, phi_mid, A_LAT, VS, VL)
        psi = psi * np.exp(-0.5j * DT * v)
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = psi * np.exp(-0.5j * DT * v)
        return psi

    targets = list(PHI_PI_TARGETS)
    snaps = []  # (phi_pi, phiL_t, phiR_t)
    ti = 0
    t = 0.0
    phi_now = gr.phi_at(t, PERIOD, ss, phis)
    while ti < len(targets):
        if phi_now / np.pi >= targets[ti] - 1e-12:
            snaps.append((targets[ti], phiL.copy(), phiR.copy(), phi_now))
            ti += 1
            continue
        phiL = step(phiL, t)
        phiR = step(phiR, t)
        t += DT
        phi_now = gr.phi_at(t, PERIOD, ss, phis)
        if t > 1.05 * PERIOD:
            break
    return snaps


# ----------------------------------------------------------------------------
# Candidate pool
# ----------------------------------------------------------------------------

def canon_key(b0, r0, b1, r1):
    a = (round(b0, 6), round(r0, 6))
    b = (round(b1, 6), round(r1, 6))
    return tuple(sorted([a, b]))


def build_pool():
    """Return list of dicts: {b0,r0,b1,r1, tier, forced, u0(optional)}."""
    pool = []
    seen = {}

    def add(b0, r0, b1, r1, tier, forced=False, u0=0.0):
        key = canon_key(b0, r0, b1, r1)
        if key in seen:
            return
        seen[key] = len(pool)
        pool.append(dict(b0=b0, r0=r0, b1=b1, r1=r1, tier=tier,
                         forced=forced, u0=u0, key=key))

    # P-A: the 32 existing terms (forced, never dropped)
    for (u, b0, r0, b1, r1) in load_k32():
        add(b0, r0, b1, r1, "P-A", forced=True, u0=u)

    widths = [(BROAD, BROAD), (BROAD, NARROW), (NARROW, BROAD), (NARROW, NARROW)]

    # P-B: COM-diagonal double-move centers
    for shift in (-2.0, -4.0, -6.0, -8.0):
        for w0, w1 in widths:
            add(w0, shift, w1, shift + A_LAT, "P-B")

    # P-C1: first-transition staggered (one to new well, one at barrier)
    for c0, c1 in ((-4.0, 6.0), (-2.0, 4.0)):
        for w0, w1 in widths:
            add(w0, c0, w1, c1, "P-C1")

    # P-C3: missing joint-width combos of existing one-move endpoint centers
    c3_centers = [(-2.0, 8.0), (-4.0, 8.0), (-6.0, 8.0), (-8.0, 8.0),
                  (0.0, 6.0), (0.0, 4.0), (0.0, 2.0), (0.0, 0.0)]
    for c0, c1 in c3_centers:
        for w0, w1 in widths:
            add(w0, c0, w1, c1, "P-C3")

    return pool


# ----------------------------------------------------------------------------
# Overlap tables (all exact via 1D factorization)
# ----------------------------------------------------------------------------

def build_overlap_tables(pool, snaps):
    # distinct 1D gaussians across the pool
    distinct = {}
    def gidx(b, r):
        k = (round(b, 6), round(r, 6))
        if k not in distinct:
            distinct[k] = (len(distinct), gauss1d(b, r))
        return distinct[k][0]

    for t in pool:
        t["i0"] = gidx(t["b0"], t["r0"])
        t["i1"] = gidx(t["b1"], t["r1"])

    nd = len(distinct)
    G = np.zeros((nd, NG))
    for (idx, arr) in distinct.values():
        G[idx] = arr.real
    ov11 = (G @ G.T) * DX                      # nd x nd real

    nsnap = len(snaps)
    ovL = np.zeros((nd, nsnap), dtype=np.complex128)
    ovR = np.zeros((nd, nsnap), dtype=np.complex128)
    psi_norm = np.zeros(nsnap)
    for m, (_pi, pL, pR, _phi) in enumerate(snaps):
        ovL[:, m] = (G @ pL) * DX
        ovR[:, m] = (G @ pR) * DX
        nLL = np.real(ov(pL, pL)); nRR = np.real(ov(pR, pR)); nLR = ov(pL, pR)
        psi_norm[m] = np.sqrt(nLL * nRR + np.abs(nLR) ** 2)

    ncand = len(pool)
    i0 = np.array([t["i0"] for t in pool])
    i1 = np.array([t["i1"] for t in pool])
    # raw norms of symmetrized primitives S[g0,g1]
    cand_norm = np.sqrt(ov11[i0, i0] * ov11[i1, i1] + ov11[i0, i1] ** 2)

    # normalized Gram Sfull[c,c']
    S = (ov11[np.ix_(i0, i0)] * ov11[np.ix_(i1, i1)]
         + ov11[np.ix_(i0, i1)] * ov11[np.ix_(i1, i0)])
    S = S / np.outer(cand_norm, cand_norm)
    np.fill_diagonal(S, 1.0)

    # normalized overlaps with each snapshot B[c,m] = <ghat_c|psihat_m>
    B = (ovL[i0] * ovR[i1] + ovR[i0] * ovL[i1])
    B = B / cand_norm[:, None] / psi_norm[None, :]

    return S, B, cand_norm, psi_norm


# ----------------------------------------------------------------------------
# Projection infidelity / greedy selection
# ----------------------------------------------------------------------------

def infid_for_set(S, B, I):
    """1 - ||P_I psi_m||^2 for each snapshot m (regularized normal equations)."""
    SII = S[np.ix_(I, I)] + REG * np.eye(len(I))
    bI = B[np.ix_(I, range(B.shape[1]))]          # |I| x nsnap
    c = np.linalg.solve(SII, bI)                   # complex
    proj = np.real(np.sum(np.conj(bI) * c, axis=0))
    return np.clip(1.0 - proj, 0.0, 1.0), c


def greedy_select(S, B, pool, k_target):
    forced = [i for i, t in enumerate(pool) if t["forced"]]
    I = list(forced)
    nsnap = B.shape[1]
    history = []
    infid, c = infid_for_set(S, B, I)
    history.append(("start", None, float(infid.max()), None, None))

    while len(I) < k_target and infid.max() > INFID_STOP:
        SII = S[np.ix_(I, I)] + REG * np.eye(len(I))
        SII_inv = np.linalg.inv(SII)
        c = SII_inv @ B[np.ix_(I, range(nsnap))]   # |I| x nsnap

        best_j, best_score, best_nov = -1, -1.0, None
        for j in range(len(pool)):
            if j in I:
                continue
            sjI = S[j, I]
            nov = 1.0 - float(sjI @ SII_inv @ sjI)   # 1 - ||P_I ghat_j||^2
            if nov < NOVELTY_VETO:
                continue
            # residual overlap per snapshot: b_jm - S_jI c_Im
            resov = B[j, :] - sjI @ c
            score = float(np.max(np.abs(resov) ** 2) / (nov + 1e-6))
            if score > best_score:
                best_score, best_j, best_nov = score, j, nov
        if best_j < 0:
            history.append(("veto-exhausted", None, float(infid.max()), None, None))
            break
        I.append(best_j)
        infid, _ = infid_for_set(S, B, I)
        # which snapshot did it help most
        args = int(np.argmax(np.abs(B[best_j, :])))
        history.append((pool[best_j]["tier"], best_j, float(infid.max()),
                        best_nov, args))
    return I, history, infid


# ----------------------------------------------------------------------------
# Initial-state refit over selected terms (+ u floor)
# ----------------------------------------------------------------------------

def refit_initial(pool, I):
    """Least-squares fit the symmetrized adjacent-cell product target with the
    selected terms; return u vector, initial fidelity, floor-adjusted u."""
    n1 = load_n1()
    phiL = np.zeros(NG); phiR = np.zeros(NG)
    for u, b, r in n1:
        phiL += u * gauss1d(b, r).real
        phiR += u * gauss1d(b, r + A_LAT).real
    phiL /= np.sqrt(np.sum(phiL * phiL) * DX)
    phiR /= np.sqrt(np.sum(phiR * phiR) * DX)
    # target normalization constant
    nLL = np.sum(phiL * phiL) * DX; nRR = np.sum(phiR * phiR) * DX
    nLR = np.sum(phiL * phiR) * DX
    tnorm = np.sqrt(nLL * nRR + nLR ** 2)

    # design overlaps via 1D factorization
    sel = [pool[i] for i in I]
    # 1D overlaps of each selected primitive gaussian with phiL, phiR
    def g(b, r):
        return gauss1d(b, r).real
    gL = [ov(g(t["b0"], t["r0"]), phiL).real for t in sel]
    gR0 = [ov(g(t["b0"], t["r0"]), phiR).real for t in sel]
    gL1 = [ov(g(t["b1"], t["r1"]), phiL).real for t in sel]
    gR1 = [ov(g(t["b1"], t["r1"]), phiR).real for t in sel]
    gL = np.array(gL); gR0 = np.array(gR0); gL1 = np.array(gL1); gR1 = np.array(gR1)
    # <S[g0,g1] | S[phiL,phiR]> = <g0|phiL><g1|phiR> + <g0|phiR><g1|phiL>
    b_vec = (gL * gR1 + gR0 * gL1) / tnorm

    # Gram of selected raw primitives
    ii0 = [t["b0"] for t in sel]
    # rebuild via 1D overlaps between selected gaussians
    G0 = np.array([g(t["b0"], t["r0"]) for t in sel])
    G1 = np.array([g(t["b1"], t["r1"]) for t in sel])
    o00 = (G0 @ G0.T) * DX; o11 = (G1 @ G1.T) * DX
    o01 = (G0 @ G1.T) * DX; o10 = o01.T
    Gram = o00 * o11 + o01 * o10
    u = np.linalg.solve(Gram + 1e-12 * np.eye(len(sel)), b_vec)
    # fidelity of fit
    fit_norm2 = float(u @ Gram @ u)
    overlap = float(u @ b_vec)
    fidelity = overlap ** 2 / fit_norm2   # target normalized
    # normalize u to unit state
    u = u / np.sqrt(fit_norm2)
    # apply floor to parked reserves
    umax = np.max(np.abs(u))
    floor = U_FLOOR_FRAC * umax
    u_floored = np.where(np.abs(u) < floor, floor, u)
    return u, u_floored, fidelity, floor


def write_union_csv(path, pool, I, u):
    sel = [pool[i] for i in I]
    with open(path, "w") as f:
        for k, (t, uk) in enumerate(zip(sel, u)):
            f.write(f"{uk:.17g},0\n")
            for _ in range(4):
                f.write("0,0\n")            # A = 0
            f.write(f"{t['b0']:.17g},0\n")
            f.write("0,0\n"); f.write("0,0\n")
            f.write(f"{t['b1']:.17g},0\n")
            f.write(f"{t['r0']:.17g},0\n")
            f.write(f"{t['r1']:.17g},0\n")
            f.write(f"{k}\n")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    print("[1] evolving free single-particle packets (full-depth schedule)...", flush=True)
    snaps = evolve_snapshots()
    print(f"    captured {len(snaps)} snapshots at phi/pi = "
          + ", ".join(f"{s[0]:.2f}" for s in snaps))
    for s in snaps[:3] + snaps[-1:]:
        print(f"      target={s[0]:.2f}  actual phi/pi={s[3]/np.pi:.4f}")

    print("[2] building candidate pool ...", flush=True)
    pool = build_pool()
    ntier = {}
    for t in pool:
        ntier[t["tier"]] = ntier.get(t["tier"], 0) + 1
    print(f"    pool size = {len(pool)}  by tier: {ntier}")

    print("[3] overlap tables (exact 1D factorization) ...", flush=True)
    S, B, cand_norm, psi_norm = build_overlap_tables(pool, snaps)

    # K32 baseline infidelity
    forced = [i for i, t in enumerate(pool) if t["forced"]]
    infid32, _ = infid_for_set(S, B, forced)
    print("[4] K32 baseline infidelity by phi/pi:")
    for m, s in enumerate(snaps):
        print(f"      phi/pi={s[0]:.2f}  infid_K32={infid32[m]:.3e}")
    print(f"    K32 worst infidelity = {infid32.max():.3e} "
          f"at phi/pi={snaps[int(np.argmax(infid32))][0]:.2f}")

    print(f"[5] greedy selection up to K={K_TARGET} "
          f"(stop at infid<={INFID_STOP}, novelty veto<{NOVELTY_VETO}) ...")
    I, history, infidK = greedy_select(S, B, pool, K_TARGET)
    print(f"    selected K = {len(I)} terms ({len(I)-32} additions)")
    print("    selection order:")
    for step, (tier, j, worst, nov, argsnap) in enumerate(history):
        if j is None:
            print(f"      round {step:2d}: {tier:14s} worst_infid={worst:.3e}")
        else:
            t = pool[j]
            sp = snaps[argsnap][0]
            print(f"      round {step:2d}: +{tier:5s} "
                  f"p0=(b{t['b0']:.3g},r{t['r0']:+.1f}) p1=(b{t['b1']:.3g},r{t['r1']:+.1f}) "
                  f"novelty={nov:.3f} peakphi/pi={sp:.2f} -> worst_infid={worst:.3e}")

    # final infidelity comparison
    infidF, _ = infid_for_set(S, B, I)
    print("[6] infidelity K32 vs K%d:" % len(I))
    for m, s in enumerate(snaps):
        print(f"      phi/pi={s[0]:.2f}  K32={infid32[m]:.3e}  K{len(I)}={infidF[m]:.3e}")
    print(f"    worst: K32={infid32.max():.3e}  K{len(I)}={infidF.max():.3e}")

    # conditioning of selected normalized Gram
    Ssel = S[np.ix_(I, I)]
    eig = np.linalg.eigvalsh(Ssel)
    print(f"[7] selected normalized Gram spectrum: "
          f"lambda_min={eig[0]:.3e} lambda_max={eig[-1]:.3e} "
          f"cond={eig[-1]/max(eig[0],1e-300):.3e}")

    # initial-state refit + fidelity
    u, u_floored, fid0, floor = refit_initial(pool, I)
    print(f"[8] initial-state refit fidelity |<psi_K{len(I)}(0)|target>|^2 = {fid0:.10f}")
    print(f"    u floor = {floor:.3e}  (reserves parked below u_on ~ 0.02-0.05)")

    outcsv = HERE / f"initial_pathpad_N2_K{len(I)}_union.csv"
    write_union_csv(outcsv, pool, I, u_floored)
    print(f"[9] wrote union basis CSV -> {outcsv}")

    # dump machine-readable summary
    summary = dict(
        n_snapshots=len(snaps),
        phi_pi=[float(s[0]) for s in snaps],
        pool_size=len(pool),
        pool_by_tier=ntier,
        K_selected=len(I),
        infid_K32=[float(v) for v in infid32],
        infid_Ksel=[float(v) for v in infidF],
        worst_K32=float(infid32.max()),
        worst_Ksel=float(infidF.max()),
        gram_lambda_min=float(eig[0]),
        gram_lambda_max=float(eig[-1]),
        init_fidelity=float(fid0),
        u_floor=float(floor),
        selected=[dict(tier=pool[i]["tier"], b0=pool[i]["b0"], r0=pool[i]["r0"],
                       b1=pool[i]["b1"], r1=pool[i]["r1"], forced=pool[i]["forced"])
                  for i in I],
    )
    (HERE / "r5_selection_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[10] wrote summary -> {HERE/'r5_selection_summary.json'}")

    # save snapshots (reduced) for optional plotting
    np.savez(HERE / "free_snapshots_meta.npz",
             phi_pi=np.array([s[0] for s in snaps]),
             infid_K32=infid32, infid_Ksel=infidF)
    print("done.")


if __name__ == "__main__":
    main()
