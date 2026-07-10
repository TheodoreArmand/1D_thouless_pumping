#!/usr/bin/env python3
"""
R5 union selection with a WIDTH-ENRICHED pool + vectorized greedy.

Motivated by tmp/.../diagnose.py: the free evolved single-particle packet
spreads/breathes (rms 0.7 -> 1.6 -> 1.3) so the fixed 2-width structured pool
cannot span the plateau. Here the candidate pool sweeps a range of widths along
each particle's transit sub-path; the 32 existing K32 terms stay forced (P1
union-not-replacement). Everything is exact via 1D factorization of the free
2-particle state psi_m = sym[phiL(t_m) (x) phiR(t_m)].

Outputs the best K-term basis, the infidelity-vs-K curve, the infidelity-vs-phi
comparison (K32 vs chosen K), the Gram spectrum, initial-state refit fidelity,
and the union basis CSV -- all under tmp/, touching no production code.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "rice_mele_reference" / "lohes_experience"))
import n2_grid_reference as gr

# ---- config ----
A_LAT = 8.0; CELLS = 16; LENGTH = CELLS * A_LAT; NG = 2048
DX = LENGTH / NG; X = np.arange(NG) * DX - LENGTH / 2
P = 2 * np.pi * np.fft.fftfreq(NG, d=DX); P2 = P * P
ER = gr.recoil_energy(A_LAT); VS = 3.0 * ER; VL = 3.0 * ER
DT = 0.01; PERIOD = 160.0 * np.pi
SCHEDULE = REPO / "rice_mele_reference" / "Vs3Vl3_3_3" / "gap_adaptive_vs3vl3_full_depth_schedule.csv"
N1CSV = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N1_K16.csv"
K32CSV = REPO / "initial_state" / "Vs3Er_Vl3Er" / "initial_pathpad_N2_K32.csv"

PHI_PI = np.round(np.arange(0.30, 0.6001, 0.02), 4)   # 16 snapshots
GATE1_PHI = 0.52          # Gate-1 acceptance window phi <= 0.52 pi
K_MAX = 64                # explore up to here; report K=48 and early-stop K
INFID_STOP = 1e-3
NOVELTY_VETO = 5e-3
REG = 1e-10
U_FLOOR_FRAC = 1e-4

# width-enriched sweep (b = 1/(2 sigma^2); rms=sigma)
W_SET = [0.18, 0.22, 0.28, 0.35, 0.45, 0.574, 0.75, 1.0, 1.436, 2.196]
CL = list(np.round(np.arange(0.0, -4.01, -0.5), 3))   # left particle transit 0 -> -4 (half sites)
CR = list(np.round(np.arange(4.0, 8.01, 0.5), 3))     # right particle transit 8 -> 4 (half sites)


def gauss(b, r):
    return gr.gaussian(X, LENGTH, b, r)

def ov(f, h):
    return np.sum(np.conj(f) * h) * DX

def load_n1():
    rows = [l.strip() for l in N1CSV.read_text().splitlines() if l.strip()]
    return [(complex(*map(float, rows[5*k].split(","))).real,
             complex(*map(float, rows[5*k+2].split(","))).real,
             complex(*map(float, rows[5*k+3].split(","))).real)
            for k in range(len(rows)//5)]

def load_k32():
    rows = [l.strip() for l in K32CSV.read_text().splitlines() if l.strip()]
    st = 12; out = []
    for k in range(len(rows)//st):
        o = st*k
        out.append((complex(*map(float, rows[o].split(","))).real,
                    complex(*map(float, rows[o+5].split(","))).real,
                    complex(*map(float, rows[o+9].split(","))).real,
                    complex(*map(float, rows[o+8].split(","))).real,
                    complex(*map(float, rows[o+10].split(","))).real))
    return out


def evolve_snaps():
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
    snaps = []; ti = 0; t = 0.0; tg = list(PHI_PI)
    phin = gr.phi_at(t, PERIOD, ss, phis)
    while ti < len(tg):
        if phin/np.pi >= tg[ti]-1e-12:
            snaps.append((tg[ti], pL.copy(), pR.copy())); ti += 1; continue
        pL = step(pL, t); pR = step(pR, t); t += DT
        phin = gr.phi_at(t, PERIOD, ss, phis)
    return snaps


def canon(b0, r0, b1, r1):
    return tuple(sorted([(round(b0, 6), round(r0, 6)), (round(b1, 6), round(r1, 6))]))


def build_pool():
    pool = []; seen = {}
    def add(b0, r0, b1, r1, tier, forced=False, u0=0.0):
        k = canon(b0, r0, b1, r1)
        if k in seen:
            return
        seen[k] = len(pool)
        pool.append(dict(b0=b0, r0=r0, b1=b1, r1=r1, tier=tier, forced=forced, u0=u0))
    for (u, b0, r0, b1, r1) in load_k32():
        add(b0, r0, b1, r1, "P-A", forced=True, u0=u)
    for cl in CL:
        for wl in W_SET:
            for cr in CR:
                for wr in W_SET:
                    add(wl, cl, wr, cr, "rich")
    return pool


def overlap_tables(pool, snaps):
    distinct = {}
    def gidx(b, r):
        k = (round(b, 6), round(r, 6))
        if k not in distinct:
            distinct[k] = (len(distinct), gauss(b, r).real)
        return distinct[k][0]
    for t in pool:
        t["i0"] = gidx(t["b0"], t["r0"]); t["i1"] = gidx(t["b1"], t["r1"])
    nd = len(distinct); G = np.zeros((nd, NG))
    for (idx, arr) in distinct.values():
        G[idx] = arr
    ov11 = (G @ G.T) * DX
    nsnap = len(snaps)
    ovL = np.zeros((nd, nsnap), complex); ovR = np.zeros((nd, nsnap), complex)
    pn = np.zeros(nsnap)
    for m, (_pi, pL, pR) in enumerate(snaps):
        ovL[:, m] = (G @ pL) * DX; ovR[:, m] = (G @ pR) * DX
        nLL = np.real(ov(pL, pL)); nRR = np.real(ov(pR, pR)); nLR = ov(pL, pR)
        pn[m] = np.sqrt(nLL*nRR + np.abs(nLR)**2)
    i0 = np.array([t["i0"] for t in pool]); i1 = np.array([t["i1"] for t in pool])
    cn = np.sqrt(ov11[i0, i0]*ov11[i1, i1] + ov11[i0, i1]**2)
    S = (ov11[np.ix_(i0, i0)]*ov11[np.ix_(i1, i1)] + ov11[np.ix_(i0, i1)]*ov11[np.ix_(i1, i0)])
    S = S/np.outer(cn, cn); np.fill_diagonal(S, 1.0)
    B = (ovL[i0]*ovR[i1] + ovR[i0]*ovL[i1]) / cn[:, None] / pn[None, :]
    return S, B


def infid_set(S, B, I):
    SII = S[np.ix_(I, I)] + REG*np.eye(len(I))
    bI = B[np.ix_(I, range(B.shape[1]))]
    c = np.linalg.solve(SII, bI)
    proj = np.real(np.sum(np.conj(bI)*c, axis=0))
    return np.clip(1.0 - proj, 0.0, 1.0)


def greedy(S, B, pool, k_max):
    forced = [i for i, t in enumerate(pool) if t["forced"]]
    I = list(forced); nsnap = B.shape[1]
    ncand = len(pool)
    is_in = np.zeros(ncand, bool); is_in[I] = True
    curve = []; order = []
    infid = infid_set(S, B, I)
    curve.append((len(I), float(infid.max())))
    while len(I) < k_max and infid.max() > INFID_STOP:
        SII = S[np.ix_(I, I)] + REG*np.eye(len(I))
        SII_inv = np.linalg.inv(SII)
        c = SII_inv @ B[np.ix_(I, range(nsnap))]      # |I| x nsnap
        SjI = S[:, I]                                  # ncand x |I|
        proj = np.einsum('ji,jk->j', SjI @ SII_inv, SjI, optimize=True)
        nov = 1.0 - proj
        resov = B - (SjI @ c)                          # ncand x nsnap
        peak = np.max(np.abs(resov)**2, axis=1)
        score = peak / (nov + 1e-6)
        score[is_in] = -np.inf
        score[nov < NOVELTY_VETO] = -np.inf
        j = int(np.argmax(score))
        if not np.isfinite(score[j]):
            order.append(("veto-exhausted", None, float(infid.max()), None, None))
            break
        I.append(j); is_in[j] = True
        infid = infid_set(S, B, I)
        argsnap = int(np.argmax(np.abs(B[j, :])))
        order.append((pool[j]["tier"], j, float(infid.max()), float(nov[j]), argsnap))
        curve.append((len(I), float(infid.max())))
    return I, order, curve


def refit_initial(pool, I):
    n1 = load_n1()
    pL = np.zeros(NG); pR = np.zeros(NG)
    for u, b, r in n1:
        pL += u*gauss(b, r).real; pR += u*gauss(b, r+A_LAT).real
    pL /= np.sqrt(np.sum(pL*pL)*DX); pR /= np.sqrt(np.sum(pR*pR)*DX)
    nLL = np.sum(pL*pL)*DX; nRR = np.sum(pR*pR)*DX; nLR = np.sum(pL*pR)*DX
    tnorm = np.sqrt(nLL*nRR + nLR**2)
    sel = [pool[i] for i in I]
    G0 = np.array([gauss(t["b0"], t["r0"]).real for t in sel])
    G1 = np.array([gauss(t["b1"], t["r1"]).real for t in sel])
    gL0 = (G0 @ pL)*DX; gR0 = (G0 @ pR)*DX; gL1 = (G1 @ pL)*DX; gR1 = (G1 @ pR)*DX
    b_vec = (gL0*gR1 + gR0*gL1)/tnorm
    o00 = (G0@G0.T)*DX; o11 = (G1@G1.T)*DX; o01 = (G0@G1.T)*DX; o10 = o01.T
    Gram = o00*o11 + o01*o10
    u = np.linalg.solve(Gram + 1e-12*np.eye(len(sel)), b_vec)
    fit_n2 = float(u@Gram@u); overlap = float(u@b_vec); fid = overlap**2/fit_n2
    u = u/np.sqrt(fit_n2)
    umax = np.max(np.abs(u)); floor = U_FLOOR_FRAC*umax
    uf = np.where(np.abs(u) < floor, floor, u)
    return uf, fid, floor


def write_csv(path, pool, I, u):
    sel = [pool[i] for i in I]
    with open(path, "w") as f:
        for k, (t, uk) in enumerate(zip(sel, u)):
            f.write(f"{uk:.17g},0\n")
            for _ in range(4):
                f.write("0,0\n")
            f.write(f"{t['b0']:.17g},0\n"); f.write("0,0\n"); f.write("0,0\n")
            f.write(f"{t['b1']:.17g},0\n")
            f.write(f"{t['r0']:.17g},0\n"); f.write(f"{t['r1']:.17g},0\n")
            f.write(f"{k}\n")


def main():
    snaps = evolve_snaps()
    gate_mask = np.array([s[0] <= GATE1_PHI + 1e-9 for s in snaps])
    pool = build_pool()
    tiers = {}
    for t in pool:
        tiers[t["tier"]] = tiers.get(t["tier"], 0)+1
    print(f"pool size={len(pool)} tiers={tiers}  snapshots={len(snaps)}")
    S, B = overlap_tables(pool, snaps)

    forced = [i for i, t in enumerate(pool) if t["forced"]]
    infid32 = infid_set(S, B, forced)
    print(f"K32 worst infid: all={infid32.max():.3e}  "
          f"gate1(phi<=0.52)={infid32[gate_mask].max():.3e}")

    I, order, curve = greedy(S, B, pool, K_MAX)
    print(f"\ngreedy selected up to K={len(I)}")
    print("selection order (round: tier center-widths novelty peakphi -> worst_infid):")
    for step, (tier, j, worst, nov, asn) in enumerate(order):
        if j is None:
            print(f"  {step:2d}: {tier}  worst={worst:.3e}"); continue
        t = pool[j]; sp = snaps[asn][0]
        print(f"  {step:2d}: +{tier:4s} p0=(b{t['b0']:.3g},r{t['r0']:+.1f}) "
              f"p1=(b{t['b1']:.3g},r{t['r1']:+.1f}) nov={nov:.3f} "
              f"peakphi/pi={sp:.2f} -> worst={worst:.3e}")

    # infid vs K curve
    print("\ninfid-vs-K (worst over all phi, and over gate1 phi<=0.52):")
    # recompute per-K gate1 worst by re-running prefixes
    Iprefix = list(forced)
    seq = [j for (_t, j, _w, _n, _a) in order if j is not None]
    print(f"  K=32  all={infid32.max():.3e}  gate1={infid32[gate_mask].max():.3e}")
    milestones = {}
    for n, j in enumerate(seq, start=33):
        Iprefix.append(j)
        inf = infid_set(S, B, Iprefix)
        allw = inf.max(); g1 = inf[gate_mask].max()
        if n in (36, 40, 44, 48, 52, 56, 60, 64):
            lmin = np.linalg.eigvalsh(S[np.ix_(Iprefix, Iprefix)])[0]
            print(f"  K={n:3d} all={allw:.3e}  gate1={g1:.3e}  gram_lmin={lmin:.2e}")
        for thr, name in [(1e-2, "all_1e-2"), (1e-3, "all_1e-3")]:
            if name not in milestones and allw <= thr:
                milestones[name] = n
        for thr, name in [(1e-2, "gate1_1e-2"), (1e-3, "gate1_1e-3")]:
            if name not in milestones and g1 <= thr:
                milestones[name] = n
    print("milestones (K to reach threshold):", milestones)

    # choose K=48 as the requested target (or fewer if early-stop hit)
    Ksel = min(48, len(I))
    Isel = list(forced) + seq[:Ksel-32]
    infidSel = infid_set(S, B, Isel)
    print(f"\n=== CHOSEN K={len(Isel)} ===")
    print("phi/pi   infid_K32   infid_Ksel")
    for m, s in enumerate(snaps):
        print(f"  {s[0]:.2f}    {infid32[m]:.3e}   {infidSel[m]:.3e}")
    print(f"worst: all K32={infid32.max():.3e} Ksel={infidSel.max():.3e} | "
          f"gate1 K32={infid32[gate_mask].max():.3e} Ksel={infidSel[gate_mask].max():.3e}")

    Ssel = S[np.ix_(Isel, Isel)]
    eig = np.linalg.eigvalsh(Ssel)
    print(f"Gram spectrum: lambda_min={eig[0]:.3e} lambda_max={eig[-1]:.3e} cond={eig[-1]/max(eig[0],1e-300):.3e}")

    u, fid0, floor = refit_initial(pool, Isel)
    print(f"initial fidelity |<K{len(Isel)}(0)|target>|^2 = {fid0:.10f}  u_floor={floor:.2e}")
    outcsv = HERE / f"initial_pathpad_N2_K{len(Isel)}_richunion.csv"
    write_csv(outcsv, pool, Isel, u)
    print(f"wrote {outcsv}")

    summary = dict(pool_size=len(pool), tiers=tiers, K_selected=len(Isel),
                   phi_pi=[float(s[0]) for s in snaps],
                   infid_K32=[float(v) for v in infid32],
                   infid_Ksel=[float(v) for v in infidSel],
                   worst_all_K32=float(infid32.max()), worst_all_Ksel=float(infidSel.max()),
                   worst_gate1_K32=float(infid32[gate_mask].max()),
                   worst_gate1_Ksel=float(infidSel[gate_mask].max()),
                   curve=[(int(k), float(v)) for k, v in curve],
                   milestones=milestones,
                   gram_lambda_min=float(eig[0]), gram_lambda_max=float(eig[-1]),
                   init_fidelity=float(fid0), u_floor=float(floor),
                   selected=[dict(tier=pool[i]["tier"], b0=pool[i]["b0"], r0=pool[i]["r0"],
                                  b1=pool[i]["b1"], r1=pool[i]["r1"], forced=pool[i]["forced"])
                             for i in Isel])
    (HERE/"r5_rich_summary.json").write_text(json.dumps(summary, indent=2))
    print("wrote r5_rich_summary.json")


if __name__ == "__main__":
    main()
