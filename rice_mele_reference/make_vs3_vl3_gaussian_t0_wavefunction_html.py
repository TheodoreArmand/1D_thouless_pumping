#!/usr/bin/env python3
"""Build a screenshot-ready Vs=Vl=3 E_R Hamiltonian and t=0 state page.

The figure is reconstructed from the normalized K=32 ``basis_initial.csv``
written by the current Gate-1 Gaussian run.  For N=2 the physical wavefunction
is two-dimensional, so the one-dimensional overlay uses the one-body number
density n(x,0) = 2 int |Psi(x,x2;0)|^2 dx2 (integral = 2).
"""
from __future__ import annotations

import base64
import csv
import io
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
TASK = (
    REPO
    / "out/vs3_n2_gate1_shortwin_rcond1em5"
    / "gaussFrozenA_gate_sigma1p000"
    / "a8p000_K32_tmax175p929_VsER3p000_VlER3p000_full_depth_schedule"
)
BASIS = TASK / "basis_initial.csv"
PROGRESS = TASK / "progress.csv"
OUT_HTML = REPO / "vs3_vl3_gaussian_t0_wavefunction.html"

N = 2
K = 32
A_LAT = 8.0
HBAR = 1.0
MASS = 1.0
ER = HBAR**2 * (2.0 * np.pi / A_LAT) ** 2 / (2.0 * MASS)
G_OVER_ER = 0.3
SIGMA_G = 1.0
GRID_POINTS = 2048
GRID_RANGE = (-64.0, 64.0)
PLOT_RANGE = (-4.0, 12.0)


def parse_complex(row: str) -> complex:
    real, imag = map(float, row.split(","))
    return complex(real, imag)


def load_basis(path: Path) -> list[tuple[complex, np.ndarray, np.ndarray, np.ndarray]]:
    rows = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    stride = 1 + N * N + N * N + N + 1
    if len(rows) != K * stride:
        raise RuntimeError(f"expected {K * stride} nonempty rows in {path}, got {len(rows)}")

    terms = []
    for k in range(K):
        offset = k * stride
        u = parse_complex(rows[offset])
        cursor = offset + 1
        A = np.array(
            [parse_complex(rows[cursor + j]) for j in range(N * N)],
            dtype=np.complex128,
        ).reshape(N, N)
        cursor += N * N
        B = np.array(
            [parse_complex(rows[cursor + j]) for j in range(N * N)],
            dtype=np.complex128,
        ).reshape(N, N)
        cursor += N * N
        R = np.array(
            [parse_complex(rows[cursor + j]) for j in range(N)],
            dtype=np.complex128,
        )
        terms.append((u, A, B, R))

    # These assertions lock the reconstruction to the actual current K32 state.
    max_imag = max(
        max(abs(u.imag), np.max(np.abs(A.imag)), np.max(np.abs(B.imag)), np.max(np.abs(R.imag)))
        for u, A, B, R in terms
    )
    max_A = max(np.max(np.abs(A)) for _u, A, _B, _R in terms)
    max_B_offdiag = max(
        max(abs(B[0, 1]), abs(B[1, 0])) for _u, _A, B, _R in terms
    )
    if max_imag > 1.0e-14 or max_A > 1.0e-14 or max_B_offdiag > 1.0e-14:
        raise RuntimeError("current t=0 basis is no longer real with A=0 and diagonal B")
    return terms


def reconstruct_state(terms):
    x0, x1 = GRID_RANGE
    dx = (x1 - x0) / GRID_POINTS
    x = x0 + np.arange(GRID_POINTS) * dx
    psi = np.zeros((GRID_POINTS, GRID_POINTS), dtype=np.complex128)
    cache: dict[tuple[float, float], np.ndarray] = {}

    def gaussian(width: float, center: float) -> np.ndarray:
        key = (width, center)
        if key not in cache:
            # The ECG engine integrates Gaussians on the real line: no ring wrap.
            cache[key] = np.exp(-width * (x - center) ** 2)
        return cache[key]

    for u, _A, B, R in terms:
        g0 = gaussian(float(B[0, 0].real), float(R[0].real))
        g1 = gaussian(float(B[1, 1].real), float(R[1].real))
        psi += u * (
            g0[:, None] * g1[None, :] + g1[:, None] * g0[None, :]
        ) / np.sqrt(2.0)

    prob = np.abs(psi) ** 2
    norm = float(prob.sum() * dx * dx)
    if abs(norm - 1.0) > 1.0e-11:
        raise RuntimeError(f"grid reconstruction norm mismatch: {norm:.17g}")
    prob /= norm
    density = (prob.sum(axis=0) + prob.sum(axis=1)) * dx
    density_sum = float(density.sum() * dx)
    if abs(density_sum - N) > 2.0e-12:
        raise RuntimeError(f"one-body density sum-rule mismatch: {density_sum:.17g}")

    x_mean_sum = float((x * density).sum() * dx)
    x_diff = x[:, None] - x[None, :]
    r12_sq = float((x_diff * x_diff * prob).sum() * dx * dx)
    v_gauss = float(
        (G_OVER_ER * ER * np.exp(-(x_diff * x_diff) / SIGMA_G**2) * prob).sum()
        * dx
        * dx
    )
    return x, density, norm, x_mean_sum, r12_sq, v_gauss


def validate_against_run(norm: float, x_mean_sum: float, r12_sq: float, v_gauss: float) -> None:
    with PROGRESS.open(newline="") as stream:
        row0 = next(csv.DictReader(stream))
    expected = {
        "norm": float(row0["norm"]),
        "x_mean": float(row0["x_mean"]),
        "r12_sq": float(row0["r12_sq"]),
        "V_gauss": float(row0["V_gauss"]),
    }
    actual = {
        "norm": norm,
        "x_mean": x_mean_sum,
        "r12_sq": r12_sq,
        "V_gauss": v_gauss,
    }
    tolerances = {"norm": 2e-12, "x_mean": 2e-11, "r12_sq": 2e-10, "V_gauss": 2e-18}
    for key in expected:
        error = abs(actual[key] - expected[key])
        if error > tolerances[key]:
            raise RuntimeError(
                f"{key} mismatch against t=0 C++ progress: "
                f"grid={actual[key]:.17g}, C++={expected[key]:.17g}, error={error:.3e}"
            )


def make_figure(x: np.ndarray, density: np.ndarray) -> str:
    mask = (x >= PLOT_RANGE[0]) & (x <= PLOT_RANGE[1])
    xp = x[mask]
    np_density = density[mask]
    potential = -3.0 * np.cos(4.0 * np.pi * xp / A_LAT) - 3.0 * np.cos(
        2.0 * np.pi * xp / A_LAT
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11.5,
            "axes.titleweight": "semibold",
            "axes.edgecolor": "#334155",
            "axes.labelcolor": "#0f172a",
            "xtick.color": "#334155",
            "ytick.color": "#334155",
        }
    )
    fig, ax_v = plt.subplots(figsize=(12.0, 5.4), dpi=220)
    fig.patch.set_facecolor("white")
    ax_v.set_facecolor("#fbfdff")

    potential_line = ax_v.plot(
        xp,
        potential,
        color="#b45309",
        lw=2.35,
        label=r"$V_{\mathrm{lat}}(x,0)/E_R$",
        zorder=3,
    )[0]
    ax_v.axhline(0.0, color="#94a3b8", lw=0.8, zorder=1)
    ax_v.set_xlim(PLOT_RANGE)
    ax_v.set_ylim(-6.8, 4.25)
    ax_v.set_xlabel(r"position $x$  (code length)")
    ax_v.set_ylabel(r"lattice potential  $V_{\mathrm{lat}}/E_R$", color="#92400e")
    ax_v.tick_params(axis="y", colors="#92400e")
    ax_v.grid(axis="x", color="#cbd5e1", lw=0.65, alpha=0.7)
    ax_v.set_xticks(np.arange(-4.0, 12.1, 2.0))

    for center, label in [(0.0, "deep well"), (4.0, "shallow well"), (8.0, "deep well")]:
        ax_v.axvline(center, color="#64748b", lw=0.85, ls=(0, (2, 3)), alpha=0.7)
        y_text = -6.35 if center != 4.0 else 0.38
        ax_v.text(
            center,
            y_text,
            label,
            ha="center",
            va="bottom",
            color="#475569",
            fontsize=9.2,
        )

    ax_n = ax_v.twinx()
    density_fill = ax_n.fill_between(
        xp,
        0.0,
        np_density,
        color="#2563eb",
        alpha=0.23,
        label=r"$n(x,0)$",
        zorder=2,
    )
    density_line = ax_n.plot(xp, np_density, color="#1d4ed8", lw=2.15, zorder=4)[0]
    ax_n.set_ylim(0.0, max(0.78, float(np_density.max()) * 1.14))
    ax_n.set_ylabel(r"one-body number density  $n(x,0)$", color="#1d4ed8")
    ax_n.tick_params(axis="y", colors="#1d4ed8")

    ax_v.set_title(
        r"Current $N=2$, $K=32$ ECG initial state on the $V_s=V_l=3E_R$ superlattice",
        pad=13,
        fontsize=14.0,
    )
    ax_v.legend(
        [potential_line, (density_fill, density_line)],
        [r"$V_{\mathrm{lat}}(x,0)/E_R$", r"$n(x,0)$  (integral $=2$)"],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        frameon=True,
        framealpha=0.94,
        edgecolor="#cbd5e1",
        ncol=2,
    )
    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def build_html(image_b64: str, norm: float, x_mean_sum: float, r12_sq: float, v_gauss: float) -> str:
    return rf'''<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Vs=Vl=3ER Gaussian 两体哈密顿量与 t=0 初态</title>
  <script>
  window.MathJax = {{
    tex: {{ inlineMath: [['\\(', '\\)']], displayMath: [['\\[', '\\]']] }},
    svg: {{ fontCache: 'global' }}
  }};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
  <style>
    :root {{
      --ink: #0f172a;
      --muted: #64748b;
      --line: #dbe4ef;
      --blue: #1d4ed8;
      --teal: #0f766e;
      --paper: #ffffff;
      --wash: #f8fafc;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: #eef2f7;
      color: var(--ink);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
      line-height: 1.55;
    }}
    main {{
      width: min(1180px, calc(100vw - 32px));
      margin: 24px auto;
      padding: 34px 40px 28px;
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 16px;
      box-shadow: 0 16px 44px rgba(15, 23, 42, 0.08);
    }}
    h1 {{ margin: 0; font-size: 1.8rem; letter-spacing: -0.02em; }}
    .lead {{ margin: 7px 0 22px; color: var(--muted); font-size: 1rem; }}
    .tags {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 0 0 20px; }}
    .tag {{
      padding: 4px 10px;
      border: 1px solid #bfdbfe;
      border-radius: 999px;
      background: #eff6ff;
      color: #1e40af;
      font-size: 0.84rem;
      font-weight: 650;
    }}
    .equation {{
      margin: 13px 0 18px;
      padding: 16px 20px;
      border: 1px solid var(--line);
      border-left: 5px solid var(--teal);
      border-radius: 10px;
      background: #fbfcfd;
      font-size: 1.02rem;
      overflow-x: auto;
      white-space: normal;
      line-height: 1.45;
    }}
    .equation h2 {{ margin: 0 0 7px; font-size: 1rem; color: #115e59; }}
    .equation mjx-container[jax="SVG"][display="true"] {{ margin: 0.4em 0; }}
    .equation .note {{
      margin-top: 8px;
      color: var(--muted);
      font-size: 0.9rem;
    }}
    figure {{
      margin: 23px 0 0;
      padding: 12px 12px 10px;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: var(--wash);
    }}
    figure img {{ display: block; width: 100%; height: auto; border-radius: 7px; background: white; }}
    figcaption {{ margin: 8px 8px 2px; color: #475569; font-size: 0.91rem; }}
    .audit {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 9px;
      margin-top: 17px;
    }}
    .audit div {{ padding: 9px 11px; border: 1px solid var(--line); border-radius: 8px; background: #fff; }}
    .audit span {{ display: block; color: var(--muted); font-size: 0.73rem; }}
    .audit b {{ font-size: 0.9rem; font-variant-numeric: tabular-nums; }}
    footer {{ margin-top: 18px; color: var(--muted); font-size: 0.78rem; }}
    code {{ padding: 1px 4px; border-radius: 4px; background: #eef2f7; }}
    @media (max-width: 760px) {{
      main {{ padding: 24px 18px; }}
      .audit {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .equation {{ font-size: 0.93rem; padding: 14px 13px; }}
    }}
    @media print {{
      body {{ background: white; }}
      main {{ width: 100%; margin: 0; border: 0; box-shadow: none; }}
    }}
  </style>
</head>
<body>
<main>
  <h1>3/3 超晶格中的两体 Gaussian Thouless 泵</h1>
  <p class="lead">当前 full-depth 协议；公式与图均对应正在使用的 \(N=2,\ K=32,\ \sigma_{{\mathrm G}}=1\) Gaussian 分支。</p>
  <div class="tags">
    <span class="tag">two bosons</span>
    <span class="tag">full-depth 3/3</span>
    <span class="tag">K = 32 ECG</span>
    <span class="tag">t = 0, φ = 0</span>
  </div>

  <section class="equation">
    <h2>1 · \(V_s=V_l=3E_R\) 的两体哈密顿量</h2>
    \[
      \boxed{{
      \hat H(t)=
      \sum_{{j=1}}^{{2}}
      \left[
        -\frac{{\hbar^2}}{{2m}}\frac{{\partial^2}}{{\partial x_j^2}}
        -3E_R\cos\!\left(\frac{{4\pi x_j}}{{a}}\right)
        -3E_R\cos\!\left(\frac{{2\pi x_j}}{{a}}+\phi(t)\right)
      \right]
      +V_{{\mathrm G}}(x_1-x_2)
      }}
    \]
    <div class="note">
      这里 \(a=8,\ \hbar=m=1\)，且 \(E_R=\hbar^2(2\pi/a)^2/(2m)=\pi^2/32\)。短、长晶格周期分别为 \(a/2=4\) 与 \(a=8\)；当前相位采用 \(+\phi(t)\) 约定。
    </div>
  </section>

  <section class="equation">
    <h2>2 · 实际使用的 Gaussian 两体相互作用</h2>
    \[
      \boxed{{
      V_{{\mathrm G}}(x_1-x_2)
      =g_{{\mathrm G}}
      \exp\!\left[-\frac{{(x_1-x_2)^2}}{{\sigma_{{\mathrm G}}^2}}\right],
      \qquad
      g_{{\mathrm G}}=0.3E_R,
      \quad \sigma_{{\mathrm G}}=1
      }}
    \]
    <div class="note">
      这里 \(\sigma_{{\mathrm G}}=1\) 是 code length（即 \(a/8\)）。这是未归一化的 Gaussian 势：没有 \(1/(\sqrt{{\pi}}\sigma_{{\mathrm G}})\) 前因子，指数分母也不是 \(2\sigma_{{\mathrm G}}^2\)。因此 \(g_{{\mathrm G}}\) 就是两粒子重合时的峰值势能。
    </div>
  </section>

  <figure>
    <img src="data:image/png;base64,{image_b64}" alt="N=2 K=32 t=0 one-body density over the Vs=Vl=3ER lattice potential">
    <figcaption>
      蓝色为完整两体 ECG 波函数导出的一体粒子数密度 \(n(x,0)=2\int |\Psi(x,x_2;0)|^2\,dx_2\)，其积分为 2；棕色为 \(t=0\)、\(\phi(0)=0\) 时的真实 3/3 晶格势。两体波函数本身依赖 \((x_1,x_2)\)，因此这里没有把它误标成一维的 \(|\psi(x)|^2\)。
    </figcaption>
  </figure>

  <div class="audit" aria-label="numerical validation">
    <div><span>grid reconstruction</span><b>‖Ψ(0)‖² = {norm:.12f}</b></div>
    <div><span>density sum rule</span><b>∫ n(x,0) dx = 2</b></div>
    <div><span>C++ matched ⟨x₁+x₂⟩</span><b>{x_mean_sum:.12f}</b></div>
    <div><span>C++ matched ⟨V_G⟩ / E_R</span><b>{v_gauss / ER:.6e}</b></div>
  </div>

  <footer>
    数据源：当前 Gate-1 Gaussian 运行写出的归一化 <code>basis_initial.csv</code>；Gaussian 与 free 分支共用这个 K32 初态，它不是重新优化的相互作用基态。图使用全实线 ECG Gaussian，不做 periodic/minimum-image 包裹。重构同时核对 \(\langle(x_1-x_2)^2\rangle={r12_sq:.12f}\)。
  </footer>
</main>
</body>
</html>
'''


def main() -> None:
    if not BASIS.exists() or not PROGRESS.exists():
        raise RuntimeError(f"missing current run data under {TASK}")
    terms = load_basis(BASIS)
    x, density, norm, x_mean_sum, r12_sq, v_gauss = reconstruct_state(terms)
    validate_against_run(norm, x_mean_sum, r12_sq, v_gauss)
    image_b64 = make_figure(x, density)
    html = build_html(image_b64, norm, x_mean_sum, r12_sq, v_gauss)
    tmp = OUT_HTML.with_suffix(OUT_HTML.suffix + ".tmp")
    tmp.write_text(html)
    tmp.replace(OUT_HTML)
    print(f"wrote {OUT_HTML}")
    print(
        "validated: "
        f"norm={norm:.16g}, int_density=2, x_mean_sum={x_mean_sum:.16g}, "
        f"r12_sq={r12_sq:.16g}, V_gauss/Er={v_gauss / ER:.16g}"
    )


if __name__ == "__main__":
    main()
