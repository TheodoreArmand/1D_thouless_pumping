#!/usr/bin/env python3
"""Build an HTML comparison of the current Vs3/Vl3 N=2 ECG and grid runs."""
from __future__ import annotations

import base64
import csv
import io
import math
import os
import re
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
OUT_HTML = REPO / "vs3_n2_ecg_grid_comparison.html"
TOTAL_STEPS = 50266
TOTAL_TIME = 502.6548245743669
IDEAL_DELTA_P = -2.0

CASES = {
    "free": {
        "title": "Free, no pair interaction",
        "job": "603717",
        "grid_job": "604164_0",
        "ecg": REPO / "out/vs3_n2_dt0p01_T160pi/free/a8p000_K24_tmax502p655_VsER3p000_VlER3p000/progress.csv",
        "grid_csv": REPO / "out/vs3_n2_dt0p01_T160pi/grid_ref/free/g1024_dt0p01_T160pi/n2_grid_reference_vs3_free_g1024_dt0p01_T160pi.csv",
        "grid_log": REPO / "slurm/vs3_n2_dt0p01_T160pi/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/vs3n2grid_604164_0.out",
        "color_ecg": "#2267a8",
        "color_grid": "#111111",
    },
    "gauss": {
        "title": "Gaussian pair interaction, g = 0.3 Er",
        "job": "603718",
        "grid_job": "604164_1",
        "ecg": REPO / "out/vs3_n2_dt0p01_T160pi/gauss/a8p000_K24_tmax502p655_VsER3p000_VlER3p000/progress.csv",
        "grid_csv": REPO / "out/vs3_n2_dt0p01_T160pi/grid_ref/gauss/g1024_dt0p01_T160pi/n2_grid_reference_vs3_gauss_g1024_dt0p01_T160pi.csv",
        "grid_log": REPO / "slurm/vs3_n2_dt0p01_T160pi/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/vs3n2grid_604164_1.out",
        "color_ecg": "#b03a2e",
        "color_grid": "#111111",
    },
}

PROGRESS_RE = re.compile(
    r"progress grid N=2 step (?P<step>\d+)/(?P<total>\d+) "
    r"(?P<pct>[0-9.]+)% t=(?P<t>[-+0-9.eE]+) eta_s=(?P<eta>[-+0-9.eE]+) "
    r"P=(?P<P>[-+0-9.eE]+) dP=(?P<dP>[-+0-9.eE]+) "
    r"r12=(?P<r12>[-+0-9.eE]+) Vg=(?P<Vg>[-+0-9.eE]+) "
    r"norm=(?P<norm>[-+0-9.eE]+)"
)


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def fnum(x: float, digits: int = 6) -> str:
    if x is None or not math.isfinite(x):
        return "n/a"
    if x == 0:
        return "0"
    if abs(x) < 1e-4 or abs(x) >= 1e5:
        return f"{x:.{digits}e}"
    return f"{x:.{digits}f}"


def pct(row: dict) -> float:
    return 100.0 * row["step"] / TOTAL_STEPS


def read_ecg(path: Path) -> list[dict]:
    rows = []
    with path.open() as f:
        for row in csv.DictReader(f):
            rows.append({
                "step": int(float(row["step"])),
                "t": float(row["t"]),
                "phi": float(row["phi"]),
                "P": float(row["polarization_cell"]),
                "dP": float(row["delta_polarization"]),
                "r12": float(row["r12_rms"]),
                "Vg": float(row["V_gauss"]),
                "norm": float(row["norm"]),
                "resid": float(row["relative_raw_residual"]),
                "minAB": float(row["min_re_AplusB"]),
                "minB": float(row["min_re_B"]),
                "rank": int(float(row["actual_solve_rank"])),
                "eta": float(row["eta_seconds"]),
            })
    return rows


def read_grid_log(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        m = PROGRESS_RE.search(line)
        if not m:
            continue
        item = {k: float(v) for k, v in m.groupdict().items()
                if k not in {"step", "total"}}
        item["step"] = int(m.group("step"))
        item["total"] = int(m.group("total"))
        rows.append(item)
    return rows


def read_grid_csv(path: Path) -> list[dict]:
    rows = []
    with path.open() as f:
        for row in csv.DictReader(f):
            t = float(row["t"])
            step = TOTAL_STEPS if abs(t - TOTAL_TIME) < 1e-9 else int(round(t / 0.01))
            rows.append({
                "step": step,
                "total": TOTAL_STEPS,
                "pct": 100.0 * t / TOTAL_TIME,
                "t": t,
                "phi": float(row["phi"]),
                "eta": float("nan"),
                "P": float(row["polarization_cell"]),
                "dP": float(row["delta_polarization"]),
                "r12": float(row["r12_rms"]),
                "Vg": float(row["V_gauss"]),
                "norm": float(row["norm"]),
            })
    return rows


def read_grid(cfg: dict) -> tuple[list[dict], Path, str]:
    csv_path = cfg.get("grid_csv")
    if isinstance(csv_path, Path) and csv_path.exists():
        return read_grid_csv(csv_path), csv_path, "CSV"
    log_path = cfg["grid_log"]
    return read_grid_log(log_path), log_path, "Slurm log"


def common_by_step(ecg_rows: list[dict], grid_rows: list[dict]) -> list[dict]:
    ecg_map = {r["step"]: r for r in ecg_rows}
    grid_map = {r["step"]: r for r in grid_rows}
    common = []
    for step in sorted(set(ecg_map) & set(grid_map)):
        e = ecg_map[step]
        g = grid_map[step]
        common.append({
            "step": step,
            "t": e["t"],
            "phi": e.get("phi", g.get("phi", float("nan"))),
            "x": e["t"] / TOTAL_TIME,
            "ecg_P": e["P"],
            "grid_P": g["P"],
            "diff_P": e["P"] - g["P"],
            "ecg_dP": e["dP"],
            "grid_dP": g["dP"],
            "diff_dP": e["dP"] - g["dP"],
            "ecg_r12": e["r12"],
            "grid_r12": g["r12"],
            "diff_r12": e["r12"] - g["r12"],
            "ecg_resid": e["resid"],
            "ecg_minAB": e["minAB"],
            "ecg_rank": e["rank"],
        })
    return common


def max_abs(rows: list[dict], key: str) -> tuple[float, dict | None]:
    if not rows:
        return float("nan"), None
    row = max(rows, key=lambda r: abs(r[key]))
    return abs(row[key]), row


def first_crossing(rows: list[dict], key: str, threshold: float) -> dict | None:
    for row in rows:
        value = row[key]
        if math.isfinite(value) and value < threshold:
            return row
    return None


def fig_to_data_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=155, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def style_axis(ax) -> None:
    ax.grid(True, color="#e2e8f0", lw=0.65)
    ax.set_axisbelow(True)
    ax.tick_params(colors="#334155", labelsize=7.8, width=0.7)
    for spine in ax.spines.values():
        spine.set_color("#cbd5e1")
        spine.set_linewidth(0.8)


def plot_delta(case_name: str, cfg: dict, ecg_rows: list[dict], grid_rows: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(3.75, 2.45))
    ex = [r["t"] / TOTAL_TIME for r in ecg_rows]
    gx = [r["t"] / TOTAL_TIME for r in grid_rows]
    ax.plot(gx, [r["dP"] for r in grid_rows], color="#475569", lw=1.45, ls="--", label="grid")
    ax.plot(ex, [r["dP"] for r in ecg_rows], color=cfg["color_ecg"], lw=1.55, label="ECG")
    ax.axhline(IDEAL_DELTA_P, color="#94a3b8", lw=0.9, ls=":")
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"$t/T$", fontsize=8.3)
    ax.set_ylabel(r"$\Delta P$", fontsize=8.3)
    ax.set_title(f"{case_name}: displacement", fontsize=9.2)
    style_axis(ax)
    ax.legend(frameon=False, loc="best", fontsize=7.4, handlelength=2.2)
    fig.tight_layout()
    return fig_to_data_uri(fig)


def plot_delta_phase(case_name: str, cfg: dict, common: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(3.75, 2.45))
    phase = [r["phi"] / math.pi for r in common if math.isfinite(r["phi"])]
    ecg_dp = [r["ecg_dP"] for r in common if math.isfinite(r["phi"])]
    grid_dp = [r["grid_dP"] for r in common if math.isfinite(r["phi"])]
    ax.plot(phase, grid_dp, color="#475569", lw=1.45, ls="--", label="grid")
    ax.plot(phase, ecg_dp, color=cfg["color_ecg"], lw=1.55, label="ECG")
    ax.axhline(IDEAL_DELTA_P, color="#94a3b8", lw=0.9, ls=":")
    ax.set_xlim(0.0, 2.0)
    ax.set_xlabel(r"$\phi/\pi$", fontsize=8.3)
    ax.set_ylabel(r"$\Delta P$", fontsize=8.3)
    ax.set_title(f"{case_name}: phase axis", fontsize=9.2)
    style_axis(ax)
    ax.legend(frameon=False, loc="best", fontsize=7.4, handlelength=2.2)
    fig.tight_layout()
    return fig_to_data_uri(fig)


def plot_difference(cfg: dict, common: list[dict]) -> str:
    fig, axes = plt.subplots(2, 1, figsize=(3.75, 2.65), sharex=True)
    x = [r["x"] for r in common]
    axes[0].plot(x, [r["diff_dP"] for r in common], color=cfg["color_ecg"], lw=1.35)
    axes[0].axhline(0.0, color="#64748b", lw=0.75)
    axes[0].set_ylabel(r"$\delta\Delta P$", fontsize=8.0)
    axes[0].set_title("ECG - grid at common steps", fontsize=9.2)
    axes[1].plot(x, [r["diff_r12"] for r in common], color="#0d9488", lw=1.35)
    axes[1].axhline(0.0, color="#64748b", lw=0.75)
    axes[1].set_ylabel(r"$\delta r_{12}$", fontsize=8.0)
    axes[1].set_xlabel(r"$t/T$", fontsize=8.3)
    for ax in axes:
        ax.set_xlim(0.0, 1.0)
        style_axis(ax)
    fig.tight_layout(h_pad=0.35)
    return fig_to_data_uri(fig)


def plot_health(case_name: str, cfg: dict, ecg_rows: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(3.75, 2.45))
    ex = [r["t"] / TOTAL_TIME for r in ecg_rows]
    min_ab = [r["minAB"] if r["minAB"] > 0 else float("nan") for r in ecg_rows]
    resid = [r["resid"] for r in ecg_rows]
    ln1 = ax.plot(ex, min_ab, color="#334155", lw=1.35, label=r"$\min\mathrm{Re}(A+B)$")
    ax.axhline(0.0227284, color="#d97706", lw=0.85, ls="--", label="successful-case min")
    ax.set_yscale("log")
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"$t/T$", fontsize=8.3)
    ax.set_ylabel("width", fontsize=8.0)
    ax.set_title(f"{case_name}: ECG health", fontsize=9.2)
    style_axis(ax)

    ax2 = ax.twinx()
    ln2 = ax2.plot(ex, resid, color=cfg["color_ecg"], lw=1.15, alpha=0.9, label="raw residual")
    ax2.set_ylabel("residual", fontsize=8.0, color=cfg["color_ecg"])
    ax2.tick_params(colors=cfg["color_ecg"], labelsize=7.8, width=0.7)
    for spine in ax2.spines.values():
        spine.set_color("#cbd5e1")
        spine.set_linewidth(0.8)
    lines = ln1 + ln2
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, frameon=False, loc="best", fontsize=7.0, handlelength=1.8)
    fig.tight_layout()
    return fig_to_data_uri(fig)


def build_case(case_name: str, cfg: dict) -> dict:
    ecg_rows = read_ecg(cfg["ecg"])
    grid_rows, grid_source, grid_source_kind = read_grid(cfg)
    common = common_by_step(ecg_rows, grid_rows)
    latest_ecg = ecg_rows[-1]
    latest_grid = grid_rows[-1]
    latest_common = common[-1] if common else None
    max_dp, max_dp_row = max_abs(common, "diff_dP")
    max_p, max_p_row = max_abs(common, "diff_P")
    max_r12, max_r12_row = max_abs(common, "diff_r12")
    unhealthy = first_crossing(ecg_rows, "minAB", 0.0227284)
    severe = first_crossing(ecg_rows, "minAB", 1e-3)
    return {
        "name": case_name,
        "cfg": cfg,
        "ecg": ecg_rows,
        "grid": grid_rows,
        "grid_source": grid_source,
        "grid_source_kind": grid_source_kind,
        "common": common,
        "latest_ecg": latest_ecg,
        "latest_grid": latest_grid,
        "latest_common": latest_common,
        "max_dp": (max_dp, max_dp_row),
        "max_p": (max_p, max_p_row),
        "max_r12": (max_r12, max_r12_row),
        "unhealthy": unhealthy,
        "severe": severe,
        "delta_plot": plot_delta(case_name, cfg, ecg_rows, grid_rows),
        "delta_phase_plot": plot_delta_phase(case_name, cfg, common) if common else "",
        "diff_plot": plot_difference(cfg, common) if common else "",
        "health_plot": plot_health(case_name, cfg, ecg_rows),
    }


def summary_table(cases: list[dict]) -> str:
    rows = []
    for item in cases:
        lc = item["latest_common"]
        le = item["latest_ecg"]
        lg = item["latest_grid"]
        max_dp, max_dp_row = item["max_dp"]
        rows.append(f"""
          <tr>
            <td>{item['cfg']['title']}</td>
            <td><code>{item['cfg']['job']}</code></td>
            <td>{le['step']}/{TOTAL_STEPS} ({pct(le):.1f}%)</td>
            <td>{lg['step']}/{TOTAL_STEPS} ({pct(lg):.1f}%)</td>
            <td>{lc['step'] if lc else 'n/a'}</td>
            <td class="num">{fnum(lc['ecg_dP'] if lc else float('nan'))}</td>
            <td class="num">{fnum(lc['grid_dP'] if lc else float('nan'))}</td>
            <td class="num">{fnum(lc['diff_dP'] if lc else float('nan'))}</td>
            <td class="num">{fnum(max_dp)}</td>
            <td class="num">{fnum(le['minAB'])}</td>
            <td class="num">{fnum(le['resid'])}</td>
          </tr>
        """)
    return "\n".join(rows)


def crossing_text(row: dict | None) -> str:
    if row is None:
        return "not reached"
    return f"step {row['step']}, t/T={row['t'] / TOTAL_TIME:.3f}, value={fnum(row['minAB'])}"


def progress_bar(label: str, value: float, color: str = "#0d9488") -> str:
    clipped = min(100.0, max(0.0, value))
    return f"""
      <div class="progline">
        <span>{label}</span>
        <div class="prog"><span style="width:{clipped:.1f}%;background:{color}"></span></div>
        <b>{clipped:.1f}%</b>
      </div>
    """


def health_class(row: dict) -> str:
    if row["minAB"] < 1e-3:
        return "bad-banner"
    if row["minAB"] < 0.0227284:
        return "warn"
    return "ok"


def kpi_cards(items: list[dict]) -> str:
    cards = []
    for item in items:
        name = item["name"]
        le = item["latest_ecg"]
        lg = item["latest_grid"]
        lc = item["latest_common"]
        diff = fnum(lc["diff_dP"]) if lc else "n/a"
        cards.append(
            f"<div><b>{pct(le):.1f}% / {pct(lg):.1f}%</b>{name} ECG/grid done</div>"
        )
        cards.append(
            f"<div><b>{diff}</b>{name} latest \\(\\delta\\Delta P\\)</div>"
        )
    return "\n".join(cards)


def case_section(item: dict) -> str:
    cfg = item["cfg"]
    le = item["latest_ecg"]
    lg = item["latest_grid"]
    lc = item["latest_common"]
    max_dp, max_dp_row = item["max_dp"]
    max_p, max_p_row = item["max_p"]
    max_r12, max_r12_row = item["max_r12"]
    if lc:
        common_summary = f"""
        <table>
          <tr><th>latest common step</th><td>{lc['step']} ({100.0 * lc['step'] / TOTAL_STEPS:.1f}% of cycle)</td></tr>
          <tr><th>ECG \\(\\Delta P\\)</th><td class="num">{fnum(lc['ecg_dP'])}</td></tr>
          <tr><th>grid \\(\\Delta P\\)</th><td class="num">{fnum(lc['grid_dP'])}</td></tr>
          <tr><th>ECG-grid \\(\\Delta P\\)</th><td class="num">{fnum(lc['diff_dP'])}</td></tr>
          <tr><th>ECG-grid \\(r_{{12}}\\)</th><td class="num">{fnum(lc['diff_r12'])}</td></tr>
          <tr><th>max \\(|\\Delta P_{{\\mathrm{{ECG}}}}-\\Delta P_{{\\mathrm{{grid}}}}|\\)</th><td class="num">{fnum(max_dp)} at step {max_dp_row['step'] if max_dp_row else 'n/a'}</td></tr>
          <tr><th>max \\(|P_{{\\mathrm{{ECG}}}}-P_{{\\mathrm{{grid}}}}|\\)</th><td class="num">{fnum(max_p)} at step {max_p_row['step'] if max_p_row else 'n/a'}</td></tr>
          <tr><th>max \\(|r_{{12,\\mathrm{{ECG}}}}-r_{{12,\\mathrm{{grid}}}}|\\)</th><td class="num">{fnum(max_r12)} at step {max_r12_row['step'] if max_r12_row else 'n/a'}</td></tr>
        </table>
        """
    else:
        common_summary = "<p>No common sampled steps yet.</p>"
    diff_status = (
        f"At the latest common step {lc['step']}, "
        f"\\(\\Delta P_{{\\mathrm{{ECG}}}}-\\Delta P_{{\\mathrm{{grid}}}}={fnum(lc['diff_dP'])}\\)."
        if lc else
        "No common ECG/grid step has been logged yet."
    )
    return f"""
      <section class="section">
        <h2>{cfg['title']}</h2>
        <div class="{health_class(le)}">
          <b>Current read.</b> ECG is at {pct(le):.1f}% of the cycle, grid is at {pct(lg):.1f}%.
          {diff_status}
          ECG latest \\(\\min \\mathrm{{Re}}(A+B)={fnum(le['minAB'])}\\), residual {fnum(le['resid'])}.
        </div>
        <div class="progpair">
          {progress_bar('ECG progress', pct(le), cfg['color_ecg'])}
          {progress_bar('grid progress', pct(lg), '#475569')}
        </div>
        <div class="grid2">
          <div>
            <h3>Run State From Current Files</h3>
            <table>
              <tr><th>ECG job</th><td><code>{cfg['job']}</code></td></tr>
              <tr><th>grid job</th><td><code>{cfg['grid_job']}</code></td></tr>
              <tr><th>ECG latest</th><td>step {le['step']}/{TOTAL_STEPS}, t/T={le['t'] / TOTAL_TIME:.3f}, \\(\\Delta P={fnum(le['dP'])}\\)</td></tr>
              <tr><th>grid latest</th><td>step {lg['step']}/{TOTAL_STEPS}, t/T={lg['t'] / TOTAL_TIME:.3f}, \\(\\Delta P={fnum(lg['dP'])}\\)</td></tr>
              <tr><th>ECG latest residual</th><td class="num">{fnum(le['resid'])}</td></tr>
              <tr><th>ECG latest \\(\\min \\mathrm{{Re}}(A+B)\\)</th><td class="num">{fnum(le['minAB'])}</td></tr>
              <tr><th>ECG latest rank</th><td class="num">{le['rank']}</td></tr>
              <tr><th>first below successful-case width</th><td>{crossing_text(item['unhealthy'])}</td></tr>
              <tr><th>first below \\(10^{{-3}}\\)</th><td>{crossing_text(item['severe'])}</td></tr>
            </table>
          </div>
          <div>
            <h3>ECG vs Grid At Common Steps</h3>
            {common_summary}
          </div>
        </div>
        <div class="figgrid">
          <figure>
            <img src="{item['delta_plot']}" alt="{cfg['title']} displacement comparison">
            <figcaption><b>Displacement.</b> ECG trajectory against the split-step grid reference; dotted line marks \\(\\Delta P=-2\\).</figcaption>
          </figure>
          <figure>
            <img src="{item['diff_plot']}" alt="{cfg['title']} ECG grid difference">
            <figcaption><b>Difference.</b> Exact common-step differences for \\(\\Delta P\\) and \\(r_{{12}}\\).</figcaption>
          </figure>
          <figure>
            <img src="{item['health_plot']}" alt="{cfg['title']} ECG health metrics">
            <figcaption><b>ECG health.</b> Width floor and raw residual; the dashed line is the successful-case width floor.</figcaption>
          </figure>
          <figure>
            <img src="{item['delta_phase_plot']}" alt="{cfg['title']} displacement phase-axis comparison">
            <figcaption><b>Phase axis.</b> Same common-step displacement comparison plotted against \\(\\phi/\\pi\\).</figcaption>
          </figure>
        </div>
        <p class="source">ECG source: <code>{rel(cfg['ecg'])}</code><br>
        Grid source ({item['grid_source_kind']}): <code>{rel(item['grid_source'])}</code></p>
      </section>
    """


def main() -> None:
    items = [build_case(name, cfg) for name, cfg in CASES.items()]
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Vs3/Vl3 N=2 ECG vs Grid Comparison</title>
  <script>
  window.MathJax = {{
    tex: {{ inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']] }},
    svg: {{ fontCache: 'global' }}
  }};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
  <style>
    :root {{
      --line: #e2e8f0;
      --muted: #64748b;
      --teal: #0d9488;
      --blue: #2563eb;
      --ink: #0f172a;
    }}
    body {{
      font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
      max-width: 1240px;
      margin: 0 auto;
      padding: 24px 18px 62px;
      color: var(--ink);
      line-height: 1.55;
      background: #ffffff;
    }}
    h1 {{ font-size: 1.58rem; line-height: 1.32; margin: 0 0 6px; }}
    h2 {{
      font-size: 1.18rem;
      margin: 28px 0 10px;
      border-bottom: 2px solid var(--line);
      padding-bottom: 6px;
    }}
    h3 {{ font-size: 0.98rem; margin: 0 0 8px; }}
    p {{ margin: 9px 0; }}
    .sub, .meta, .source {{ color: var(--muted); }}
    .meta {{ margin-bottom: 14px; }}
    .banner {{
      border-radius: 8px;
      padding: 12px 16px;
      margin: 14px 0;
      border: 1px solid #bfdbfe;
      background: #eff6ff;
      border-left: 5px solid var(--blue);
    }}
    .ok {{
      background: #ecfdf5;
      border: 1px solid #a7f3d0;
      border-left: 5px solid #059669;
      border-radius: 8px;
      padding: 10px 14px;
      margin: 12px 0;
    }}
    .bad-banner {{
      background: #fef2f2;
      border: 1px solid #fecaca;
      border-left: 5px solid #dc2626;
      border-radius: 8px;
      padding: 10px 14px;
      margin: 12px 0;
    }}
    .warn {{
      background: #fffbeb;
      border: 1px solid #fde68a;
      border-left: 5px solid #d97706;
      border-radius: 8px;
      padding: 10px 14px;
      margin: 12px 0;
    }}
    .kpi {{
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      margin: 12px 0 4px;
    }}
    .kpi div {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 8px 13px;
      background: #f8fafc;
      min-width: 150px;
      flex: 1 1 180px;
      font-size: 0.9rem;
    }}
    .kpi b {{ display: block; font-size: 1.08rem; }}
    .equation {{
      margin: 12px 0 16px;
      padding: 14px 18px;
      border: 1px solid var(--line);
      border-left: 5px solid var(--teal);
      border-radius: 8px;
      background: #fbfcfd;
      font-size: 1.0rem;
      overflow-x: auto;
      white-space: normal;
      line-height: 1.45;
    }}
    .equation mjx-container[jax="SVG"][display="true"] {{ margin: 0.3em 0; }}
    table {{
      border-collapse: collapse;
      font-size: 0.86rem;
      margin: 9px 0;
      width: 100%;
    }}
    td, th {{
      border: 1px solid var(--line);
      padding: 5px 8px;
      text-align: center;
      vertical-align: top;
    }}
    th {{ background: #f8fafc; font-weight: 650; }}
    .num {{
      font-variant-numeric: tabular-nums;
      text-align: right;
      white-space: nowrap;
    }}
    .grid2 {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
      gap: 14px;
      align-items: start;
    }}
    .figgrid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px;
      align-items: start;
      margin-top: 14px;
    }}
    figure {{ text-align: center; margin: 0; }}
    img {{
      width: 100%;
      height: auto;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #ffffff;
    }}
    figcaption {{
      color: var(--muted);
      font-size: 0.82rem;
      margin-top: 6px;
      text-align: left;
      line-height: 1.38;
    }}
    code {{
      background: #f1f5f9;
      padding: 1px 5px;
      border-radius: 4px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 0.9em;
    }}
    .progpair {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
      gap: 12px;
      margin: 10px 0 14px;
    }}
    .progline {{
      display: grid;
      grid-template-columns: auto minmax(90px, 1fr) 50px;
      gap: 8px;
      align-items: center;
      color: var(--muted);
      font-size: 0.86rem;
    }}
    .progline b {{
      color: var(--ink);
      text-align: right;
      font-variant-numeric: tabular-nums;
    }}
    .prog {{
      height: 13px;
      background: #f1f5f9;
      border-radius: 8px;
      overflow: hidden;
      border: 1px solid var(--line);
    }}
    .prog > span {{ display: block; height: 100%; }}
    @media (max-width: 900px) {{
      .grid2, .figgrid, .progpair {{ grid-template-columns: 1fr; }}
      body {{ padding: 20px 13px 42px; }}
    }}
  </style>
</head>
<body>
<main>
  <h1>Vs3/Vl3 N=2 ECG vs Grid Comparison</h1>
  <p class="meta">Generated {generated} from existing progress files and grid CSV/log data in <code>{REPO}</code>.</p>
  <div class="banner"><b>Scope.</b> This is a compact live comparison of the two Slurm ECG jobs against the split-step grid reference. The grid reference is drawn as the grey dashed curve; ECG is the colored solid curve.</div>
  <div class="kpi">
    {kpi_cards(items)}
  </div>

  <section class="section">
    <h2>Definitions</h2>
    <p>This report compares the two running ECG-TDVP jobs with the currently available split-step grid-reference data. Grid curves are read from the finished grid <code>.csv</code> files when present, otherwise from Slurm progress lines.</p>
    <div class="equation">
      \\[
      P(t)=\\frac{{\\langle x_0+x_1\\rangle}}{{a}},
      \\qquad
      \\Delta P(t)=P(t)-P(0),
      \\qquad
      \\Delta P_{{\\mathrm{{ideal}}}}\\simeq -2 .
      \\]
    </div>
    <div class="banner">The comparison at common steps uses exact shared step numbers. When grid CSV files are available, this uses the same 25-step cadence as the ECG progress files; Slurm logs are used only as a fallback.</div>
  </section>

  <section class="section">
    <h2>Current Summary</h2>
    <table>
      <tr>
        <th>case</th>
        <th>ECG job</th>
        <th>ECG latest</th>
        <th>grid latest</th>
        <th>latest common step</th>
        <th>ECG \\(\\Delta P\\)</th>
        <th>grid \\(\\Delta P\\)</th>
        <th>difference</th>
        <th>max common |difference|</th>
        <th>ECG latest min width</th>
        <th>ECG latest residual</th>
      </tr>
      {summary_table(items)}
    </table>
  </section>

  {''.join(case_section(item) for item in items)}
</main>
</body>
</html>
"""
    OUT_HTML.write_text(html)
    print(OUT_HTML)


if __name__ == "__main__":
    main()
