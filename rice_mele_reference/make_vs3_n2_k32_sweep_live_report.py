#!/usr/bin/env python3
"""Build a compact live HTML report for the Vs3/Vl3 N=2 K32 width sweep."""
from __future__ import annotations

import base64
import csv
import io
import argparse
import math
import os
import re
import time
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parent.parent
OUT_HTML = REPO / "vs3_n2_k32_sweep_live_report.html"
SCHEDULE_CSV = REPO / "rice_mele_reference" / "Vs3Vl3_3_3" / "gap_adaptive_vs3vl3_maincpp_schedule.csv"
REPORT_TITLE = "Vs3/Vl3 N=2 K32 live sweep"
REPORT_COMMAND = "python3 rice_mele_reference/make_vs3_n2_k32_sweep_live_report.py"
REPORT_NOTE = ""

TOTAL_TIME = 502.6548245743669
GRID_STEPS = 50266
IDEAL_DELTA_P = -2.0
SUCCESS_WIDTH_FLOOR = 0.0227284

ECG_ROOT = REPO / "out" / "vs3_n2_dt0p01_T160pi_K32_sweep"
GRID_ROOT = ECG_ROOT / "grid_ref"
GRID_LOG_DIR = (
    REPO
    / "slurm"
    / "vs3_n2_dt0p01_T160pi_K32_sweep"
    / "grid_ref"
    / "ecg1d_vs3_n2_grid_k32_sweep_20260709"
)

OLD_GRID_FREE = (
    REPO
    / "out"
    / "vs3_n2_dt0p01_T160pi"
    / "grid_ref"
    / "free"
    / "g1024_dt0p01_T160pi"
    / "n2_grid_reference_vs3_free_g1024_dt0p01_T160pi.csv"
)
OLD_GRID_GAUSS_SIGMA1 = (
    REPO
    / "out"
    / "vs3_n2_dt0p01_T160pi"
    / "grid_ref"
    / "gauss"
    / "g1024_dt0p01_T160pi"
    / "n2_grid_reference_vs3_gauss_g1024_dt0p01_T160pi.csv"
)

GRID_PROGRESS_RE = re.compile(
    r"progress grid N=2 step (?P<step>\d+)/(?P<total>\d+) "
    r"(?P<pct>[0-9.]+)% t=(?P<t>[-+0-9.eE]+) eta_s=(?P<eta>[-+0-9.eE]+) "
    r"P=(?P<P>[-+0-9.eE]+) dP=(?P<dP>[-+0-9.eE]+) "
    r"r12=(?P<r12>[-+0-9.eE]+) Vg=(?P<Vg>[-+0-9.eE]+) "
    r"norm=(?P<norm>[-+0-9.eE]+)"
)

_SCHEDULE: tuple[list[float], list[float]] | None = None


def load_schedule() -> tuple[list[float], list[float]]:
    global _SCHEDULE
    if _SCHEDULE is not None:
        return _SCHEDULE
    ss: list[float] = []
    phis: list[float] = []
    if SCHEDULE_CSV.exists():
        with SCHEDULE_CSV.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                s = parse_float(row.get("s"))
                phi = parse_float(row.get("phi"))
                if math.isfinite(s) and math.isfinite(phi):
                    ss.append(s)
                    phis.append(phi)
    _SCHEDULE = (ss, phis)
    return _SCHEDULE


def phi_at_time(t: float) -> float:
    if not math.isfinite(t):
        return float("nan")
    ss, phis = load_schedule()
    if not ss:
        return 2.0 * math.pi * max(0.0, min(1.0, t / TOTAL_TIME))
    s = max(0.0, min(1.0, t / TOTAL_TIME))
    if s <= ss[0]:
        return phis[0]
    if s >= ss[-1]:
        return phis[-1]
    lo = 0
    hi = len(ss) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ss[mid] <= s:
            lo = mid
        else:
            hi = mid
    span = ss[hi] - ss[lo]
    if span <= 0:
        return phis[lo]
    alpha = (s - ss[lo]) / span
    return phis[lo] * (1.0 - alpha) + phis[hi] * alpha


def tag(x: float) -> str:
    s = f"{x:.3f}"
    return s.replace("-", "m").replace(".", "p")


CASES = [
    {
        "key": "free",
        "title": "free",
        "mode": "free",
        "sigma": None,
        "grid_index": 0,
        "color": "#2267a8",
    },
    *[
        {
            "key": f"gauss_sigma{tag(sigma)}",
            "title": f"gauss sigma={sigma:g}",
            "mode": "gauss",
            "sigma": sigma,
            "grid_index": idx,
            "color": color,
        }
        for idx, (sigma, color) in enumerate(
            [
                (0.5, "#b03a2e"),
                (0.75, "#a855f7"),
                (1.0, "#dc2626"),
                (1.25, "#ea580c"),
                (1.5, "#ca8a04"),
                (2.0, "#059669"),
                (3.0, "#0891b2"),
            ],
            start=1,
        )
    ],
]


def rel(path: Path | None) -> str:
    if path is None:
        return "n/a"
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def parse_float(value: str | None, default: float = float("nan")) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def parse_int(value: str | None, default: int = -1) -> int:
    x = parse_float(value)
    if not math.isfinite(x):
        return default
    return int(round(x))


def fnum(x: float | None, digits: int = 4) -> str:
    if x is None or not math.isfinite(x):
        return "n/a"
    if x == 0:
        return "0"
    if abs(x) < 1e-4 or abs(x) >= 1e5:
        return f"{x:.{digits}e}"
    return f"{x:.{digits}f}"


def fmt_eta(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return "n/a"
    hours = seconds / 3600.0
    if hours < 1:
        return f"{seconds / 60.0:.0f} min"
    return f"{hours:.1f} h"


def ecg_progress_path(case: dict) -> Path:
    if case["mode"] == "free":
        root = ECG_ROOT / "free"
    else:
        root = ECG_ROOT / case["key"]
    return root / "a8p000_K32_tmax502p655_VsER3p000_VlER3p000" / "progress.csv"


def grid_csv_path(case: dict) -> Path:
    if case["mode"] == "free":
        subdir = "free"
        name = "n2_grid_reference_vs3_free_g1024_dt0p01_T160pi.csv"
    else:
        subdir = case["key"]
        name = f"n2_grid_reference_vs3_{case['key']}_g1024_dt0p01_T160pi.csv"
    return GRID_ROOT / subdir / "g1024_dt0p01_T160pi" / name


def old_grid_fallback(case: dict) -> Path | None:
    if case["mode"] == "free":
        return OLD_GRID_FREE if OLD_GRID_FREE.exists() else None
    if case["sigma"] == 1.0:
        return OLD_GRID_GAUSS_SIGMA1 if OLD_GRID_GAUSS_SIGMA1.exists() else None
    return None


def grid_log_candidates(case: dict) -> list[Path]:
    pattern = f"vs3n2gk32_*_{case['grid_index']}.out"
    return sorted(GRID_LOG_DIR.glob(pattern), key=lambda p: p.stat().st_mtime)


def read_ecg_csv(path: Path) -> list[dict]:
    rows: list[dict] = []
    if not path.exists():
        return rows
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                step = parse_int(row.get("step"))
                if step < 0:
                    continue
                total = parse_int(row.get("steps_total_est"), 0)
                t = parse_float(row.get("t"))
                phi = parse_float(row.get("phi"))
                rows.append(
                    {
                        "step": step,
                        "total": total,
                        "progress": 100.0 * step / total if total > 0 else float("nan"),
                        "t": t,
                        "phi": phi,
                        "phase": phi / math.pi if math.isfinite(phi) else float("nan"),
                        "P": parse_float(row.get("polarization_cell")),
                        "dP": parse_float(row.get("delta_polarization")),
                        "r12": parse_float(row.get("r12_rms")),
                        "Vg": parse_float(row.get("V_gauss")),
                        "norm": parse_float(row.get("norm")),
                        "eta": parse_float(row.get("eta_seconds")),
                        "resid": parse_float(row.get("relative_raw_residual")),
                        "discarded": parse_float(row.get("discarded_rhs_fraction")),
                        "minAB": parse_float(row.get("min_re_AplusB")),
                        "minB": parse_float(row.get("min_re_B")),
                        "rank": parse_int(row.get("actual_solve_rank")),
                        "sv_max": parse_float(row.get("sv_max")),
                    }
                )
            except (KeyError, ValueError):
                continue
    return rows


def read_grid_csv(path: Path, source_kind: str) -> tuple[list[dict], Path | None, str]:
    rows: list[dict] = []
    if not path.exists():
        return rows, None, source_kind
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t = parse_float(row.get("t"))
                if not math.isfinite(t):
                    continue
                step = GRID_STEPS if abs(t - TOTAL_TIME) < 1e-7 else int(round(t / 0.01))
                phi = parse_float(row.get("phi"))
                rows.append(
                    {
                        "step": step,
                        "total": GRID_STEPS,
                        "progress": 100.0 * step / GRID_STEPS,
                        "t": t,
                        "phi": phi,
                        "phase": phi / math.pi if math.isfinite(phi) else float("nan"),
                        "P": parse_float(row.get("polarization_cell")),
                        "dP": parse_float(row.get("delta_polarization")),
                        "r12": parse_float(row.get("r12_rms")),
                        "Vg": parse_float(row.get("V_gauss")),
                        "norm": parse_float(row.get("norm")),
                        "eta": 0.0 if step >= GRID_STEPS else float("nan"),
                    }
                )
            except (KeyError, ValueError):
                continue
    return rows, path, source_kind


def read_grid_log(path: Path) -> tuple[list[dict], Path | None, str]:
    rows: list[dict] = []
    if not path.exists():
        return rows, None, "grid log"
    for line in path.read_text(errors="ignore").splitlines():
        m = GRID_PROGRESS_RE.search(line)
        if not m:
            continue
        t = float(m.group("t"))
        phi = phi_at_time(t)
        rows.append(
            {
                "step": int(m.group("step")),
                "total": int(m.group("total")),
                "progress": float(m.group("pct")),
                "t": t,
                "phi": phi,
                "phase": phi / math.pi if math.isfinite(phi) else float("nan"),
                "P": float(m.group("P")),
                "dP": float(m.group("dP")),
                "r12": float(m.group("r12")),
                "Vg": float(m.group("Vg")),
                "norm": float(m.group("norm")),
                "eta": float(m.group("eta")),
            }
        )
    return rows, path, "grid slurm log"


def read_grid(case: dict) -> tuple[list[dict], Path | None, str]:
    new_csv = grid_csv_path(case)
    if new_csv.exists():
        return read_grid_csv(new_csv, "K32 grid CSV")

    logs = grid_log_candidates(case)
    if logs:
        rows, path, kind = read_grid_log(logs[-1])
        if rows:
            return rows, path, kind

    old = old_grid_fallback(case)
    if old is not None:
        return read_grid_csv(old, "old K24 grid CSV fallback")

    return [], None, "missing"


def interpolate_grid_by_phase(grid: list[dict], phase: float, key: str) -> float:
    if not math.isfinite(phase):
        return float("nan")
    rows = [
        r for r in grid
        if math.isfinite(r.get("phase", float("nan"))) and math.isfinite(r.get(key, float("nan")))
    ]
    if not rows:
        return float("nan")
    rows.sort(key=lambda r: r["phase"])
    if phase < rows[0]["phase"] or phase > rows[-1]["phase"]:
        return float("nan")
    if phase == rows[0]["phase"]:
        return rows[0][key]
    for lo, hi in zip(rows, rows[1:]):
        plo = lo["phase"]
        phi = hi["phase"]
        if phase <= phi:
            if phi <= plo:
                return lo[key]
            alpha = (phase - plo) / (phi - plo)
            return lo[key] * (1.0 - alpha) + hi[key] * alpha
    return rows[-1][key]


def matched_phase_rows(ecg: list[dict], grid: list[dict]) -> list[dict]:
    matched: list[dict] = []
    for e in ecg:
        phase = e["phase"]
        if not math.isfinite(phase):
            continue
        grid_dP = interpolate_grid_by_phase(grid, phase, "dP")
        grid_r12 = interpolate_grid_by_phase(grid, phase, "r12")
        if not (math.isfinite(grid_dP) and math.isfinite(grid_r12)):
            continue
        matched.append(
            {
                "step": e["step"],
                "phase": phase,
                "ecg_dP": e["dP"],
                "grid_dP": grid_dP,
                "diff_dP": e["dP"] - grid_dP,
                "ecg_r12": e["r12"],
                "grid_r12": grid_r12,
                "diff_r12": e["r12"] - grid_r12,
            }
        )
    return matched


def fig_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def style_axis(ax) -> None:
    ax.grid(True, color="#e2e8f0", lw=0.55)
    ax.set_axisbelow(True)
    ax.tick_params(colors="#334155", labelsize=7.0, width=0.65)
    for spine in ax.spines.values():
        spine.set_color("#cbd5e1")
        spine.set_linewidth(0.7)


def phase_values(rows: list[dict]) -> list[float]:
    values: list[float] = []
    for r in rows:
        phase = r["phase"]
        if not math.isfinite(phase) and math.isfinite(r["t"]):
            phi = phi_at_time(r["t"])
            phase = phi / math.pi if math.isfinite(phi) else float("nan")
        values.append(phase)
    return values


def plot_transport(case: dict, ecg: list[dict], grid: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(3.18, 2.05))
    if grid:
        ax.plot(phase_values(grid), [r["dP"] for r in grid], color="#475569", lw=1.15, ls="--", label="grid")
    if ecg:
        ax.plot(phase_values(ecg), [r["dP"] for r in ecg], color=case["color"], lw=1.25, label="ECG")
    ax.axhline(IDEAL_DELTA_P, color="#94a3b8", lw=0.75, ls=":")
    ax.set_xlim(0.0, 2.0)
    ax.set_xlabel("phi / pi", fontsize=7.5)
    ax.set_ylabel("Delta P", fontsize=7.5)
    ax.set_title("transport", fontsize=8.4)
    style_axis(ax)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=False, fontsize=6.8, loc="best", handlelength=1.7)
    fig.tight_layout()
    return fig_uri(fig)


def plot_difference(case: dict, common: list[dict]) -> str:
    fig, axes = plt.subplots(2, 1, figsize=(3.18, 2.25), sharex=True)
    x = [r["phase"] for r in common if math.isfinite(r["phase"])]
    rows = [r for r in common if math.isfinite(r["phase"])]
    if rows:
        axes[0].plot(x, [r["diff_dP"] for r in rows], color=case["color"], lw=1.1)
        axes[1].plot(x, [r["diff_r12"] for r in rows], color="#0f766e", lw=1.1)
    for ax in axes:
        ax.axhline(0.0, color="#94a3b8", lw=0.65)
        ax.set_xlim(0.0, 2.0)
        style_axis(ax)
    axes[0].set_ylabel("dDP", fontsize=7.2)
    axes[1].set_ylabel("dr12", fontsize=7.2)
    axes[1].set_xlabel("phi / pi", fontsize=7.5)
    axes[0].set_title("ECG - grid", fontsize=8.4)
    fig.tight_layout(h_pad=0.2)
    return fig_uri(fig)


def plot_health(case: dict, ecg: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(3.18, 2.05))
    x = phase_values(ecg)
    min_ab = [r["minAB"] if r["minAB"] > 0 else float("nan") for r in ecg]
    resid = [r["resid"] if r["resid"] > 0 else float("nan") for r in ecg]
    ln1 = ax.plot(x, min_ab, color="#334155", lw=1.05, label="lambda_min Re(A+B)")
    ax.axhline(SUCCESS_WIDTH_FLOOR, color="#d97706", lw=0.75, ls="--")
    ax.axhline(1e-3, color="#dc2626", lw=0.65, ls=":")
    ax.set_yscale("log")
    ax.set_xlim(0.0, 2.0)
    ax.set_xlabel("phi / pi", fontsize=7.5)
    ax.set_ylabel("width", fontsize=7.2)
    ax.set_title("health", fontsize=8.4)
    style_axis(ax)

    ax2 = ax.twinx()
    ln2 = ax2.plot(x, resid, color=case["color"], lw=0.95, alpha=0.9, label="res")
    ax2.set_yscale("log")
    ax2.set_ylabel("res", fontsize=7.2, color=case["color"])
    ax2.tick_params(colors=case["color"], labelsize=7.0, width=0.65)
    for spine in ax2.spines.values():
        spine.set_color("#cbd5e1")
        spine.set_linewidth(0.7)
    lines = ln1 + ln2
    ax.legend(lines, [line.get_label() for line in lines], frameon=False, fontsize=6.3, loc="best")
    fig.tight_layout()
    return fig_uri(fig)


def latest(rows: list[dict]) -> dict | None:
    return rows[-1] if rows else None


def max_abs(rows: list[dict], key: str) -> float:
    if not rows:
        return float("nan")
    vals = [abs(r[key]) for r in rows if math.isfinite(r[key])]
    return max(vals) if vals else float("nan")


def build_case(case: dict) -> dict:
    ecg_path = ecg_progress_path(case)
    ecg_rows = read_ecg_csv(ecg_path)
    grid_rows, grid_source, grid_kind = read_grid(case)
    common = matched_phase_rows(ecg_rows, grid_rows)
    le = latest(ecg_rows)
    lg = latest(grid_rows)
    return {
        **case,
        "ecg_path": ecg_path,
        "grid_source": grid_source,
        "grid_kind": grid_kind,
        "ecg": ecg_rows,
        "grid": grid_rows,
        "common": common,
        "latest_ecg": le,
        "latest_grid": lg,
        "transport_plot": plot_transport(case, ecg_rows, grid_rows),
        "diff_plot": plot_difference(case, common),
        "health_plot": plot_health(case, ecg_rows),
        "max_diff_dP": max_abs(common, "diff_dP"),
        "max_diff_r12": max_abs(common, "diff_r12"),
    }


def progress_bar(label: str, pct: float, eta: float | None, color: str) -> str:
    value = max(0.0, min(100.0, pct)) if math.isfinite(pct) else 0.0
    return f"""
      <div class="progline">
        <span>{label}</span>
        <div class="prog"><i style="width:{value:.2f}%;background:{color}"></i></div>
        <b>{value:.1f}%</b>
        <em>{fmt_eta(eta)}</em>
      </div>
    """


def case_section(item: dict) -> str:
    le = item["latest_ecg"]
    lg = item["latest_grid"]
    ecg_pct = le["progress"] if le else float("nan")
    grid_pct = lg["progress"] if lg else float("nan")
    ecg_eta = le["eta"] if le else float("nan")
    grid_eta = lg["eta"] if lg else float("nan")
    common = item["common"]
    lc = common[-1] if common else None
    minab = le["minAB"] if le else float("nan")
    resid = le["resid"] if le else float("nan")
    status_class = "ok"
    if math.isfinite(minab) and minab < 1e-3:
        status_class = "bad"
    elif math.isfinite(minab) and minab < SUCCESS_WIDTH_FLOOR:
        status_class = "warn"

    if lc:
        compare = (
            f"phase {fnum(lc['phase'], 3)} pi: dDP={fnum(lc['diff_dP'])}, "
            f"dr12={fnum(lc['diff_r12'])}"
        )
    else:
        compare = "phase-matched comparison: waiting for grid data"

    return f"""
    <section class="case">
      <div class="case-head">
        <div>
          <h2>{item['title']}</h2>
          <p>{compare}</p>
        </div>
        <div class="{status_class}">
          \\( \\lambda_{{\\min}}[\\operatorname{{Re}}(A+B)] \\)={fnum(minab)}<br>
          res={fnum(resid)}
        </div>
      </div>
      <div class="progress-grid">
        {progress_bar("ECG", ecg_pct, ecg_eta, item["color"])}
        {progress_bar("grid", grid_pct, grid_eta, "#475569")}
      </div>
      <div class="mini">
        <span>ECG step {le['step'] if le else 'n/a'}/{le['total'] if le else 'n/a'}</span>
        <span>grid source: {item['grid_kind']}</span>
        <span>max |dDP|={fnum(item['max_diff_dP'])}</span>
        <span>max |dr12|={fnum(item['max_diff_r12'])}</span>
      </div>
      <div class="figs">
        <figure><img src="{item['transport_plot']}" alt="{item['title']} transport"><figcaption>Delta P vs phase; grid is dashed.</figcaption></figure>
        <figure><img src="{item['diff_plot']}" alt="{item['title']} difference"><figcaption>ECG minus phase-interpolated grid.</figcaption></figure>
        <figure><img src="{item['health_plot']}" alt="{item['title']} health"><figcaption>Width floor and residual monitor.</figcaption></figure>
      </div>
      <p class="source">ECG: <code>{rel(item['ecg_path'])}</code><br>grid: <code>{rel(item['grid_source'])}</code></p>
    </section>
    """


def summary_rows(items: list[dict]) -> str:
    rows = []
    for item in items:
        le = item["latest_ecg"]
        lg = item["latest_grid"]
        lc = item["common"][-1] if item["common"] else None
        rows.append(
            f"""
            <tr>
              <td>{item['title']}</td>
              <td>{fnum(le['progress'] if le else float('nan'), 2)}%</td>
              <td>{fmt_eta(le['eta'] if le else float('nan'))}</td>
              <td>{fnum(lg['progress'] if lg else float('nan'), 2)}%</td>
              <td>{fmt_eta(lg['eta'] if lg else float('nan'))}</td>
              <td>{fnum(lc['phase'] if lc else float('nan'), 3)}</td>
              <td class="num">{fnum(lc['diff_dP'] if lc else float('nan'))}</td>
              <td class="num">{fnum(le['minAB'] if le else float('nan'))}</td>
              <td class="num">{fnum(le['resid'] if le else float('nan'))}</td>
            </tr>
            """
        )
    return "\n".join(rows)


def write_report() -> Path:
    items = [build_case(case) for case in CASES]
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{REPORT_TITLE}</title>
  <script>
  window.MathJax = {{
    tex: {{ inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']] }},
    svg: {{ fontCache: 'global' }}
  }};
  </script>
  <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
  <style>
    body {{
      font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
      max-width: 1320px;
      margin: 0 auto;
      padding: 18px 16px 52px;
      color: #0f172a;
      background: #fff;
    }}
    h1 {{ font-size: 1.45rem; margin: 0 0 4px; }}
    h2 {{ font-size: 1.02rem; margin: 0 0 2px; }}
    p {{ margin: 5px 0; }}
    .meta, .source, figcaption, .mini {{ color: #64748b; }}
    .banner {{
      border: 1px solid #bfdbfe;
      border-left: 5px solid #2563eb;
      background: #eff6ff;
      border-radius: 8px;
      padding: 9px 12px;
      margin: 10px 0 14px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
      margin: 10px 0 18px;
    }}
    th, td {{
      border: 1px solid #e2e8f0;
      padding: 5px 7px;
      text-align: center;
    }}
    th {{ background: #f8fafc; }}
    .num {{ text-align: right; font-variant-numeric: tabular-nums; }}
    .case {{
      border-top: 2px solid #e2e8f0;
      padding-top: 11px;
      margin-top: 15px;
    }}
    .case-head {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 12px;
      align-items: start;
    }}
    .ok, .warn, .bad {{
      border-radius: 8px;
      padding: 7px 10px;
      font-size: 0.83rem;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }}
    .ok {{ background: #ecfdf5; border: 1px solid #a7f3d0; }}
    .warn {{ background: #fffbeb; border: 1px solid #fde68a; }}
    .bad {{ background: #fef2f2; border: 1px solid #fecaca; }}
    .progress-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 9px;
      margin: 7px 0;
    }}
    .progline {{
      display: grid;
      grid-template-columns: 38px minmax(90px, 1fr) 52px 58px;
      gap: 7px;
      align-items: center;
      font-size: 0.8rem;
      color: #64748b;
    }}
    .progline b, .progline em {{
      color: #0f172a;
      font-style: normal;
      text-align: right;
      font-variant-numeric: tabular-nums;
    }}
    .prog {{
      height: 11px;
      border-radius: 999px;
      background: #f1f5f9;
      border: 1px solid #e2e8f0;
      overflow: hidden;
    }}
    .prog i {{ display: block; height: 100%; }}
    .mini {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
      font-size: 0.8rem;
      margin-bottom: 7px;
    }}
    .figs {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
    }}
    figure {{ margin: 0; }}
    img {{
      width: 100%;
      height: auto;
      border: 1px solid #e2e8f0;
      border-radius: 7px;
      background: #fff;
    }}
    figcaption {{ font-size: 0.76rem; line-height: 1.3; margin-top: 4px; }}
    code {{
      background: #f1f5f9;
      padding: 1px 4px;
      border-radius: 4px;
      font-size: 0.9em;
    }}
    .source {{ font-size: 0.76rem; line-height: 1.35; }}
    @media (max-width: 900px) {{
      .figs, .progress-grid, .case-head {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <h1>{REPORT_TITLE}</h1>
  <p class="meta">Generated {generated}. Rerun <code>{REPORT_COMMAND}</code> to refresh from current files.</p>
  <div class="banner">
    Eight cases are shown below. Each case has three compact phase-axis plots: transport, ECG-grid difference, and ECG health.
    New grid references use <code>initial_pathpad_N2_K32.csv</code>; old K24 grid CSV is used only as a temporary fallback until the new grid sweep finishes.
    {REPORT_NOTE}
  </div>
  <table>
    <tr>
      <th>case</th>
      <th>ECG progress</th>
      <th>ECG ETA</th>
      <th>grid progress</th>
      <th>grid ETA</th>
      <th>matched phase/pi</th>
      <th>ECG-grid Delta P</th>
      <th>\\(\\lambda_{{\\min}}[\\operatorname{{Re}}(A+B)]\\)</th>
      <th>res</th>
    </tr>
    {summary_rows(items)}
  </table>
  {''.join(case_section(item) for item in items)}
</body>
</html>
"""
    OUT_HTML.write_text(html)
    return OUT_HTML


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--watch",
        type=float,
        default=0.0,
        help="Regenerate the report every WATCH seconds until interrupted.",
    )
    args = parser.parse_args()

    while True:
        path = write_report()
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {path}", flush=True)
        if args.watch <= 0:
            break
        time.sleep(args.watch)


if __name__ == "__main__":
    main()
