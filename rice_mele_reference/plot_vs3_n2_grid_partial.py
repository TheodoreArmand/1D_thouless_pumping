#!/usr/bin/env python3
"""Plot partial Vs3/Vl3 N=2 grid-reference progress from Slurm logs."""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
REF_ROOT = REPO / "rice_mele_reference"
LOG_DIR = REPO / "slurm" / "vs3_n2_dt0p01_T160pi" / "grid_ref" / "ecg1d_vs3_n2_grid_ref_20260708"
TOTAL_STEPS = 50266
TOTAL_TIME = 502.6548245743669

CASES = {
    "free": {
        "log": LOG_DIR / "vs3n2grid_604164_0.out",
        "out_dir": REF_ROOT / "Vs3Vl3_3_3_2noin",
        "color": "#1f77b4",
        "label": "grid free",
    },
    "gauss": {
        "log": LOG_DIR / "vs3n2grid_604164_1.out",
        "out_dir": REF_ROOT / "Vs3Vl3_3_3_2gauss",
        "color": "#d62728",
        "label": "grid gauss",
    },
}

PROGRESS_RE = re.compile(
    r"progress grid N=2 step (?P<step>\d+)/(?P<total>\d+) "
    r"(?P<pct>[0-9.]+)% t=(?P<t>[-+0-9.eE]+) eta_s=(?P<eta>[-+0-9.eE]+) "
    r"P=(?P<P>[-+0-9.eE]+) dP=(?P<dP>[-+0-9.eE]+) "
    r"r12=(?P<r12>[-+0-9.eE]+) Vg=(?P<Vg>[-+0-9.eE]+) "
    r"norm=(?P<norm>[-+0-9.eE]+)"
)


def read_progress(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        m = PROGRESS_RE.search(line)
        if not m:
            continue
        row = {k: float(v) for k, v in m.groupdict().items()
               if k not in {"step", "total"}}
        row["step"] = int(m.group("step"))
        row["total"] = int(m.group("total"))
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no progress rows found in {path}")
    return rows


def write_csv(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["step", "total", "pct", "t", "eta", "P", "dP", "r12", "Vg", "norm"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in fields})


def plot_case(case_name: str, cfg, rows) -> None:
    out_dir = cfg["out_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    t_over_T = [r["t"] / TOTAL_TIME for r in rows]
    step = [r["step"] for r in rows]
    p = [r["P"] for r in rows]
    dp = [r["dP"] for r in rows]
    r12 = [r["r12"] for r in rows]
    vg = [r["Vg"] for r in rows]
    norm = [r["norm"] for r in rows]
    latest = rows[-1]

    fig, axes = plt.subplots(4, 1, figsize=(9.2, 9.5), sharex=True, constrained_layout=True)
    axes[0].plot(t_over_T, p, color=cfg["color"], lw=1.8)
    axes[0].set_ylabel("P = <x0+x1>/a")
    axes[0].grid(True, color="#d8d8d8", lw=0.7)

    axes[1].plot(t_over_T, dp, color=cfg["color"], lw=1.8)
    axes[1].axhline(-1.0, color="#444444", lw=0.9, ls="--")
    axes[1].set_ylabel("Delta P")
    axes[1].grid(True, color="#d8d8d8", lw=0.7)

    axes[2].plot(t_over_T, r12, color=cfg["color"], lw=1.8)
    axes[2].set_ylabel("r12 rms")
    axes[2].grid(True, color="#d8d8d8", lw=0.7)

    axes[3].plot(t_over_T, vg, color=cfg["color"], lw=1.8)
    axes[3].plot(t_over_T, norm, color="#555555", lw=1.0, alpha=0.7)
    axes[3].set_ylabel("Vg / norm")
    axes[3].set_xlabel("t / T")
    axes[3].grid(True, color="#d8d8d8", lw=0.7)

    axes[0].set_title(
        f"Partial grid reference ({case_name}): step {latest['step']}/{TOTAL_STEPS} "
        f"({100.0 * latest['step'] / TOTAL_STEPS:.1f}%), t/T={latest['t'] / TOTAL_TIME:.3f}"
    )
    fig.savefig(out_dir / "grid_reference_partial_from_log.png", dpi=170)
    plt.close(fig)

    write_csv(out_dir / "grid_reference_partial_from_log.csv", rows)


def plot_comparison(data) -> None:
    out_path = REF_ROOT / "Vs3Vl3_3_3_grid_reference_partial_compare.png"
    fig, axes = plt.subplots(3, 1, figsize=(9.6, 8.4), sharex=True, constrained_layout=True)
    for name, cfg, rows in data:
        t_over_T = [r["t"] / TOTAL_TIME for r in rows]
        axes[0].plot(t_over_T, [r["P"] for r in rows], label=cfg["label"], color=cfg["color"], lw=1.7)
        axes[1].plot(t_over_T, [r["dP"] for r in rows], label=cfg["label"], color=cfg["color"], lw=1.7)
        axes[2].plot(t_over_T, [r["r12"] for r in rows], label=cfg["label"], color=cfg["color"], lw=1.7)
        latest = rows[-1]
        axes[0].scatter([latest["t"] / TOTAL_TIME], [latest["P"]], color=cfg["color"], s=22)
    axes[1].axhline(-1.0, color="#444444", lw=0.9, ls="--")
    axes[0].set_ylabel("P")
    axes[1].set_ylabel("Delta P")
    axes[2].set_ylabel("r12 rms")
    axes[2].set_xlabel("t / T")
    for ax in axes:
        ax.grid(True, color="#d8d8d8", lw=0.7)
        ax.legend(frameon=False)
    axes[0].set_title("Partial Vs3/Vl3 N=2 grid reference from current Slurm logs")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def main() -> None:
    data = []
    for name, cfg in CASES.items():
        rows = read_progress(cfg["log"])
        plot_case(name, cfg, rows)
        data.append((name, cfg, rows))
        latest = rows[-1]
        print(
            f"{name}: step {latest['step']}/{latest['total']} "
            f"({latest['pct']:.1f}%), t/T={latest['t'] / TOTAL_TIME:.3f}, "
            f"P={latest['P']:.8f}, dP={latest['dP']:.8f}, "
            f"r12={latest['r12']:.8f}, Vg={latest['Vg']:.8g}"
        )
    plot_comparison(data)
    print(REF_ROOT / "Vs3Vl3_3_3_grid_reference_partial_compare.png")


if __name__ == "__main__":
    main()
