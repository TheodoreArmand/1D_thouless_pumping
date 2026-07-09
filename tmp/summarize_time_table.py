#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_CSV = ROOT / "tmp" / "time_summary_all.csv"
OUT_MD = ROOT / "tmp" / "time_summary.md"


def parse_float(value):
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def parse_int(value):
    x = parse_float(value)
    if not math.isfinite(x):
        return None
    return int(round(x))


def clean_value(value):
    return value.split("#", 1)[0].strip()


def parse_key_value_file(path):
    data = {}
    if not path.exists():
        return data
    for line in path.read_text(errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = clean_value(value)
    return data


def read_last_progress(path):
    if not path.exists():
        return {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        last = None
        for row in reader:
            last = row
    return last or {}


def short_run_name(run_dir):
    rel = run_dir.relative_to(ROOT)
    parts = rel.parts
    if len(parts) >= 3 and parts[0] == "out":
        return "/".join(parts[1:-1])
    if len(parts) >= 3 and parts[0] == "tmp":
        return "/".join(parts[:-1])
    return str(rel)


def infer_family(run_dir):
    rel = run_dir.relative_to(ROOT)
    parts = rel.parts
    if parts and parts[0] == "out" and len(parts) >= 2:
        return parts[1]
    if parts and parts[0] == "tmp":
        return "/".join(parts[:2])
    return parts[0] if parts else "."


def collect_runs():
    dirs = set()
    for base in [ROOT / "out", ROOT / "tmp"]:
        if not base.exists():
            continue
        for path in base.rglob("summary.txt"):
            dirs.add(path.parent)
        for path in base.rglob("progress.csv"):
            dirs.add(path.parent)

    rows = []
    for run_dir in sorted(dirs):
        summary = parse_key_value_file(run_dir / "summary.txt")
        config = parse_key_value_file(run_dir / "config.txt")
        progress = read_last_progress(run_dir / "progress.csv")

        merged = {}
        merged.update(config)
        merged.update(summary)

        summary_done = bool(summary)
        steps_done = parse_int(summary.get("steps_total")) if summary_done else None
        if steps_done is None:
            steps_done = parse_int(progress.get("step"))

        steps_est = parse_int(summary.get("estimated_time_steps"))
        if steps_est is None:
            steps_est = parse_int(progress.get("steps_total_est"))

        if summary_done:
            wall_s = parse_float(summary.get("evolution_wall_seconds"))
            sec_per_step = parse_float(summary.get("evolution_seconds_per_step"))
            eta_s = 0.0 if math.isfinite(wall_s) else math.nan
        else:
            wall_s = parse_float(progress.get("wall_seconds"))
            sec_per_step = parse_float(progress.get("seconds_per_step"))
            eta_s = parse_float(progress.get("eta_seconds"))

        progress_frac = parse_float(progress.get("progress_frac"))
        if not math.isfinite(progress_frac) and steps_done is not None and steps_est:
            progress_frac = steps_done / steps_est
        if summary_done:
            progress_frac = 1.0

        accepted_dt = parse_float(progress.get("accepted_dt"))
        min_dt = parse_float(summary.get("min_accepted_dt"))
        max_dt = parse_float(summary.get("max_accepted_dt"))
        if not math.isfinite(min_dt):
            min_dt = accepted_dt
        if not math.isfinite(max_dt):
            max_dt = accepted_dt

        row = {
            "status": "done" if summary_done else "progress",
            "family": infer_family(run_dir),
            "run": short_run_name(run_dir),
            "config_name": merged.get("config_name", ""),
            "model": merged.get("model", ""),
            "potential": merged.get("potential", ""),
            "N": merged.get("N", ""),
            "K": merged.get("K", ""),
            "param_dim": merged.get("param_dim", progress.get("param_dim", "")),
            "time_step_mode": merged.get("time_step_mode", ""),
            "dt": merged.get("dt", ""),
            "steps_done": "" if steps_done is None else steps_done,
            "steps_est": "" if steps_est is None else steps_est,
            "progress_pct": progress_frac * 100.0 if math.isfinite(progress_frac) else math.nan,
            "wall_hours": wall_s / 3600.0 if math.isfinite(wall_s) else math.nan,
            "sec_per_step": sec_per_step,
            "eta_hours": eta_s / 3600.0 if math.isfinite(eta_s) else math.nan,
            "min_dt": min_dt,
            "max_dt": max_dt,
            "path": str(run_dir.relative_to(ROOT)),
        }
        rows.append(row)
    return rows


def fmt_num(value, digits=3):
    if value is None:
        return ""
    if isinstance(value, str):
        value = parse_float(value)
    if not math.isfinite(value):
        return ""
    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.{digits}f}"


def fmt_cell(value):
    if isinstance(value, float):
        return fmt_num(value)
    return str(value)


def markdown_table(rows, columns):
    if not rows:
        return "_No rows._\n"
    out = []
    out.append("| " + " | ".join(title for _, title in columns) + " |")
    out.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        out.append("| " + " | ".join(fmt_cell(row.get(key, "")) for key, _ in columns) + " |")
    return "\n".join(out) + "\n"


def write_csv(rows):
    fields = [
        "status", "family", "run", "config_name", "model", "potential",
        "N", "K", "param_dim", "time_step_mode", "dt",
        "steps_done", "steps_est", "progress_pct", "wall_hours",
        "sec_per_step", "eta_hours", "min_dt", "max_dt", "path",
    ]
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def aggregate_families(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["family"]].append(row)

    out = []
    for family, items in sorted(groups.items()):
        progress_items = [r for r in items if r["status"] == "progress"]
        done_items = [r for r in items if r["status"] == "done"]
        sps_values = [r["sec_per_step"] for r in items if math.isfinite(r["sec_per_step"])]
        wall_values = [r["wall_hours"] for r in items if math.isfinite(r["wall_hours"])]
        eta_values = [r["eta_hours"] for r in progress_items if math.isfinite(r["eta_hours"])]
        progress_values = [r["progress_pct"] for r in items if math.isfinite(r["progress_pct"])]
        out.append({
            "family": family,
            "runs": len(items),
            "done": len(done_items),
            "progressing": len(progress_items),
            "mean_progress_pct": sum(progress_values) / len(progress_values) if progress_values else math.nan,
            "mean_sec_per_step": sum(sps_values) / len(sps_values) if sps_values else math.nan,
            "sum_wall_hours": sum(wall_values) if wall_values else math.nan,
            "max_eta_hours": max(eta_values) if eta_values else 0.0 if progress_items else "",
        })
    return out


def write_markdown(rows):
    bench = [
        r for r in rows
        if r["status"] == "done" and (
            r["family"].startswith("bench_") or r["family"].startswith("dtscan_")
        )
    ]
    bench.sort(key=lambda r: (r["family"], r["run"]))

    k32_progress = [
        r for r in rows
        if r["status"] == "progress"
        and "vs3_n2_dt0p01_T160pi_K32" in r["family"]
        and "/grid_ref/" not in r["path"]
    ]
    k32_progress.sort(key=lambda r: (r["family"], r["run"]))

    families = aggregate_families([
        r for r in rows
        if "vs3_n2_dt0p01_T160pi_K32" in r["family"]
        or r["family"].startswith("bench_")
        or r["family"].startswith("dtscan_")
    ])

    bench_cols = [
        ("run", "run"),
        ("steps_done", "steps"),
        ("param_dim", "dim"),
        ("dt", "dt"),
        ("wall_hours", "wall h"),
        ("sec_per_step", "s/step"),
    ]
    progress_cols = [
        ("run", "run"),
        ("steps_done", "step"),
        ("steps_est", "est"),
        ("progress_pct", "%"),
        ("wall_hours", "elapsed h"),
        ("sec_per_step", "s/step"),
        ("eta_hours", "ETA h"),
        ("min_dt", "dt"),
    ]
    family_cols = [
        ("family", "family"),
        ("runs", "runs"),
        ("done", "done"),
        ("progressing", "progress"),
        ("mean_progress_pct", "mean %"),
        ("mean_sec_per_step", "mean s/step"),
        ("sum_wall_hours", "sum wall h"),
        ("max_eta_hours", "max ETA h"),
    ]

    text = []
    text.append("# Time Summary\n")
    text.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")
    text.append(f"Rows scanned: {len(rows)}\n")
    text.append("Full per-run data: `tmp/time_summary_all.csv`.\n")
    text.append("## Family Aggregate\n")
    text.append(markdown_table(families, family_cols))
    text.append("## Completed Benchmarks\n")
    text.append(markdown_table(bench, bench_cols))
    text.append("## K32 Long-Run Progress\n")
    text.append(markdown_table(k32_progress, progress_cols))
    OUT_MD.write_text("\n".join(text))


def main():
    rows = collect_runs()
    write_csv(rows)
    write_markdown(rows)
    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"runs={len(rows)}")


if __name__ == "__main__":
    main()
