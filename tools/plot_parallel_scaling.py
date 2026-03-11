#!/usr/bin/env python3
"""Plot OpenMP parallel scaling from two run_summary.csv files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt


def read_summary(path: Path) -> List[Tuple[int, float]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")

    rows: List[Tuple[int, float]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        required = {"threads", "elapsed_seconds", "exit_code"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"Invalid header in {path}, need fields: {sorted(required)}")

        for row in reader:
            try:
                exit_code = int(row["exit_code"])
                if exit_code != 0:
                    continue
                threads = int(row["threads"])
                elapsed = float(row["elapsed_seconds"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Bad row in {path}: {row}") from exc
            rows.append((threads, elapsed))

    if not rows:
        raise ValueError(f"No successful rows (exit_code==0) found in {path}")

    rows.sort(key=lambda x: x[0])
    return rows


def compute_speedup(rows: List[Tuple[int, float]]) -> Tuple[List[int], List[float], List[float]]:
    threads = [r[0] for r in rows]
    elapsed = [r[1] for r in rows]
    t1 = elapsed[0]
    speedup = [t1 / t for t in elapsed]
    efficiency = [s / n for s, n in zip(speedup, threads)]
    return threads, speedup, efficiency


def plot_scaling(
    csv_100nm: Path,
    csv_1000nm: Path,
    output_png: Path,
    label_100nm: str,
    label_1000nm: str,
) -> None:
    data_100nm = read_summary(csv_100nm)
    data_1000nm = read_summary(csv_1000nm)

    th_100nm = [r[0] for r in data_100nm]
    tm_100nm = [r[1] for r in data_100nm]
    th_1000nm = [r[0] for r in data_1000nm]
    tm_1000nm = [r[1] for r in data_1000nm]

    sp_th_100nm, sp_100nm, ef_100nm = compute_speedup(data_100nm)
    sp_th_1000nm, sp_1000nm, ef_1000nm = compute_speedup(data_1000nm)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    axes[0].plot(th_100nm, tm_100nm, marker="o", linewidth=2, label=label_100nm)
    axes[0].plot(th_1000nm, tm_1000nm, marker="s", linewidth=2, label=label_1000nm)
    axes[0].set_title("Runtime vs Threads")
    axes[0].set_xlabel("Threads")
    axes[0].set_ylabel("Elapsed time (s)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(sp_th_100nm, sp_100nm, marker="o", linewidth=2, label=label_100nm)
    axes[1].plot(sp_th_1000nm, sp_1000nm, marker="s", linewidth=2, label=label_1000nm)
    max_threads = max(max(sp_th_100nm), max(sp_th_1000nm))
    axes[1].plot([1, max_threads], [1, max_threads], "--", linewidth=1.2, label="Ideal")
    axes[1].set_title("Speedup")
    axes[1].set_xlabel("Threads")
    axes[1].set_ylabel("Speedup (T1 / Tn)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    axes[2].plot(sp_th_100nm, ef_100nm, marker="o", linewidth=2, label=label_100nm)
    axes[2].plot(sp_th_1000nm, ef_1000nm, marker="s", linewidth=2, label=label_1000nm)
    axes[2].axhline(1.0, linestyle="--", linewidth=1.2, color="gray", label="Ideal")
    axes[2].set_title("Parallel Efficiency")
    axes[2].set_xlabel("Threads")
    axes[2].set_ylabel("Efficiency (Speedup / Threads)")
    axes[2].set_ylim(0, 1.1)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()

    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot parallel scaling using two run_summary.csv files."
    )
    parser.add_argument(
        "--csv-100nm",
        type=Path,
        default=Path("cross_100nm_test_threads/run_summary.csv"),
        help="Path to 100nm run_summary.csv",
    )
    parser.add_argument(
        "--csv-1000nm",
        type=Path,
        default=Path("cross_1000nm_test_threads/run_summary.csv"),
        help="Path to 1000nm run_summary.csv",
    )
    parser.add_argument(
        "--label-100nm",
        default="cross_100nm",
        help="Legend label for 100nm dataset",
    )
    parser.add_argument(
        "--label-1000nm",
        default="cross_1000nm",
        help="Legend label for 1000nm dataset",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("cross_thread_parallel_scaling.png"),
        help="Output PNG path",
    )
    args = parser.parse_args()

    plot_scaling(
        csv_100nm=args.csv_100nm,
        csv_1000nm=args.csv_1000nm,
        output_png=args.output,
        label_100nm=args.label_100nm,
        label_1000nm=args.label_1000nm,
    )
    print(f"Saved plot: {args.output.resolve()}")


if __name__ == "__main__":
    main()
