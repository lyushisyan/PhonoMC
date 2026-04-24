#!/usr/bin/env python3
"""Read convergence.txt and plot selected temperature traces versus time."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Sequence

import numpy as np

import matplotlib

if not os.environ.get("DISPLAY") and not os.environ.get("MPLBACKEND"):
    matplotlib.use("Agg")

import matplotlib.pyplot as plt


def _parse_header(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                return s[1:].strip().split()
            return s.split()
    raise ValueError(f"Empty file: {path}")


def _pick_evenly(indices: Sequence[int], n_pick: int) -> np.ndarray:
    if not indices:
        return np.array([], dtype=int)
    n = min(max(1, int(n_pick)), len(indices))
    pos = np.linspace(0, len(indices) - 1, n)
    picked = np.unique(np.round(pos).astype(int))
    while picked.size < n:
        for i in range(len(indices)):
            if i not in picked:
                picked = np.append(picked, i)
                if picked.size == n:
                    break
    picked.sort()
    return np.array([indices[i] for i in picked], dtype=int)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot selected temperature columns from convergence.txt")
    parser.add_argument("--input", default="convergence.txt", help="Path to convergence.txt")
    parser.add_argument("--n", type=int, default=20, help="Number of temperature columns to plot")
    parser.add_argument(
        "--out",
        default="",
        help="Output figure path (default: <input_stem>_temps20.png next to input)",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output figure DPI")
    parser.add_argument("--show-legend", action="store_true", help="Show legend for selected temperature IDs")
    args = parser.parse_args()

    in_path = Path(args.input).resolve()
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    if args.out:
        out_path = Path(args.out).resolve()
    else:
        out_path = in_path.with_name(f"{in_path.stem}_temps{int(args.n)}.png")

    header = _parse_header(in_path)
    if "time_ps" in header:
        time_idx = header.index("time_ps")
        time_label = "Time (ps)"
    elif "time" in header:
        time_idx = header.index("time")
        time_label = "Time"
    else:
        time_idx = None
        time_label = "Step"

    temp_indices = [i for i, c in enumerate(header) if c.startswith("T_")]
    if not temp_indices:
        raise ValueError(f"No temperature columns (T_*) found in: {in_path}")

    selected_temp_indices = _pick_evenly(temp_indices, args.n)
    selected_temp_names = [header[i] for i in selected_temp_indices]

    required_max_idx = max([i for i in selected_temp_indices] + ([time_idx] if time_idx is not None else [0]))

    time_vals: List[float] = []
    temp_series: List[List[float]] = [[] for _ in selected_temp_indices]

    with in_path.open("r", encoding="utf-8") as f:
        step = 0
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            arr = np.fromstring(s, sep=" ")
            if arr.size <= required_max_idx:
                continue

            if time_idx is not None:
                time_vals.append(float(arr[time_idx]))
            else:
                time_vals.append(float(step))

            for j, col_idx in enumerate(selected_temp_indices):
                temp_series[j].append(float(arr[col_idx]))
            step += 1

    if not time_vals:
        raise ValueError(f"No valid data rows found in: {in_path}")

    time_arr = np.asarray(time_vals, dtype=float)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=160)
    cmap = plt.get_cmap("tab20")
    for i, (name, y_list) in enumerate(zip(selected_temp_names, temp_series)):
        y = np.asarray(y_list, dtype=float)
        ax.plot(time_arr, y, lw=1.2, color=cmap(i % 20), label=name)

    ax.set_xlabel(time_label)
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"{in_path.name}: {len(selected_temp_indices)} Temperature Traces")
    ax.grid(alpha=0.3)

    if args.show_legend:
        ax.legend(ncol=2, fontsize=8, frameon=False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=int(args.dpi))
    plt.close(fig)

    print(f"[ok] input: {in_path}")
    print(f"[ok] selected_temperatures: {', '.join(selected_temp_names)}")
    print(f"[ok] output: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
