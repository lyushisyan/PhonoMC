#!/usr/bin/env python3
"""Plot 1D PhonoMC results for cross-plane and in-plane cases.

Outputs three PNG files:
- temperature_vs_time.png
- heatflux_vs_time.png
- kappa_vs_time.png
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import matplotlib

if not os.environ.get("DISPLAY") and not os.environ.get("MPLBACKEND"):
    matplotlib.use("Agg")

import matplotlib.pyplot as plt


def _read_convergence(path: Path) -> Tuple[List[str], np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    header_line = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                header_line = s
                break

    if header_line is None:
        raise ValueError(f"No header line found in: {path}")

    headers = header_line.lstrip("#").strip().split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.size == 0:
        raise ValueError(f"No numeric rows in: {path}")

    if data.shape[1] != len(headers):
        n = min(data.shape[1], len(headers))
        headers = headers[:n]
        data = data[:, :n]

    return headers, data


def _columns(headers: List[str], data: np.ndarray) -> Dict[str, np.ndarray]:
    return {h: data[:, i] for i, h in enumerate(headers)}


def _temperature_columns(cols: Dict[str, np.ndarray]) -> List[str]:
    names = [c for c in cols if c.startswith("T_")]

    def key(name: str) -> int:
        try:
            return int(name[2:])
        except ValueError:
            return 10**9

    return sorted(names, key=key)


def _time_axis(cols: Dict[str, np.ndarray]) -> Tuple[np.ndarray, str]:
    if "time_ps" in cols:
        return cols["time_ps"], "Time (ps)"
    if "timestep" in cols:
        return cols["timestep"], "Timestep"
    first = next(iter(cols.values()))
    return np.arange(first.size, dtype=float), "Record"


def _pick_evenly(names: List[str], count: int) -> List[str]:
    if count <= 0 or len(names) <= count:
        return names
    idx = np.unique(np.round(np.linspace(0, len(names) - 1, count)).astype(int))
    return [names[int(i)] for i in idx]


def _plot_temperatures(
    x: np.ndarray,
    xlabel: str,
    cols: Dict[str, np.ndarray],
    temp_cols: List[str],
    out_path: Path,
    max_lines: int,
    show_legend: bool,
    dpi: int,
) -> None:
    selected = _pick_evenly(temp_cols, max_lines)
    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=dpi)
    cmap = plt.get_cmap("viridis")
    denom = max(1, len(selected) - 1)
    for i, name in enumerate(selected):
        ax.plot(x, cols[name], lw=1.2, color=cmap(i / denom), label=name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature vs Time")
    ax.grid(alpha=0.28, linewidth=0.6)
    if show_legend:
        ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def _plot_single_series(
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    out_path: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=dpi)
    ax.plot(x, y, lw=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.28, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def _plot_kappa(
    x: np.ndarray,
    xlabel: str,
    cols: Dict[str, np.ndarray],
    out_path: Path,
    dpi: int,
) -> None:
    missing = [name for name in ("kappa_int", "kappa_eff") if name not in cols]
    if missing:
        raise ValueError(f"Missing thermal conductivity column(s): {', '.join(missing)}")

    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=dpi)
    ax.plot(x, cols["kappa_int"], lw=1.5, label="kappa_int")
    ax.plot(x, cols["kappa_eff"], lw=1.5, label="kappa_eff")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Thermal Conductivity (W/mK)")
    ax.set_title("Thermal Conductivity vs Time")
    ax.grid(alpha=0.28, linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_1d(result_dir: Path, out_dir: Path, max_temp_lines: int, show_legend: bool, dpi: int) -> None:
    headers, data = _read_convergence(result_dir / "convergence.txt")
    cols = _columns(headers, data)
    x, xlabel = _time_axis(cols)

    temp_cols = _temperature_columns(cols)
    if not temp_cols:
        raise ValueError("No temperature columns (T_*) found in convergence.txt.")
    if "heatflux" not in cols:
        raise ValueError("Missing heatflux column in convergence.txt.")

    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_temperatures(
        x=x,
        xlabel=xlabel,
        cols=cols,
        temp_cols=temp_cols,
        out_path=out_dir / "temperature_vs_time.png",
        max_lines=max_temp_lines,
        show_legend=show_legend,
        dpi=dpi,
    )
    _plot_single_series(
        x=x,
        y=cols["heatflux"],
        xlabel=xlabel,
        ylabel="Heat Flux (W/m^2)",
        title="Heat Flux vs Time",
        out_path=out_dir / "heatflux_vs_time.png",
        dpi=dpi,
    )
    _plot_kappa(
        x=x,
        xlabel=xlabel,
        cols=cols,
        out_path=out_dir / "kappa_vs_time.png",
        dpi=dpi,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot 1D cross-plane/in-plane PhonoMC results from convergence.txt."
    )
    parser.add_argument("result_dir", help="Result directory containing convergence.txt")
    parser.add_argument("--out-dir", default="", help="Output directory (default: <result_dir>/plots_1d)")
    parser.add_argument("--max-temp-lines", type=int, default=80, help="Maximum temperature traces to draw")
    parser.add_argument("--show-legend", action="store_true", help="Show temperature trace legend")
    parser.add_argument("--dpi", type=int, default=220, help="Figure DPI")
    args = parser.parse_args()

    result_dir = Path(args.result_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else result_dir / "plots_1d"
    plot_1d(
        result_dir=result_dir,
        out_dir=out_dir,
        max_temp_lines=args.max_temp_lines,
        show_legend=args.show_legend,
        dpi=args.dpi,
    )
    print(f"[ok] temperature: {out_dir / 'temperature_vs_time.png'}")
    print(f"[ok] heatflux: {out_dir / 'heatflux_vs_time.png'}")
    print(f"[ok] kappa: {out_dir / 'kappa_vs_time.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
