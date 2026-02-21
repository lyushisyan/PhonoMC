#!/usr/bin/env python3
"""Plot convergence diagnostics from ntmc convergence.txt files.

Per-case outputs:
- heatflux_convergence.png
- kappa_convergence.png
- temperature_convergence.png
- temperature_steady_tailN.png
- temperature_steady_tailN.csv

If multiple result dirs are passed, summary comparison figures are also produced.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
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

    if data.shape[1] != len(headers):
        n = min(data.shape[1], len(headers))
        headers = headers[:n]
        data = data[:, :n]

    if data.shape[0] == 0:
        raise ValueError(f"No numeric rows in: {path}")

    return headers, data


def _to_dict(headers: List[str], data: np.ndarray) -> Dict[str, np.ndarray]:
    return {h: data[:, i] for i, h in enumerate(headers)}


def _pick_column(cols: Dict[str, np.ndarray], candidates: List[str]) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


def _temp_columns(cols: Dict[str, np.ndarray]) -> List[str]:
    tcols = [c for c in cols if c.startswith("T_sv_")]

    def _key(name: str) -> int:
        try:
            return int(name.split("_")[-1])
        except ValueError:
            return 10**9

    tcols.sort(key=_key)
    return tcols


def _plot_case(result_dir: Path, tail: int, err_mode: str, dpi: int) -> Dict[str, np.ndarray]:
    conv = result_dir / "convergence.txt"
    headers, data = _read_convergence(conv)
    cols = _to_dict(headers, data)

    out_dir = result_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    tcol = _pick_column(cols, ["timestep", "Timestep"])
    if tcol is None:
        raise ValueError(f"Cannot find timestep column in {conv}")
    x = cols[tcol]

    hf_col = _pick_column(cols, ["heatflux", "Hflux"])
    if hf_col is None:
        raise ValueError(f"Cannot find heatflux column in {conv}")

    kfit_col = _pick_column(cols, ["kappa_fit"])
    kend_col = _pick_column(cols, ["kappa_end", "Kappa", "kappa"])
    temp_cols = _temp_columns(cols)
    if not temp_cols:
        raise ValueError(f"Cannot find temperature columns (T_sv_*) in {conv}")

    # 1) Heatflux convergence
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    ax.plot(x, cols[hf_col], lw=1.4)
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Heat Flux (W/m^2)")
    ax.set_title(f"Heat Flux Convergence: {result_dir.name}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "heatflux_convergence.png")
    plt.close(fig)

    # 2) Kappa convergence
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    if kfit_col is not None:
        ax.plot(x, cols[kfit_col], lw=1.4, label=kfit_col)
        # --- 修改 1: y轴最大值设置为 kappa_fit 最终值的 1.5 倍 ---
        if cols[kfit_col].size > 0:
            final_kappa = cols[kfit_col][-1]
            if not np.isnan(final_kappa) and final_kappa > 0:
                ax.set_ylim(0, final_kappa * 1.5)
                
    if kend_col is not None:
        ax.plot(x, cols[kend_col], lw=1.4, label=kend_col)
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Thermal Conductivity")
    ax.set_title(f"Kappa Convergence: {result_dir.name}")
    if (kfit_col is not None) or (kend_col is not None):
        ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "kappa_convergence.png")
    plt.close(fig)

    # 3) Temperature convergence
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    # --- 修改 2: 使用渐变色 (如 viridis) 代替 tab10 ---
    cmap = plt.get_cmap("viridis") 
    n_lines = len(temp_cols)
    for i, c in enumerate(temp_cols):
        # 根据子体积索引分配渐变色
        color = cmap(i / (n_lines - 1)) if n_lines > 1 else cmap(0.5)
        ax.plot(x, cols[c], lw=1.0, color=color, label=c)
    
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"Subvolume Temperature Convergence: {result_dir.name}")
    ax.grid(alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "temperature_convergence.png")
    plt.close(fig)

    # 4) Steady temperature distribution from last N records
    n_tail = max(1, min(tail, data.shape[0]))
    temp_matrix = np.column_stack([cols[c][-n_tail:] for c in temp_cols])
    mean = temp_matrix.mean(axis=0)
    if temp_matrix.shape[0] > 1:
        std = temp_matrix.std(axis=0, ddof=1)
    else:
        std = np.zeros_like(mean)
    if err_mode == "sem":
        err = std / math.sqrt(temp_matrix.shape[0])
    else:
        err = std

    # --- 修改 3: 横坐标索引加 1 (从 1 开始) ---
    idx = np.arange(1, len(temp_cols) + 1)
    
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    ax.errorbar(idx, mean, yerr=err, fmt="o-", capsize=4, lw=1.2, mfc='none', mec='C0', mew=1.2)
    ax.set_xlabel("Subvolume Index")
    ax.set_ylabel("Temperature (K)")
    
    # --- 修改 4: 纵坐标设置默认 299 到 301 ---
    ax.set_ylim(299, 301)
    
    ax.set_title(f"Steady Temperature (last {n_tail} steps): {result_dir.name}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"temperature_steady_tail{n_tail}.png")
    plt.close(fig)

    csv_path = out_dir / f"temperature_steady_tail{n_tail}.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("subvolume,mean_T,error,std\n")
        for i in range(len(temp_cols)):
            f.write(f"{i+1},{mean[i]:.10g},{err[i]:.10g},{std[i]:.10g}\n")

    return {
        "name": np.array([result_dir.name]),
        "x": x,
        "heatflux": cols[hf_col],
        "kappa_fit": cols[kfit_col] if kfit_col is not None else np.array([]),
        "kappa_end": cols[kend_col] if kend_col is not None else np.array([]),
        "steady_mean": mean,
        "steady_err": err,
    }


def _plot_summary(cases: List[Dict[str, np.ndarray]], summary_dir: Path, dpi: int) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)

    # heatflux comparison
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=dpi)
    for c in cases:
        ax.plot(c["x"], c["heatflux"], lw=1.3, label=c["name"][0])
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Heat Flux (W/m^2)")
    ax.set_title("Heat Flux Convergence Comparison")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "compare_heatflux.png")
    plt.close(fig)

    # kappa comparison
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=dpi)
    for c in cases:
        if c["kappa_end"].size:
            ax.plot(c["x"], c["kappa_end"], lw=1.3, label=f"{c['name'][0]}: kappa_end")
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Kappa")
    ax.set_title("Kappa_end Convergence Comparison")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "compare_kappa_end.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=dpi)
    any_fit = False
    for c in cases:
        if c["kappa_fit"].size:
            any_fit = True
            ax.plot(c["x"], c["kappa_fit"], lw=1.3, label=f"{c['name'][0]}: kappa_fit")
    if any_fit:
        ax.set_xlabel("Timestep")
        ax.set_ylabel("Kappa")
        ax.set_title("Kappa_fit Convergence Comparison")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(summary_dir / "compare_kappa_fit.png")
    plt.close(fig)

    # steady temperature profile comparison
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=dpi)
    for c in cases:
        idx = np.arange(c["steady_mean"].size)
        ax.errorbar(idx, c["steady_mean"], yerr=c["steady_err"], fmt="o-", capsize=3, lw=1.1, label=c["name"][0])
    ax.set_xlabel("Subvolume Index")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Steady Temperature Profile Comparison")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(summary_dir / "compare_steady_temperature.png")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot convergence and steady-state temperature from ntmc results.")
    parser.add_argument("result_dirs", nargs="+", help="Result directories containing convergence.txt")
    parser.add_argument("--tail", type=int, default=50, help="Number of last records for steady-state averaging")
    parser.add_argument("--error", choices=["sem", "std"], default="sem", help="Error bar type for steady temperature")
    parser.add_argument("--dpi", type=int, default=220, help="Figure DPI")
    parser.add_argument("--summary-dir", default="", help="Directory for cross-case comparison figures")
    args = parser.parse_args()

    case_data: List[Dict[str, np.ndarray]] = []
    for rd in args.result_dirs:
        result_dir = Path(rd).expanduser().resolve()
        d = _plot_case(result_dir, tail=args.tail, err_mode=args.error, dpi=args.dpi)
        case_data.append(d)
        print(f"[ok] plotted case: {result_dir}")

    if len(case_data) > 1:
        if args.summary_dir:
            summary_dir = Path(args.summary_dir).expanduser().resolve()
        else:
            summary_dir = (Path(args.result_dirs[0]).expanduser().resolve().parent / "plots_summary")
        _plot_summary(case_data, summary_dir=summary_dir, dpi=args.dpi)
        print(f"[ok] summary plots: {summary_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
