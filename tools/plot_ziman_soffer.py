#!/usr/bin/env python3
"""绘制 Ziman–Soffer 镜面反射概率，依赖 numpy 和 matplotlib。

    p = exp[-(4*pi*eta*cos(theta)/lambda)^2]

eta 是 RMS 粗糙度，与 lambda 均以 nm 为单位；theta 是相对法线的角度。
默认输出 eta = 0, 0.1, 1 nm 的三联图（PNG 和 PDF）。

用法：
    python3 tools/plot_ziman_soffer.py
    python3 tools/plot_ziman_soffer.py --lambda-min 0.1 --lambda-max 10 --xscale linear
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def specularity(wavelength_nm, theta_deg, eta_nm):
    """计算 p，支持 NumPy 广播；波长必须为正，角度范围为 [0, 90] 度。"""
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    theta_deg = np.asarray(theta_deg, dtype=float)
    eta_nm = np.asarray(eta_nm, dtype=float)
    if np.any(~np.isfinite(wavelength_nm)) or np.any(wavelength_nm <= 0):
        raise ValueError("Wavelength must be finite and positive.")
    if np.any(~np.isfinite(theta_deg)) or np.any((theta_deg < 0) | (theta_deg > 90)):
        raise ValueError("Theta must be finite and between 0 and 90 degrees.")
    if np.any(~np.isfinite(eta_nm)) or np.any(eta_nm < 0):
        raise ValueError("RMS roughness must be finite and nonnegative.")
    # 令掠入射处 cos(90°) 精确等于零，避免浮点余量。
    cosine = np.where(theta_deg == 90, 0.0, np.cos(np.deg2rad(theta_deg)))
    return np.exp(-np.square(4 * np.pi * eta_nm * cosine / wavelength_nm))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambda-min", type=float, default=0.1, help="Minimum wavelength (nm).")
    parser.add_argument("--lambda-max", type=float, default=100.0, help="Maximum wavelength (nm).")
    parser.add_argument("--xscale", choices=("linear", "log"), default="log")
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).resolve().parents[1] / "output" / "ziman_soffer",
    )
    args = parser.parse_args()
    if not (np.isfinite(args.lambda_min) and np.isfinite(args.lambda_max)
            and 0 < args.lambda_min < args.lambda_max):
        parser.error("Require finite 0 < --lambda-min < --lambda-max.")

    sample = np.geomspace if args.xscale == "log" else np.linspace
    wavelengths = sample(args.lambda_min, args.lambda_max, 801)
    angles = np.linspace(0, 90, 721)
    wavelength_grid, theta_grid = np.meshgrid(wavelengths, angles)

    with plt.rc_context({
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "pdf.fonttype": 42,
    }):
        fig, axes = plt.subplots(
            1, 3, figsize=(12.5, 4.2), sharex=True, sharey=True,
            layout="constrained",
        )
        for ax, eta in zip(axes, (0.0, 0.1, 1.0)):
            probability = specularity(wavelength_grid, theta_grid, eta)
            mesh = ax.pcolormesh(
                wavelengths, angles, probability, shading="auto",
                cmap="viridis", vmin=0, vmax=1, rasterized=True,
            )
            ax.set_xscale(args.xscale)
            ax.set_xlim(args.lambda_min, args.lambda_max)
            ax.set_ylim(0, 90)
            ax.set_yticks(np.arange(0, 91, 15))
            ax.set_xlabel(r"$\lambda$ (nm)")
            ax.set_title(rf"$\eta = {eta:g}\ \mathrm{{nm}}$")
            ax.tick_params(direction="out")

        axes[0].set_ylabel(r"$\theta$ (deg, from surface normal)")
        colorbar = fig.colorbar(mesh, ax=axes, pad=0.025, fraction=0.025)
        colorbar.set_label(r"Specularity $p$")
        colorbar.set_ticks(np.linspace(0, 1, 6))
        fig.suptitle(
            r"Ziman–Soffer:  $p=\exp[-(4\pi\eta\cos\theta/\lambda)^2]$",
            fontsize=15,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for extension in ("png", "pdf"):
            path = args.output_dir / f"ziman_soffer.{extension}"
            fig.savefig(path, dpi=300, facecolor="white")
            print(path.resolve())
        plt.close(fig)


if __name__ == "__main__":
    main()
