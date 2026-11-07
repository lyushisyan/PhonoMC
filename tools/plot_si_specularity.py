#!/usr/bin/env python3
"""Si 色散按 Ziman–Soffer p 着色，计算平均 p(omega) 和 p(T)。

依赖 numpy, matplotlib, h5py, phonopy, pyyaml。
默认读取配套的 Si 25^3 全网格及其 FC2；不修改力常数。
运行：python3 tools/plot_si_specularity.py --normal 0 0 1

q 是含 2*pi 的笛卡尔波矢 (nm^-1)，p=exp[-4*eta^2*(q.n)^2]。
频率平均：pbar(omega) = sum(w*p*delta)/sum(w*delta)。
温度平均：pbar_C(T) = sum(w*C(omega,T)*p)/sum(w*C(omega,T))，
其中 C/kB = x^2 exp(x)/(exp(x)-1)^2, x = hbar*omega/(kB*T)。
温度平均亦即 integral(D*pbar*C domega)/integral(D*C domega)。
默认 T=50--1000 K；使用固定谐性色散及粗糙度，不包含入射通量权重。
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import numpy as np
import phonopy
import yaml

ROOT = Path(__file__).resolve().parents[1]
ETAS = np.array([0.0, 0.1, 1.0])
H_OVER_KB_THz_K = 6.62607015e-34 * 1e12 / 1.380649e-23


def mode_heat_capacity_over_kb(frequency_thz, temperature_k):
    """C/kB；排除零频平移模态，使用 exp(-x) 形式避免大 x 溢出。"""
    frequency = np.asarray(frequency_thz, dtype=float)
    if not np.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("Temperature must be finite and positive.")
    if not np.isfinite(frequency).all() or np.any(frequency < 0):
        raise ValueError("Frequencies must be finite and nonnegative.")
    x = H_OVER_KB_THz_K * frequency / temperature_k
    capacity = np.zeros_like(x)
    active = (frequency > 1e-4) & (x < 700)
    capacity[active] = (
        (x[active] / (-np.expm1(-x[active])))**2 * np.exp(-x[active])
    )
    return capacity


def temperature_average(frequency, weights, probabilities, temperatures):
    """从未展宽的完整模态数据求和，避免 DOS 平滑影响温度平均。"""
    means = np.empty((len(probabilities), len(temperatures)))
    cv_per_cell = np.empty(len(temperatures))
    for i, temperature in enumerate(temperatures):
        capacity = mode_heat_capacity_over_kb(frequency, temperature)
        q_capacity = capacity.sum(axis=1) * weights
        total = q_capacity.sum()
        if total <= 0:
            raise ValueError("Temperature is too low for this finite q mesh.")
        means[:, i] = np.einsum("eq,q->e", probabilities, q_capacity) / total
        cv_per_cell[i] = total / weights.sum()
    np.testing.assert_allclose(means[0], 1, atol=1e-12)
    assert np.all((means >= 0) & (means <= 1 + 1e-12))
    assert np.all(np.diff(means, axis=0) <= 1e-12)
    return means, cv_per_cell


def mode_p(q_cart_nm, normal, eta):
    return np.exp(-4 * eta**2 * np.einsum("ij,j->i", q_cart_nm, normal)**2)


def smooth_histogram(values, sigma_bins):
    """归一化高斯核，镜像延拓在端点处保留积分。"""
    radius = max(1, int(np.ceil(4 * sigma_bins)))
    offsets = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (offsets / sigma_bins)**2)
    kernel /= kernel.sum()
    return np.convolve(np.pad(values, radius, mode="symmetric"), kernel, mode="valid")


def bz_averaged_p(q_fractional, reciprocal_nm, normal):
    """在第一 BZ 边界的等长代表元之间等分权重，避免方向选择偏差。

    搜索半径适用于此 Si fcc 原胞；检查最短代表元不在搜索边缘。
    """
    q = q_fractional - np.rint(q_fractional)
    shifts = np.array(list(itertools.product(range(-2, 3), repeat=3)))
    length2 = np.array([
        np.sum(np.einsum("ij,jk->ik", q + shift, reciprocal_nm)**2, axis=1)
        for shift in shifts
    ])
    minimum = length2.min(axis=0)
    tied = np.isclose(length2, minimum[None, :], rtol=1e-10, atol=1e-12)
    if np.any(tied[np.any(np.abs(shifts) == 2, axis=1)]):
        raise ValueError("First-BZ search radius is insufficient for this cell.")
    multiplicity = tied.sum(axis=0)
    probabilities = np.zeros((len(ETAS), len(q)))
    for shift, mask in zip(shifts, tied):
        if np.any(mask):
            cart = np.einsum("ij,jk->ik", q[mask] + shift, reciprocal_nm)
            for j, eta in enumerate(ETAS):
                probabilities[j, mask] += mode_p(cart, normal, eta)
    probabilities /= multiplicity[None, :]
    return probabilities, multiplicity


def spectra(frequency, weights, probabilities, bin_thz, sigma_thz):
    """先以完全相同的频率分箱和高斯核处理 D 与 Dp，再取比值。"""
    maximum = np.ceil((frequency.max() + 5 * sigma_thz) / bin_thz) * bin_thz
    edges = np.arange(round(maximum / bin_thz) + 1) * bin_thz
    f = frequency.ravel()
    mode_weights = np.broadcast_to(weights[:, None], frequency.shape).ravel()
    counts = np.histogram(f, edges, weights=mode_weights)[0]
    numerator = np.array([
        np.histogram(
            f, edges,
            weights=(weights[:, None] * p[:, None]
                     * np.ones_like(frequency)).ravel(),
        )[0] for p in probabilities
    ])
    # 反射核保留正频率积分；同一核用于分子与分母。
    smooth = smooth_histogram(counts, sigma_thz / bin_thz)
    smooth_p = np.array([
        smooth_histogram(row, sigma_thz / bin_thz) for row in numerator
    ])
    omega = 2 * np.pi * (edges[:-1] + edges[1:]) / 2
    domega = 2 * np.pi * bin_thz
    dos = smooth / weights.sum() / domega
    weighted_dos = smooth_p / weights.sum() / domega
    mean = np.full_like(weighted_dos, np.nan)
    # 屏蔽 DOS 极小或超出非零频率采样范围的区域，避免尾部比值误读。
    nonzero_min = frequency[frequency > 0.01].min()
    support = ((dos > 1e-3 * dos.max())
               & (omega >= 2 * np.pi * nonzero_min)
               & (omega <= 2 * np.pi * frequency.max()))
    np.divide(weighted_dos, dos, out=mean, where=support[None, :])
    np.testing.assert_allclose(np.sum(dos) * domega, frequency.shape[1], atol=1e-10)
    np.testing.assert_allclose(mean[0, support], 1, atol=1e-12)
    np.testing.assert_allclose(
        weighted_dos.sum(axis=1) * domega,
        (probabilities * weights).sum(axis=1) / weights.sum() * frequency.shape[1],
        atol=1e-10,
    )
    return omega, dos, weighted_dos, mean, counts


def save_figure(fig, out, name):
    for extension in ("png", "pdf"):
        fig.savefig(out / f"{name}.{extension}", dpi=240, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material-dir", type=Path,
                        default=ROOT / "material/Si_ph3py_m252525_T50-800")
    parser.add_argument("--harmonic-input", type=Path,
                        default=ROOT / "output/phono3py_20260910/Si/input")
    parser.add_argument("--normal", type=float, nargs=3, default=(0, 0, 1),
                        help="Cartesian normal; for cubic Si also the [hkl] direction.")
    parser.add_argument("--bin-thz", type=float, default=0.1)
    parser.add_argument("--sigma-thz", type=float, default=0.15)
    parser.add_argument("--temperature-min", type=float, default=50.0, help="K")
    parser.add_argument("--temperature-max", type=float, default=1000.0, help="K")
    parser.add_argument("--temperature-points", type=int, default=191)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/si_specularity")
    args = parser.parse_args()
    normal = np.array(args.normal, dtype=float)
    if not np.isfinite(normal).all() or np.linalg.norm(normal) == 0:
        parser.error("--normal must be finite and nonzero.")
    if not (np.isfinite(args.bin_thz) and np.isfinite(args.sigma_thz)
            and args.bin_thz > 0 and args.sigma_thz > 0):
        parser.error("Bin width and Gaussian sigma must be finite and positive.")
    if not (np.isfinite(args.temperature_min) and np.isfinite(args.temperature_max)
            and 0 < args.temperature_min < args.temperature_max
            and args.temperature_points >= 2):
        parser.error("Require 0 < temperature-min < temperature-max and at least 2 points.")
    normal /= np.linalg.norm(normal)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    provenance = json.loads((args.material_dir / "provenance.json").read_text())
    hashes = {}
    for name in ("fc2.hdf5", "phono3py.yaml"):
        hashes[name] = hashlib.sha256((args.harmonic_input / name).read_bytes()).hexdigest()
        if hashes[name] != provenance["input_sha256"][name]:
            raise ValueError(f"{name} does not match the material provenance.")
    source_yaml = yaml.safe_load((args.harmonic_input / "phono3py.yaml").read_text())
    factor = source_yaml["phono3py"]["frequency_unit_conversion_factor"]
    ph = phonopy.load(
        str(args.harmonic_input / "phono3py.yaml"),
        force_constants_filename=str(args.harmonic_input / "fc2.hdf5"),
        symmetrize_fc=False, is_nac=False, factor=factor,
    )
    cell = np.asarray(ph.primitive.cell)
    a = 2 * cell[0, 1]
    np.testing.assert_allclose(cell, (np.ones((3, 3)) - np.eye(3)) * a / 2, atol=1e-8)
    reciprocal_nm = 2 * np.pi * np.linalg.inv(cell).T * 10
    with h5py.File(args.material_dir / "kappa-fbz.hdf5") as h:
        q = h["qpoint"][:]
        frequency = h["frequency"][:]
        weights = h["weight"][:].astype(float)
        mesh = h["mesh"][:]
    if len(q) != np.prod(mesh) or not np.all(weights == 1):
        raise ValueError("A complete equal-weight full q mesh is required.")
    addresses = np.rint(q * mesh).astype(int) % mesh
    if len(np.unique(addresses, axis=0)) != len(q):
        raise ValueError("Full mesh has duplicated/missing q points.")
    if not np.isfinite(frequency).all() or frequency.min() < 0:
        raise ValueError("Material frequency data must be finite and nonnegative.")
    # 重新计算抽样点，核对色散和 DOS 使用相同的声子模型。
    sample = np.unique(np.linspace(0, len(q) - 1, 64).astype(int))
    ph.run_qpoints(q[sample])
    recalculated = ph.get_qpoints_dict()["frequencies"]
    # 原始 FC2 在 Gamma 有约 -0.007 THz 的微小数值虚频，导出表设为 0。
    if recalculated.min() < -0.01:
        raise ValueError("Unexpected imaginary frequencies.")
    frequency_error = float(np.max(np.abs(np.maximum(recalculated, 0) - frequency[sample])))
    np.testing.assert_allclose(np.maximum(recalculated, 0), frequency[sample], atol=5e-6)

    probabilities, multiplicity = bz_averaged_p(q, reciprocal_nm, normal)
    np.testing.assert_allclose(probabilities[0], 1, atol=1e-12)
    assert np.all((probabilities >= 0) & (probabilities <= 1))
    assert np.all(np.diff(probabilities, axis=0) <= 1e-12)
    omega, dos, weighted_dos, mean, counts = spectra(
        frequency, weights, probabilities, args.bin_thz, args.sigma_thz,
    )
    np.savetxt(
        out / "si_dos_specularity.csv",
        np.column_stack((omega, omega / (2 * np.pi), counts, dos,
                         weighted_dos.T, mean.T)),
        delimiter=",", comments="",
        header="omega_rad_ps,f_THz,raw_bin_mode_weight,DOS_per_rad_ps_per_cell,"
               "Dp_eta0,Dp_eta0.1nm,Dp_eta1nm,pmean_eta0,pmean_eta0.1nm,pmean_eta1nm",
    )
    temperatures = np.linspace(args.temperature_min, args.temperature_max, args.temperature_points)
    temperature_mean, cv_per_cell = temperature_average(
        frequency, weights, probabilities, temperatures,
    )
    np.savetxt(
        out / "si_temperature_specularity.csv",
        np.column_stack((temperatures, temperature_mean.T, cv_per_cell)),
        delimiter=",", comments="",
        header="temperature_K,pmean_C_eta0,pmean_C_eta0.1nm,pmean_C_eta1nm,Cv_per_cell_over_kB",
    )
    # 热容的高温极限为每个正频模态 kB，因此平均回到全谱态数平均。
    active_branches_per_q = np.sum(frequency > 1e-4, axis=1)
    active_weight = weights * active_branches_per_q
    high_temperature_limit = (
        np.einsum("eq,q->e", probabilities, active_weight) / active_weight.sum()
    )
    reference_temperatures = np.array([50., 100., 200., 300., 500., 800., 1000.])
    reference_means, _ = temperature_average(
        frequency, weights, probabilities, reference_temperatures,
    )

    # fcc 标准路径旋转到 X || [001]，笛卡尔坐标单位为 2*pi/a。
    nodes = np.array([[0, 0, 0], [0, 0, 1], [.5, 0, 1],
                      [.75, 0, .75], [0, 0, 0], [.5, .5, .5]])
    labels = [r"$\Gamma$", r"$X_z$", "$W$", "$K$", r"$\Gamma$", "$L$"]
    node_cart = nodes * 2 * np.pi / (a / 10)
    node_fractional = np.einsum("ij,jk->ik", node_cart, np.linalg.inv(reciprocal_nm))
    band_q = np.concatenate([
        np.linspace(left, right, 121)[:-1]
        for left, right in zip(node_fractional[:-1], node_fractional[1:])
    ] + [node_fractional[-1:]])
    ph.run_qpoints(band_q)
    band_raw = ph.get_qpoints_dict()["frequencies"]
    if band_raw.min() < -0.01:
        raise ValueError("Unexpected imaginary frequencies along the band path.")
    band_omega = 2 * np.pi * np.maximum(band_raw, 0)
    cart = np.einsum("ij,jk->ik", band_q, reciprocal_nm)
    distance = np.r_[0, np.cumsum(np.linalg.norm(np.diff(cart, axis=0), axis=1))]
    ticks = distance[np.arange(len(nodes)) * 120]
    band_p = np.array([mode_p(cart, normal, eta) for eta in ETAS])
    np.savez_compressed(
        out / "si_dispersion_specularity.npz",
        q_fractional=band_q, q_cartesian_nm_inverse=cart,
        path_distance_nm_inverse=distance, omega_rad_ps=band_omega,
        frequency_raw_THz=band_raw, p=band_p, eta_nm=ETAS,
        normal=normal, tick_positions=ticks, tick_labels=labels,
    )
    direction = "[" + " ".join(f"{v:g}" for v in args.normal) + "]"
    with plt.rc_context({"font.size": 11, "axes.labelsize": 12, "pdf.fonttype": 42}):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True, layout="constrained")
        for j, (ax, eta) in enumerate(zip(axes, ETAS)):
            for branch in range(band_omega.shape[1]):
                points = np.column_stack((distance, band_omega[:, branch]))
                lines = LineCollection(
                    np.stack((points[:-1], points[1:]), axis=1),
                    cmap="viridis", norm=Normalize(0, 1), linewidth=2.0,
                )
                lines.set_array((band_p[j, :-1] + band_p[j, 1:]) / 2)
                ax.add_collection(lines)
            ax.set(xlim=(distance[0], distance[-1]), ylim=(0, 100),
                   title=rf"$\eta={eta:g}$ nm", xlabel="Wave-vector path")
            ax.set_xticks(ticks, labels)
            ax.grid(axis="x", color="0.85", linewidth=.7)
        axes[0].set_ylabel(r"$\omega$ (rad ps$^{-1}$)")
        fig.colorbar(lines, ax=axes, fraction=.025, pad=.02, label="Specularity $p$")
        fig.suptitle(
            rf"Si | normal {direction} | $p=\exp[-4\eta^2(\mathbf{{q}}\cdot\hat{{n}})^2]$"
        )
        save_figure(fig, out, "si_dispersion_p")

        colors = ["#424242", "#167a9c", "#bf4c32"]
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), layout="constrained")
        for j, (eta, color) in enumerate(zip(ETAS, colors)):
            label = rf"$\eta={eta:g}$ nm"
            axes[0].plot(omega, weighted_dos[j], color=color, lw=1.8, label=label)
            axes[1].plot(omega, mean[j], color=color, lw=1.9, label=label)
        axes[0].set(
            title=r"$p$-weighted DOS $D_p(\omega)$  ($\eta=0$: total DOS)",
            ylabel=r"States / primitive cell / (rad ps$^{-1}$)",
        )
        axes[1].set(title=r"DOS average $\bar{p}(\omega)=D_p(\omega)/D(\omega)$",
                    ylabel=r"$\bar{p}$", ylim=(-.025, 1.05))
        for ax in axes:
            ax.set(xlim=(0, 100), xlabel=r"$\omega$ (rad ps$^{-1}$)")
            ax.grid(alpha=.2)
            ax.legend(frameon=False)
        axes[0].set_ylim(bottom=0)
        fig.suptitle(
            f"Si | normal {direction} | full {mesh[0]}×{mesh[1]}×{mesh[2]} q mesh"
            f" | Gaussian σ = {args.sigma_thz:g} THz"
        )
        save_figure(fig, out, "si_dos_p")

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), layout="constrained")
        for j, (eta, color) in enumerate(zip(ETAS, colors)):
            label = rf"$\eta={eta:g}$ nm"
            axes[0].plot(omega, mean[j], color=color, lw=2, label=label)
            axes[1].plot(temperatures, temperature_mean[j], color=color, lw=2, label=label)
        axes[0].set(
            title=r"(a) DOS average $\bar{p}(\omega)$",
            xlabel=r"$\omega$ (rad ps$^{-1}$)", ylabel=r"$\bar{p}(\omega)$",
            xlim=(0, 100),
        )
        axes[1].set(
            title=r"(b) Heat-capacity average $\bar{p}_C(T)$",
            xlabel=r"$T$ (K)", ylabel=r"$\bar{p}_C(T)$",
            xlim=(temperatures[0], temperatures[-1]),
        )
        for ax in axes:
            ax.set_ylim(-.025, 1.05)
            ax.grid(alpha=.2)
            ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(.55, .16))
        fig.suptitle(
            f"Si | normal {direction} | fixed harmonic phonons and roughness"
        )
        save_figure(fig, out, "si_mean_p_omega_temperature")

    metadata = {
        "material": str(args.material_dir.resolve()),
        "harmonic_input": str(args.harmonic_input.resolve()),
        "input_sha256": hashes, "phonopy_version": phonopy.__version__,
        "normal_cartesian_unit": normal.tolist(), "eta_nm": ETAS.tolist(),
        "mesh": mesh.tolist(), "bin_width_THz": args.bin_thz,
        "gaussian_sigma_THz": args.sigma_thz,
        "sampled_frequency_max_error_THz": frequency_error,
        "DOS_integral_states_per_primitive_cell": float(dos.sum() * 2 * np.pi * args.bin_thz),
        "BZ_boundary_q_count": int(np.count_nonzero(multiplicity > 1)),
        "BZ_boundary_convention": "Equal weights over all shortest reciprocal representatives.",
        "formula": "p(q,s)=exp[-4*eta_nm^2*(q_cart_nm^-1 dot normal)^2]",
        "average": "sum_qs w_q p_q delta(omega-omega_qs) / sum_qs w_q delta(omega-omega_qs)",
        "frequency_convention": "omega [rad/ps] = 2*pi*f [THz]",
        "temperature_average": "sum_qs w_q*C(omega_qs,T)*p_q / sum_qs w_q*C(omega_qs,T)",
        "mode_heat_capacity": "C/kB=x^2*exp(-x)/(1-exp(-x))^2; x=h*f/(kB*T)",
        "temperature_range_K": [args.temperature_min, args.temperature_max],
        "temperature_points": args.temperature_points,
        "temperature_frequency_cutoff_THz": 1e-4,
        "high_temperature_limit": high_temperature_limit.tolist(),
        "reference_temperature_K": reference_temperatures.tolist(),
        "reference_temperature_means_eta_by_T": reference_means.tolist(),
        "notes": [
            "Frequency curve: DOS-weighted average. Temperature curve: heat-capacity-weighted average.",
            "Temperature average uses unsmoothed full-grid modal sums, excluding zero-frequency translations.",
            "Fixed harmonic frequencies and eta; all T dependence comes from statistical heat-capacity weights.",
            "No incoming-flux weighting; very low temperatures require denser q sampling.",
            "Both DOS numerator and denominator use identical bins and Gaussian smoothing.",
            "Low-DOS tails and below-first-nonzero-sample averages are masked.",
            "Same q gives the same p for all six branches in this model.",
            "This is the bare Ziman-Soffer approximation, not atomistic reflection or the solver's corrected discrete reflection probability.",
            "X_z is (0,0,1) in conventional Cartesian units 2*pi/a.",
            "Tiny raw Gamma acoustic imaginary frequencies (above -0.01 THz) are displayed as zero; FC2 is unchanged.",
        ],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    print(f"Figures and data saved to {out.resolve()}")


if __name__ == "__main__":
    main()
