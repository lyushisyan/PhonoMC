#!/usr/bin/env python3
"""Si [001] 的热容加权平均 p(T)，1--300 K，含 Gamma 加密与收敛检查。

沿用 plot_si_specularity.py 的材料、FC2、eta 和热容加权定义。
在内存中对 FC2 施加平移/置换对称性；原文件不变。
使用 fcc 第一 BZ 的 1/48 楔形区，q=(t,t*u,t*u*v)*2*pi/a。
0<u,v<1, 0<t<min(1,1.5/(1+u+u*v))；Jacobian=t^2*u。
径向分成 [0,2^-8,2^-7,...,1] 的几何区间，逐层加密 Gamma。
对立方对称星，p_[001] 平均为 sum_i exp(-4*eta^2*q_i^2)/3。
Gauss-Legendre 节点具有体积积分权重；不是把非均匀节点等权平均。

运行：
  OMP_NUM_THREADS=1 python3 tools/plot_si_specularity_low_temperature.py
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
import numpy as np
import phonopy
from numpy.polynomial.legendre import leggauss

from plot_si_specularity import (
    ROOT, ETAS, H_OVER_KB_THz_K, bz_averaged_p, save_figure,
)

MATERIAL = ROOT / "material/Si_ph3py_m252525_T50-800"
INPUT = ROOT / "output/phono3py_20260910/Si/input"
RADIAL_EDGES = np.r_[0., 2. ** np.arange(-8, 1)]
LEVELS = [(12, 6), (20, 8), (32, 10), (48, 12)]


def quadrature(angular_order, radial_order):
    """返回笛卡尔 q/(2*pi/a) 及全 BZ 归一化体积权重。"""
    nodes, weights = leggauss(angular_order)
    nodes, weights = (nodes + 1) / 2, weights / 2
    rn, rw = leggauss(radial_order)
    r = np.concatenate([
        lo + (rn + 1) * (hi - lo) / 2
        for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])
    ])
    r_weight = np.concatenate([
        rw * (hi - lo) / 2
        for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])
    ])
    points, volume_weights = [], []
    for v, wv in zip(nodes, weights):
        # 在立方体面和八面体面相交处拆分，避免积分跨过 tmax 的折点。
        split = .5 / (1 + v)
        for ulo, uhi in [(0., split), (split, 1.)]:
            u = ulo + nodes * (uhi - ulo)
            wu = weights * (uhi - ulo)
            tmax = np.minimum(1., 1.5 / (1 + u + u * v))
            t = tmax[:, None] * r[None, :]
            point = np.stack((t, t * u[:, None], t * u[:, None] * v), axis=-1)
            w = (wv * wu[:, None] * u[:, None] * tmax[:, None]**3
                 * r[None, :]**2 * r_weight[None, :])
            points.append(point.reshape(-1, 3))
            volume_weights.append(w.ravel())
    points = np.concatenate(points)
    # 全 BZ 体积是 4*(2*pi/a)^3，48 个等价楔形区 -> 归一化系数 12。
    weights = 12 * np.concatenate(volume_weights)
    np.testing.assert_allclose(weights.sum(), 1., atol=1e-11)
    assert np.all(points >= 0) and np.all(points[:, 0] <= 1)
    assert np.all(points.sum(axis=1) <= 1.5)
    return points, weights


def thermal_average(frequency, weights, probability, temperatures):
    """返回 p(T) 和 C_V/kB/原胞；正频模态均参与，没有频率截断。"""
    means = np.empty((3, len(temperatures)))
    total_cv = np.empty(len(temperatures))
    for j, temperature in enumerate(temperatures):
        x = H_OVER_KB_THz_K * frequency / temperature
        capacity = np.zeros_like(x)
        active = (x > 0) & (x < 700)
        xx = x[active]
        capacity[active] = (xx / -np.expm1(-xx))**2 * np.exp(-xx)
        q_weight = weights * capacity.sum(axis=1)
        total_cv[j] = q_weight.sum()
        if total_cv[j] <= 0:
            raise ValueError("No thermally active modes.")
        means[:, j] = np.einsum("eq,q->e", probability, q_weight) / total_cv[j]
    np.testing.assert_allclose(means[0], 1, atol=2e-12)
    assert np.all((means >= 0) & (means <= 1 + 2e-12))
    return means, total_cv


def evaluate(ph, q_fractional):
    result = []
    for start in range(0, len(q_fractional), 8192):
        ph.run_qpoints(q_fractional[start:start + 8192])
        result.append(ph.get_qpoints_dict()["frequencies"].copy())
    return np.concatenate(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/si_specularity_low_temperature")
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((MATERIAL / "provenance.json").read_text())
    hashes = {
        name: hashlib.sha256((INPUT / name).read_bytes()).hexdigest()
        for name in ("fc2.hdf5", "phono3py.yaml")
    }
    for name, digest in hashes.items():
        assert digest == provenance["input_sha256"][name]
    ph = phonopy.load(
        str(INPUT / "phono3py.yaml"),
        force_constants_filename=str(INPUT / "fc2.hdf5"),
        symmetrize_fc=False, is_nac=False, factor=15.633302,
    )
    ph.run_qpoints([[0., 0., 0.]])
    raw_gamma = ph.get_qpoints_dict()["frequencies"][0].copy()
    original_fc = ph.force_constants.copy()
    ph.symmetrize_force_constants()
    fc_max_change = float(np.max(np.abs(ph.force_constants - original_fc)))
    ph.run_qpoints([[0., 0., 0.]])
    corrected_gamma = ph.get_qpoints_dict()["frequencies"][0].copy()
    assert np.max(np.abs(corrected_gamma[:3])) < 1e-5
    cell = np.array(ph.primitive.cell)
    a = 2 * cell[0, 1]
    np.testing.assert_allclose(cell, (np.ones((3, 3)) - np.eye(3)) * a / 2, atol=1e-10)
    reciprocal = 2 * np.pi * np.linalg.inv(cell).T * 10
    cart_to_fractional = np.linalg.inv(reciprocal)
    scale = 2 * np.pi / (a / 10)
    # 检查施加 ASR 后仍保留立方对称性，才能用 1/48 楔形区积分。
    representative = np.array([.23, .12, .07]) * scale
    orbit = np.array([
        representative[list(perm)] * signs
        for perm in itertools.permutations(range(3))
        for signs in itertools.product((-1, 1), repeat=3)
    ])
    orbit_f = evaluate(ph, np.einsum("ij,jk->ik", orbit, cart_to_fractional))
    symmetry_error = float(np.max(np.abs(orbit_f - orbit_f[0])))
    np.testing.assert_allclose(orbit_f, np.broadcast_to(orbit_f[0], orbit_f.shape), atol=1e-7)

    reference_t = np.array([1., 2., 5., 10., 20., 30., 50., 100., 200., 300.])
    temperatures = np.unique(np.r_[np.geomspace(1, 40, 100), np.linspace(40, 300, 66), reference_t])
    results = []
    for angular, radial in LEVELS:
        points, weights = quadrature(angular, radial)
        cart = points * scale
        q = np.einsum("ij,jk->ik", cart, cart_to_fractional)
        frequency = evaluate(ph, q)
        if np.min(frequency) <= 0:
            raise ValueError("Nonpositive frequency on the refined grid after ASR.")
        probability = np.array([
            np.exp(-4 * eta**2 * cart**2).mean(axis=1) for eta in ETAS
        ])
        mean, cv = thermal_average(frequency, weights, probability, temperatures)
        np.savez_compressed(
            out / f"quadrature_a{angular}_r{radial}.npz",
            q_fractional=q, q_cart_nm_inverse=cart, weight=weights,
            frequency_THz=frequency, probability_star_average=probability,
            temperature_K=temperatures, mean_p=mean, Cv_over_kB_per_cell=cv,
        )
        results.append({
            "angular_order": angular, "radial_order_per_shell": radial,
            "q_points_in_wedge": len(q), "min_q_nm_inverse": float(np.linalg.norm(cart, axis=1).min()),
            "min_frequency_THz": float(frequency.min()), "mean": mean, "cv": cv,
        })
        print(f"a={angular}, r={radial}, {len(q)} q: p(1 K)={mean[:,0]}, Cv/kB={cv[0]:.8g}", flush=True)

    # 对比原始 25^3 网格和同一网格的 ASR 结果，分离采样与力常数修正的影响。
    with h5py.File(MATERIAL / "kappa-fbz.hdf5") as h:
        q25 = h["qpoint"][:]
        f25 = h["frequency"][:]
    p25, _ = bz_averaged_p(q25, reciprocal, np.array([0., 0., 1.]))
    w25 = np.full(len(q25), 1 / len(q25))
    old_mean, old_cv = thermal_average(f25, w25, p25, temperatures)
    f25_asr = evaluate(ph, q25)
    f25_asr[np.linalg.norm(q25, axis=1) < 1e-10, :3] = 0
    asr25_mean, asr25_cv = thermal_average(f25_asr, w25, p25, temperatures)
    final = results[-1]
    previous = results[-2]
    p_error = np.abs(final["mean"] - previous["mean"])
    cv_error = np.abs(final["cv"] / previous["cv"] - 1)
    convergence = {
        "max_absolute_mean_p_change": float(p_error.max()),
        "max_relative_Cv_change": float(cv_error.max()),
        "mean_p_absolute_tolerance": 5e-4, "Cv_relative_tolerance": 5e-3,
        "passed": bool(p_error.max() < 5e-4 and cv_error.max() < 5e-3),
    }
    if not convergence["passed"]:
        raise RuntimeError(f"Refine quadrature further: {convergence}")
    mean = final["mean"]
    np.savetxt(
        out / "si_mean_p_T_1_300K.csv",
        np.column_stack((temperatures, mean.T, final["cv"],
                         old_mean[1:].T, asr25_mean[1:].T, old_cv, asr25_cv)),
        delimiter=",", comments="",
        header="T_K,p_eta0,p_eta0.1nm,p_eta1nm,Cv_over_kB_per_cell,"
               "raw25_p_eta0.1nm,raw25_p_eta1nm,asr25_p_eta0.1nm,asr25_p_eta1nm,"
               "raw25_Cv_over_kB_per_cell,asr25_Cv_over_kB_per_cell",
    )
    # 低温声学极限：C~T^3、1-p~T^2，检查是否被数值误差污染。
    two = int(np.flatnonzero(temperatures == 2)[0])
    low_t_checks = {
        "Cv_2K_over_Cv_1K": float(final["cv"][two] / final["cv"][0]),
        "one_minus_p_eta0.1_2K_over_1K": float((1 - mean[1, two]) / (1 - mean[1, 0])),
    }
    with plt.rc_context({"font.size": 11, "axes.labelsize": 12, "pdf.fonttype": 42}):
        colors = ["#424242", "#167a9c", "#bf4c32"]
        fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5), layout="constrained")
        for eta, values, color in zip(ETAS, mean, colors):
            for ax in axes:
                ax.plot(temperatures, values, color=color, lw=2, label=rf"$\eta={eta:g}$ nm")
        axes[0].set(xlim=(1, 300), title="1–300 K")
        axes[1].set(xlim=(1, 30), title="Low-temperature detail: 1–30 K")
        axes[0].set_xticks([1, 50, 100, 150, 200, 250, 300])
        axes[1].set_xticks([1, 5, 10, 15, 20, 25, 30])
        for ax in axes:
            ax.set(xlabel="$T$ (K)", ylabel=r"$\bar p_C(T)$", ylim=(-.025, 1.035))
            ax.grid(alpha=.2)
            ax.legend(frameon=False, loc="center right")
        axes[0].legend(frameon=False, loc="center", bbox_to_anchor=(.75, .3))
        fig.suptitle("Si | [001] normal | heat-capacity weighted | refined q integration")
        save_figure(fig, out, "si_mean_p_T_1_300K")

        fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5), layout="constrained")
        axes[0].plot(temperatures, old_mean[2], "--", color="0.55", lw=1.8, label=r"Original $25^3$")
        axes[0].plot(temperatures, asr25_mean[2], ":", color="0.15", lw=1.8, label=r"$25^3$ + ASR")
        for result in results:
            label = f"a={result['angular_order']}, r={result['radial_order_per_shell']}"
            axes[0].plot(temperatures, result["mean"][2], lw=1.1, label=label)
            axes[1].semilogx(temperatures, result["cv"] / final["cv"] - 1, lw=1.1, label=label)
        axes[0].set(xlim=(1, 30), ylim=(0, 1), title=r"Mean $p$: $\eta=1$ nm",
                    xlabel="$T$ (K)", ylabel=r"$\bar p_C$")
        axes[1].set(xlim=(1, 300), title="Heat capacity convergence",
                    xlabel="$T$ (K)", ylabel=r"$C_V/C_{V,\mathrm{finest}}-1$")
        for ax in axes:
            ax.grid(alpha=.2)
            ax.legend(frameon=False, fontsize=9)
        save_figure(fig, out, "si_low_temperature_convergence")

    metadata = {
        "input_sha256": hashes, "phonopy_version": phonopy.__version__,
        "normal": [0, 0, 1], "eta_nm": ETAS.tolist(),
        "raw_Gamma_frequency_THz": raw_gamma.tolist(),
        "ASR_Gamma_frequency_THz": corrected_gamma.tolist(),
        "max_FC2_change_eV_A2": fc_max_change,
        "FC2_change_relative_to_max": fc_max_change / float(np.max(np.abs(original_fc))),
        "ASR_cubic_symmetry_frequency_error_THz": symmetry_error,
        "quadrature": "48 symmetry-related wedges; q=(t,tu,tuv)*2pi/a, Jacobian=t^2*u; normalized weights",
        "radial_shell_edges_fraction_of_BZ_boundary": RADIAL_EDGES.tolist(),
        "levels": [{k: v for k, v in row.items() if k not in ("mean", "cv")} for row in results],
        "convergence": convergence, "low_T_checks": low_t_checks,
        "reference_T_K": reference_t.tolist(),
        "reference_mean_p_eta_by_T": mean[:, np.searchsorted(temperatures, reference_t)].tolist(),
        "definition": "pbar_C(T)=sum_qs w_q*C_qs(T)*p_qs / sum_qs w_q*C_qs(T)",
        "notes": [
            "Original FC2 files are unchanged; ASR/permutation symmetrization is in memory.",
            "Eta and harmonic phonon dispersion are fixed; temperature changes heat-capacity weights only.",
            "Nonuniform cubature has volume weights; no Debye replacement or frequency smoothing.",
            "Cubic star averaging is valid for the specified [001] normal.",
            "Uniform 25^3 curves are comparison only and are not valid at 1 K.",
        ],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
