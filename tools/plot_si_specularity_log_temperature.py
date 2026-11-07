#!/usr/bin/env python3
"""从已收敛的模态积分结果绘制 1--100 K 的平均 p，温度轴取对数。

定义：
  p_qs = exp(-4 eta^2 q_z^2), [001] 法线；
  C_qs/kB = x^2 exp(x)/(exp(x)-1)^2, x = hbar*omega_qs/(kB*T)；
  pbar_C(T) = sum_qs(w_q*C_qs*p_qs) / sum_qs(w_q*C_qs)。

仅重绘已收敛的热容加权数据，不重新计算力常数或 q 网格。
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter
import numpy as np


def main():
    out = Path(__file__).resolve().parents[1] / "output/si_specularity_low_temperature"
    data = np.genfromtxt(out / "si_mean_p_T_1_300K.csv", delimiter=",", names=True)
    data = data[(data["T_K"] >= 1) & (data["T_K"] <= 100)]
    temperature = data["T_K"]
    if temperature[0] != 1 or temperature[-1] != 100:
        raise ValueError("Converged source data must include both 1 K and 100 K.")
    values = np.array([data[name] for name in ("p_eta0", "p_eta01nm", "p_eta1nm")])
    with plt.rc_context({"font.size": 12, "axes.labelsize": 14, "pdf.fonttype": 42}):
        fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
        for eta, probability, color in zip(
            (0, .1, 1), values, ("#424242", "#167a9c", "#bf4c32"),
        ):
            ax.plot(temperature, probability, lw=2.3, color=color,
                    label=rf"$\eta={eta:g}$ nm")
        ax.set(xscale="log", xlim=(1, 100), ylim=(0, 1.035),
               xlabel=r"$T$ (K, log scale)", ylabel=r"$\bar p_C(T)$",
               title="Si | [001] normal | heat-capacity weighted")
        ax.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50, 100]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(3, 4, 6, 7, 8, 9)))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.grid(which="major", alpha=.22)
        ax.grid(which="minor", axis="x", alpha=.08)
        ax.legend(frameon=False, loc="lower left")
        for extension in ("png", "pdf"):
            path = out / f"si_mean_p_T_1_100K_log.{extension}"
            fig.savefig(path, dpi=300, facecolor="white")
            print(path)
        plt.close(fig)
    np.savetxt(
        out / "si_mean_p_T_1_100K_log.csv",
        np.column_stack((temperature, values.T)), delimiter=",", comments="",
        header="T_K,p_eta0,p_eta0.1nm,p_eta1nm",
    )


if __name__ == "__main__":
    main()
