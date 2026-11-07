# 中文手册配套输入

入口：[中文用户手册](../PhonoMC_中文用户手册.md) · [截图指南](../PhonoMC_截图与审阅指南.md)

所有输入从项目根目录使用以下形式运行：

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/01_cross_smoke.toml
```

输入材料和几何路径已按本目录调整；结果写入 `docs/manual_runs/案例名_序号`。以运行打印的路径为准。运行前按手册准备依赖与可执行文件。

| 输入 | 用途 |
|---|---|
| [01_cross_smoke.toml](01_cross_smoke.toml) | 20 ps 短流程检查 |
| [02_cross_reference.toml](02_cross_reference.toml) | 跨面 Si 参考设置 |
| [03_inplane_finite.toml](03_inplane_finite.toml) | 有限长度面内薄膜 |
| [04_inplane_gradient.toml](04_inplane_gradient.toml) | 无限重复薄膜梯度驱动 |
| [05_bilayer_dmm.toml](05_bilayer_dmm.toml) | Si/Ge DMM 界面 |
| [06_bilayer_mmm.toml](06_bilayer_mmm.toml) | Si/Ge MMM 界面 |
| [07_total_power.toml](07_total_power.toml) | 定总功率热源与账本 |
| [08_gaussian_spectrum.toml](08_gaussian_spectrum.toml) | Gaussian 频谱热源 |
| [09_finfet.toml](09_finfet.toml) | FinFET 三维局域加热 |
| [10_callaway_linearized.toml](10_callaway_linearized.toml) | 合成线性 Callaway 演示 |
| [11_gradient_stationarity.toml](11_gradient_stationarity.toml) | 梯度热启动/平稳性筛查 |
| [12_square_source.toml](12_square_source.toml) | 方波热源 |
| [13_callaway_nonlinear.toml](13_callaway_nonlinear.toml) | 固定 300 K N/R 率的合成非线性 Callaway 演示 |

Callaway 两例需先生成合成数据：

```bash
python tools/create_callaway_demo.py --out docs/manual_runs/callaway_demo_material
```

已存在材料时不重复生成。合成输出不代表真实材料，所有教学算例均未声明物理收敛。保留输入快照和实际结果路径；修改参数时使用不同的输出基名。

输入校验范围见 [VALIDATION.md](VALIDATION.md)。
