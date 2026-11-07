# 中文手册本地检查记录

检查日期：2026-10-23。对象为本次新建中文手册、截图指南及 13 个配套 TOML。没有执行版本控制或上传操作。

## 已完成

- 对照本地 `include/SimulationConfig.h` 和 `src/SimulationConfig.cpp` 核对默认值、识别键与组合限制；手册覆盖全部 76 个识别键，明确标注内部保留/公开不可用项。
- 使用 Python 标准库 `tomllib` 解析全部 13 个输入和手册的 TOML 代码块。
- 将当前 `src/SimulationConfig.cpp` 与临时小型驱动程序编译到系统临时目录，逐一调用真实的 `load_simulation_config`；13 个输入全部通过当前 C++ 配置校验。
- 检查材料和几何相对路径；真实材料与 FinFET 模型路径存在。两份 Callaway 演示共同依赖的 `docs/manual_runs/callaway_demo_material` 需要读者先按手册命令生成，此项明确作为前置条件。
- 检查所有配套输出基名解析到 `docs/manual_runs`。
- 核对快速入门完整 TOML 与 `01_cross_smoke.toml` 的配置一致。
- 检查内嵌 Python 示例语法、Markdown 代码围栏、章节锚点与本地链接。
- 对照当前结果写入代码核实表头、网格 nm 与热源账本 Å 的区别、单步诊断与累计账本的区别。

| 输入 | C++ 配置校验 | 最大模拟时间 |
|---|---|---:|
| 01_cross_smoke | 通过 | 20 ps |
| 02_cross_reference | 通过 | 1000 ps |
| 03_inplane_finite | 通过 | 1000 ps |
| 04_inplane_gradient | 通过 | 2000 ps |
| 05_bilayer_dmm | 通过 | 1000 ps |
| 06_bilayer_mmm | 通过 | 1000 ps |
| 07_total_power | 通过 | 100 ps |
| 08_gaussian_spectrum | 通过 | 10 ps |
| 09_finfet | 通过 | 100 ps |
| 10_callaway_linearized | 通过 | 2 ps |
| 11_gradient_stationarity | 通过 | 15000 ps，允许平稳性提前停止 |
| 12_square_source | 通过 | 50 ps |
| 13_callaway_nonlinear | 通过 | 2 ps |

## 验证范围

配置校验不会加载完整材料、创建运行几何或推进求解器，因此不替代运行时材料/面片/界面/模态兼容性检查。未为编写本手册启动长时间输运模拟，未认证任何教学输入的物理收敛，也未生成用户操作截图。

安装、构建、实际运行和图形结果由读者按截图指南逐步验证。若本地编译环境或数据发生变化，应重新检查；这里的结果不自动适用于另一个源码版本。
