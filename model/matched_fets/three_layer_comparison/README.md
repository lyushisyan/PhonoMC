# FinFET 与三层 GAAFET 的统一热输运模型

2026-10-15：按用户选择，以已有 52×24×50 nm FinFET 为基准建立三层 GAAFET。
参考 NanoBTE（arXiv:2607.03226v3）的三层结构与热浴布置，但采用本项目的匹配尺寸和高斯热源；这不是该文算例的逐参数复现。

## 几何与匹配量

X 沿沟道，Y 为宽度，Z 为高度；所有几何单位 nm。

| 项目 | FinFET | GAAFET |
|---|---|---|
| 衬底 | 52×24×20 | 相同 |
| 源/漏（各） | 10×24×30 | 相同 |
| 沟道长度 | 32 | 32 |
| 沟道截面 | 6×20，单鳍 | 10×4，三层 |
| 总沟道截面积 | 120 nm² | 120 nm² |
| 总沟道体积 | 3840 nm³ | 3840 nm³ |
| 总硅体积 | 43200 nm³ | 43200 nm³ |
| 外包围盒 | 52×24×50 | 相同 |

源极 X=0–10，漏极 X=42–52；二者均覆盖 Y=0–24、Z=20–50。
FinFET 沟道为 X=10–42、Y=9–15、Z=20–40。
GAAFET 沟道为 X=10–42、Y=7–17，三层 Z 分别为：

- 下层：24.5–28.5；
- 中层：33–37；
- 上层：41.5–45.5。

衬底到下层、层间以及上层到源漏顶面的间隙均为 4.5 nm。
选择 4 nm 片厚是为了在增加到三层后保留与 FinFET 相同的总沟道截面积。
FinFET STL 与父目录原模型逐字节相同；三层 GAAFET 是新的连通、封闭布尔并集。
旧双层模型和旧热源片段保留，当前三层比较以本目录为准。

## 统一边界

| 表面 | 两种器件采用的相同条件 |
|---|---|
| 衬底底面 Z=0 | T=300 K，面积 1248 nm² |
| 源/漏顶面 Z=50 | T=300 K，各 240 nm²，总计 480 nm² |
| 衬底四个侧面，0<Z<20 | P，X/Y 周期配对，总面积 3040 nm² |
| 其余暴露表面（含片下表面、衬底露出顶面及源漏侧面） | R，η=1 nm |

基底 X=0↔52、Y=0↔24 的相对侧面周期配对，周期分别为 52、24 nm。
这是基底连接的周期重复单元，允许声子跨单元传播；周期面不作为外部热沉。
上部源漏侧面仍采用 R=1 nm，未设为周期面，因此它们不代表跨单元连续的硅连接。
热源是体积区，不是恒温边界，也不是内部反射界面。
两种碰撞模型均设置 `conservative_boundaries=true`，采用一致的固定入射热浴通量和弹性频率壳反射处理。

输入中的 R=1 nm 是粗糙度依赖的反射模型，不等于论文的严格完全漫反射。
为避免连接的共面侧壁跨越 Z=20 后被合并为一个 facet，显式使用
`merge_coplanar_facets=false`。边界选择覆盖全部三角面，再覆盖衬底侧面及热浴。

## 统一热源和总功率

基准为每个完整器件 **1 μW**，并非每层 1 μW。GAAFET 每层为 1/3 μW。
沿 X 的规定体积源为 exp[-(X−39)²/(2×1.5²)]，限于 X=36–42（±2σ）；
各沟道截面内均匀。FinFET 的 Y/Z 支撑区为 [9,15]×[20,40]，GAAFET 为
[7,17]×[24.5,45.5]，只对实际硅网格加热，片间空隙不参与归一化。
相等的总沟道截面积也使两者的高斯加权热源体积相等。
采用 `total_power` 而非固定峰值发热率，thermal 源谱、恒定时间包络。
这不是电子输运求得的焦耳热分布。

## 材料与输入

两者均采用 `material/Si_ph3py_m252525_T50-800` 的同一份 Si 25³ 全区模态数据，
含 50–800 K 的三声子 N/U 率；未加入栅极、氧化层、空气或异质材料界面。
材料哈希与 N+U 一致性检查在 `model.json`。

- `finfet_rta.toml`、`gaafet_rta.toml`：局域温度 RTA。
- `finfet_callaway.toml`、`gaafet_callaway.toml`：非线性 Callaway，局域温度 N/R 率。
- 两者初始和固定参考温度均为 300 K；10M 初始载体，dt=0.05 ps，40000 步（2 ns）。
- 网格 [52,24,100]，单元尺寸 1×1×0.5 nm³；准确对齐半整数的片层位置。
- 每 5 ps 输出，包含网格热流。以上是待验证的运行基准，不代表已收敛设置。

## 已完成的检查

`validation/solver_audit.json` 由本仓库实际 C++ 配置、域构建和热源模块生成，
覆盖四套输入的解析、实际边界分配、热浴面积、活动网格体积、热源支撑及分层功率。
周期检查包含双向配对、纯 X/Y 平移，以及每对三角面各顶点在平移后的精确对应；
不只检查面积或形状边长，避免三角剖分造成额外的切向位移。
每个器件有 86400 个活动单元，均为 0.5 nm³；加热单元数均为 1440。
GAAFET 的空隙中没有活动单元或热源功率；非接触沟道表面未被设为热浴。
STL 检查含封闭性、连通性、法向一致性及解析体积。

这些是模型前置检查，未运行载体输运、材料反射核初始化或稳态/收敛模拟。
正式运行前还需完成实际材料的零源平衡、源开启短测、反射核与能量账本检查。

## 重生成与复核

在仓库根目录执行（Python 需要 numpy、h5py、matplotlib）：

```sh
python3 tools/prepare_three_layer_fets.py --power-w 1e-6
c++ -std=c++17 -O2 -Iinclude tools/validate_three_layer_fets.cpp \
  src/SimulationConfig.cpp src/config/ConfigParsing.cpp \
  src/config/ConfigValidation.cpp src/config/OutputDirectory.cpp \
  src/SurfaceMesh.cpp src/SimulationDomain.cpp \
  src/geometry/MeshReader.cpp src/geometry/SurfaceSampling.cpp \
  src/geometry/DomainOutput.cpp src/solver/PrescribedHeatSource.cpp \
  -o /tmp/phonomc_validate_three_layer_fets
/tmp/phonomc_validate_three_layer_fets model/matched_fets/three_layer_comparison/validation/solver_audit.json model/matched_fets/three_layer_comparison/finfet_rta.toml model/matched_fets/three_layer_comparison/finfet_callaway.toml model/matched_fets/three_layer_comparison/gaafet_rta.toml model/matched_fets/three_layer_comparison/gaafet_callaway.toml
```

重生成或改变功率后必须重新执行复核；旧 audit 不能验证新输入。
`comparison.png` 与可编辑 `comparison.svg` 显示几何、热源及边界，均非温度计算结果。
