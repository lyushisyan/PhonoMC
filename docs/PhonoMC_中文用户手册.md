# PhonoMC 中文用户手册

**本地审阅版 · 2026-10-23**

适用对象：需要使用 PhonoMC 计算薄膜、层状材料与器件声子热输运的研究人员和学生。

本手册依据当前工作目录中的源代码、输入示例和英文文档编写。它描述的是这份本地代码的行为，不代表某个已发布版本。操作以 macOS/Linux 终端为主。文中的命令需要由读者逐步执行；给出命令、通过输入校验与完成物理验证是三件不同的事。

**审阅与截图约定：**正文使用 `S01` 等编号标记截图点，每处说明操作、截图时机和通过标准。请配合 [截图与审阅指南](PhonoMC_截图与审阅指南.md) 使用。示例输入放在 [user_manual_examples](user_manual_examples/) 中。截图由实际操作产生，不用示意图替代运行证据。

## 目录

1. [软件简介与适用范围](#chapter-1)
2. [安装与编译](#chapter-2)
3. [第一个计算：快速入门](#chapter-3)
4. [基本概念与模型选择](#chapter-4)
5. [材料数据准备](#chapter-5)
6. [几何、网格与边界条件](#chapter-6)
7. [输入文件参数详解](#chapter-7)
8. [热源与材料界面](#chapter-8)
9. [典型算例教程](#chapter-9)
10. [输出文件与后处理](#chapter-10)
11. [收敛性与结果验证](#chapter-11)
12. [性能与故障排查](#chapter-12)
13. [附录：单位、术语、公式与复现记录](#appendix)

<a id="chapter-1"></a>
## 1. 软件简介与适用范围

### 1.1 PhonoMC 能解决什么问题

PhonoMC 是使用 C++17 编写的声子蒙特卡洛热输运求解器。程序从材料文件读取逐模式频率、群速度与散射线宽，跟踪计算载体在几何域中的运动、碰撞、边界反射和热浴交换，重建温度、热流及适用条件下的热导率。

输入材料可以来自第一性原理计算，也可以来自经过验证的其他晶格动力学数据源。随项目提供的 3C-SiC 数据来自 SiC 专用 NEP 势，因此不能把所有附带材料统一称为直接 DFT 数据。

| 研究目标 | 推荐设置 | 主要观察量 |
|---|---|---|
| 跨面薄膜输运 | 输运方向两端热浴，横向周期边界 | 温度分布、热流、有效热导率 |
| 有限长度面内输运 | 两端热浴，厚度方向粗糙表面 | 接触温降、长度与厚度效应 |
| 无限长度薄膜面内输运 | 周期单胞、施加温度梯度 | 面内热导率、厚度方向热流分布 |
| 局域加热与散热 | 给定体热源、热浴和反射边界 | 瞬态温升、热点、热源能量收支 |
| 层状异质材料 | 轴向分层、DMM/AMM/MMM 界面 | 界面温降、跨界面热流 |
| N/R 散射作用比较 | 具备独立 N/U 数据的 Callaway 模型 | 热流演化、与匹配 RTA 模型的差异 |

### 1.2 使用边界

- 公开输入散射模型为 `rta` 和 `callaway`。目录中即使存在全矩阵历史示例，也不意味着 `full_matrix` 是可用的公开输入模式。
- Callaway 是碰撞积分的模型近似，不等同于求解完整 LBTE 碰撞矩阵。
- 热源给定的是进入晶格的功率。程序不自动求解电子输运、电子—声子耦合或真实电路耗散。
- 多材料域目前按 x、y 或 z 方向分层；一般三维任意曲面材料界面不属于该输入工作流。
- 普通表面网格单元使用近似体积；特殊等高挤出几何可启用精确裁切体积，见第 6 章。
- AMM 是标量色散声学通道近似；MMM 的混合系数不是表面粗糙度。
- 示例是起点。运行成功、曲线平滑、能量守恒均不能单独证明物理模型正确或结果已收敛。

### 1.3 推荐阅读路线

首次使用者按第 2、3 章完成安装与第一例，再读第 4、6、7、10、11 章。需要材料准备、器件热源或界面计算时，再阅读第 5、8、9 章。参数表用于查询，不要求初次阅读时逐项记忆。

<a id="chapter-2"></a>
## 2. 安装与编译

### 2.1 目录与命令约定

本手册假定你已有本地 PhonoMC 源码目录，不包含获取源码或提交、上传操作。除特别说明，所有命令均在项目根目录执行。终端中不要输入代码块外的文字，也不需要在命令前加 `$`。

先进入项目目录。当前审阅机器可以使用：

```bash
cd "/Users/sxliu/Library/CloudStorage/OneDrive-个人/Software/GITHUB/PhonoMC"
pwd
ls
```

其他机器应替换为实际目录；路径中含空格时保留双引号。

**截图 S01：项目位置。** 等 `ls` 输出完成，截取完整命令、`pwd` 路径以及包含 `CMakeLists.txt`、`src`、`include`、`material`、`example`、`docs` 的列表。通过标准是当前目录正确，不是只看到一个名为 PhonoMC 的上级文件夹。

### 2.2 依赖

| 依赖 | 用途 | 要求 |
|---|---|---|
| CMake | 配置构建 | 3.20 或更高 |
| C++ 编译器 | 编译求解器 | 支持 C++17 |
| HDF5 | 读取材料数据 | 包含 C++ 接口 |
| OpenMP | 多线程并行 | 可选；未发现时可构建串行版 |
| Qhull 的 `qdelaunay` | 一般表面网格工作流 | 长方体基本算例不需要 |
| Python | 绘图和材料辅助工具 | 本手册建议 3.11 或更高，以直接使用 `tomllib` |
| NumPy、Matplotlib | 后处理 | 使用绘图工具时安装 |
| h5py | 检查 HDF5、生成合成演示材料 | 使用相关工具时安装 |

在 macOS 上，若已经安装 Homebrew，可按需要执行：

```bash
brew install cmake hdf5 qhull libomp
```

若缺少 Apple 命令行编译工具，可以执行 `xcode-select --install`，按系统窗口完成安装后再继续。未安装 Homebrew 时，应先按所用软件包管理器的安装说明配置环境；不要把 `brew: command not found` 当作 PhonoMC 错误。

Ubuntu/Debian 常见安装方式：

```bash
sudo apt update
sudo apt install build-essential cmake libhdf5-dev qhull-bin python3 python3-venv python3-pip
```

GCC 通常配套提供 OpenMP 支持；使用 Clang 且缺少运行库时，可按发行版配置安装 `libomp-dev`。包管理命令会改变系统环境，已有依赖可跳过，学校或集群机器按管理员提供的模块环境操作。

检查基础工具：

```bash
cmake --version
c++ --version
python3 --version
command -v qdelaunay
```

**截图 S02：环境检查。** 截取版本号和 `qdelaunay` 路径。只做长方体案例时允许暂缺 Qhull，并在审阅记录注明。HDF5 是否可供 C++ 构建使用，以后面的 CMake 检测结果为准。

### 2.3 配置独立构建目录

本手册使用 `build-manual`，避免混淆项目中其他用途的构建目录。

```bash
cmake -S . -B build-manual \
  -DCMAKE_BUILD_TYPE=Release \
  -DPHONOMC_ENABLE_OPENMP=ON \
  -DBUILD_TESTING=OFF
```

成功时应看到 HDF5 查找结果以及 `Configuring done`、`Generating done`。启用并行时还应看到 `OpenMP enabled for PhonoMC`。出现“未找到 OpenMP、构建串行版”通常不是构建失败，但后续设置线程数不会使串行可执行文件变成并行程序。

**截图 S03：配置结果。** 保留完整配置命令、HDF5/OpenMP 检测信息和结尾成功信息；一屏放不下可保存 `S03a`、`S03b`。

macOS Apple Clang 未自动找到 OpenMP 时，可以重新配置同一目录，并显式指定运行库：

```bash
cmake -S . -B build-manual \
  -DCMAKE_BUILD_TYPE=Release \
  -DPHONOMC_ENABLE_OPENMP=ON \
  -DBUILD_TESTING=OFF \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib"
```

仅在已安装 Homebrew `libomp` 的情况下使用该命令。更换编译器或 HDF5 安装来源时，选择另一个新构建目录。明确需要串行版时将 `PHONOMC_ENABLE_OPENMP` 设为 `OFF`。

### 2.4 编译

```bash
cmake --build build-manual -j 4
ls -l build-manual/PhonoMC
```

`-j 4` 指编译阶段并行任务数，与模拟阶段的线程数是不同设置。通过标准是构建成功且可执行文件存在。

**截图 S04：编译成功。** 截取构建结尾、`Built target PhonoMC` 或等价成功信息，以及可执行文件列表。不要把只有中间百分比的截图当作完成证据。

### 2.5 配置后处理环境

```bash
python3 -m venv .venv-manual
source .venv-manual/bin/activate
python -m pip install numpy matplotlib h5py
python -c "import numpy, matplotlib, h5py; print('后处理依赖检查通过')"
```

以下 Python 命令默认在该环境中执行。新开终端后需要再次 `source .venv-manual/bin/activate`。求解器本身不依赖这个 Python 虚拟环境。

**截图 S05：后处理环境。** 截取激活命令、依赖导入命令和“后处理依赖检查通过”。安装日志不必逐行截图。

<a id="chapter-3"></a>
## 3. 第一个计算：快速入门

### 3.1 计算目标

用一个 `100 × 100 × 100 nm³` Si 长方体演示跨面热输运。这里的跨面方向是 **x**，不是固定约定为 z：x 两端设置 301 K 和 299 K 热浴，y/z 方向采用周期边界。沿 x 划分 10 个单元。

```text
             沿 x 输运，长度 100 nm
  301 K 热浴 ┃ [1][2][3]……[9][10] ┃ 299 K 热浴
             y、z 两方向周期配对
```

手册配套的 `01_cross_smoke.toml` 使用 10,000 个载体、0.1 ps 时间步、200 步，总模拟时间仅 20 ps，用于检查流程。`02_cross_reference.toml` 保留原始算例的 100,000 个载体和 10,000 步，总时间为 1,000 ps；它也不是预先认证的收敛结果。

### 3.2 阅读完整输入

在编辑器中打开 `docs/user_manual_examples/01_cross_smoke.toml`。该文件是可直接运行的完整输入，核心内容如下：

```toml
[geometry]
model = "box"
sizes = [100, 100, 100]          # nm，顺序为 x/y/z

[simulation]
particle_count = 10000
time_step = 0.1                 # ps
iterations = 200                # 短流程检查，不作为收敛结果
random_seed = 12345
compute_kappa = true
transport_axis = "x"
initial_temperature = "linear"
background_temperature = 300
lifetime_temperature = "local"
grid_xyz = [10, 1, 1]
convergence_write_interval = 10
progress_temperature_summary_only = true

[scattering]
model = "rta"

[boundary]
boundary_position = [
  [-0.01,0,0,0.01,1,1], [0.99,0,0,1.01,1,1],
  [0,-0.01,0,1,0.01,1], [0,0.99,0,1,1.01,1],
  [0,0,-0.01,1,1,0.01], [0,0,0.99,1,1,1.01]
]
boundary_conditions = ["T","T","P","P","P","P"]
boundary_values = [301,299,0,0,0,0]
periodic_pair = [
  [0,-0.01,0,1,0.01,1], [0,0.99,0,1,1.01,1],
  [0,0,-0.01,1,1,0.01], [0,0,0.99,1,1,1.01]
]

[io]
material_folder = "../../material/Si"
output_folder = "../manual_runs/cross_smoke"
```

这里的材料路径可从输入文件所在目录解析，输出路径始终相对于输入文件所在目录解析。因此结果写入 `docs/manual_runs/cross_smoke_序号`，而不是项目根目录下的 `manual_runs`。

**截图 S06：输入上半部分。** 显示文件名、几何、模拟参数和散射模型。检查载体数、步长、步数、网格、输运方向和温度。

**截图 S07：输入下半部分。** 显示三组边界数组、周期配对和 `[io]`。检查六个边界条件一一对应，周期边界值为 0。

### 3.3 执行短流程检查

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/01_cross_smoke.toml
```

不要同时启动多个相同算例；按顺序运行便于记录实际输出目录。初始化时会读取材料、建立几何和缓存，不能仅凭暂时没有步进输出判断死机。

**截图 S08：启动。** 当终端出现 `Results folder:` 和初始化信息时，截取运行命令、结果目录及材料/几何信息。先截图再等待结束；短算例已结束时可滚动回看。

**截图 S09：结束。** 等终端恢复命令提示符，截取最后一段进度与结束信息。若出现异常，保留完整报错，不继续按成功步骤解释结果。物理上应检查热流方向与温差一致，但 20 ps 短算例的瞬时热导率不要求等于任何文献数值。

### 3.4 找到本次结果

将下面的路径替换为 S08 中实际打印的目录。第一次运行可能是 `_0`，再次运行会递增，不能始终复制 `_0`。

```bash
RESULT_DIR="docs/manual_runs/cross_smoke_0"
ls "$RESULT_DIR"
head -n 60 "$RESULT_DIR/summary.txt"
head -n 6 "$RESULT_DIR/grid_centers.csv"
head -n 2 "$RESULT_DIR/convergence.txt"
```

至少应有 `summary.txt`、`grid_centers.csv` 和 `convergence.txt`。10 个网格温度列为 `T_1` 到 `T_10`；网格文件索引则从 0 开始。线性初始化的单元中心温度不应被误认为必须等于边界温度。

**截图 S10：结果与表头。** 可分成 `S10a` 摘要、`S10b` 网格与收敛表头。保留文件名，长表头可在编辑器中打开查看，不必强行缩小到看不清。

### 3.5 绘制结果

确保 Python 环境已激活：

```bash
python tools/plot_convergence.py "$RESULT_DIR" --max-temp-lines 10 --show-legend
ls "$RESULT_DIR/plots_1d"
```

打开 `plots_1d` 中的 `temperature_vs_time.png`、`heatflux_vs_time.png`、`kappa_vs_time.png`。macOS 可以在 Finder 中打开，也可以执行：

```bash
open "$RESULT_DIR/plots_1d/temperature_vs_time.png"
open "$RESULT_DIR/plots_1d/heatflux_vs_time.png"
open "$RESULT_DIR/plots_1d/kappa_vs_time.png"
```

Linux 请用桌面文件管理器或 `xdg-open`。无图形界面的服务器只需生成 PNG，再在有图形界面的本地环境查看。

**截图 S11：绘图命令与生成文件。** 证明图来自本次结果目录。

**截图 S12、S13、S14：温度、热流、热导率图。** 每张保留标题、横纵轴、单位和图例。标注“短流程检查，未验证收敛”。曲线仍漂移或噪声明显时如实保留，不能用截图裁切隐去。

### 3.6 继续到参考算例

完成短例检查后执行：

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/02_cross_reference.toml
```

重复 S08—S14，编号改为 `S15a`、`S15b` 等，记录新的结果目录。比较时间跨度、波动和漂移。只有按第 11 章完成收敛性验证后，才能将热导率作为定量结果报告。

<a id="chapter-4"></a>
## 4. 基本概念与模型选择

### 4.1 声子模式与计算载体

声子模式由 q 点和分支索引共同指定。计算载体携带位置、模式、群速度、占据数及相对于背景的能量偏差；一个载体不是一个真实声子。`particle_count` 控制蒙特卡洛采样规模，而不是体系中的物理声子总数。

RTA/Callaway 碰撞与体热源通过改变现有载体的占据数工作，不通过不断创建载体来实现散射或加热。热浴属于开放边界，可以吸收和注入载体，某些输运路径中的总载体数因此会波动。闭合周期/反射域与开放热浴域应分别检查。

### 4.2 四种容易混淆的温度

| 名称 | 含义 | 配置或来源 |
|---|---|---|
| 初始温度 | 计算开始时的状态 | `initial_temperature` |
| 背景参考温度 | 表示偏差能量的固定基准 | `background_temperature` |
| 散射率查询温度 | 在材料线宽表中查询散射率所用的温度 | `lifetime_temperature`，适用性取决于模型 |
| 重建局域温度 | 由网格能量和模型的能量—温度关系得到的结果 | `T_1` 等输出 |

碰撞算子还使用满足守恒条件的辅助目标温度。它一般不等于背景温度，也不应直接当作输出的局域温度。

普通有限域 RTA 常见配置：

| 方案 | `background_temperature` | `lifetime_temperature` | 含义 |
|---|---:|---|---|
| 全分布采样 | `0` | `"local"` | 相对于零温背景，局域温度查询散射率 |
| 固定散射率偏差计算 | `300` | `300` | 固定参考和固定线宽，仍可使用非线性 Bose 目标 |
| 温度相关偏差计算 | `300` | `"local"` | 固定参考，线宽随重建局域温度变化 |

`background_temperature = "local"` 不合法。固定散射率不等于将整个输运问题线性化。

### 4.3 模型选择表

| 模型 | 适用输入 | 温度处理 | 主要约束 |
|---|---|---|---|
| 普通 RTA | `model = "rta"`，无梯度驱动 | 非线性 Bose 目标，固定或局域查询线宽 | 总线宽数据；注意时间步与材料温区 |
| 梯度驱动 RTA | RTA 加 `[driving]` | 在固定参考温度线性化 | 单材料 box、无热浴和体热源、驱动方向周期 |
| 线性 Callaway | `model = "callaway"`，默认 `linearized` | 固定参考温度的 N/R 双弛豫 | 独立 N/U 数据，小温度扰动；`lifetime_temperature` 不改变该固定参考算子 |
| 非线性 Callaway | `callaway_formulation = "nonlinear"` | 完整普通/漂移 Bose 目标，固定或局域 N/R 率 | 有限域工作流；不能结合周期梯度驱动；温度不得超出 N/R 数据范围 |

非线性 Callaway 当前代码已有输入支持，但不应由此推断任意强非平衡问题都已经得到物理验证。需要额外检查离散模式、守恒求解、采样与时间步。

### 4.4 三种驱动方式

**热浴驱动：**两端温差建立热流。可观察接触温度跳跃，并区分端点温差与内部拟合梯度定义的热导率。

**周期梯度驱动：**给定非零梯度 G，模拟相对于仿射温度背景的周期偏差。输出温度是参考温度加周期修正，不包含完整线性温度斜坡。用于无限重复薄膜面内输运时，x/y 为周期长度，z 为物理厚度。

**局域热源驱动：**输入进入晶格的功率，观察温升和热流。热源开启后自动禁用一般热导率估计。绝热封闭体系持续加热时温度持续升高属于预期行为，不是稳态求解失败。

### 4.5 时间与误差

无提前停止时，总模拟时间为 `time_step × iterations`。输出间隔为 `time_step × convergence_write_interval`。减少时间步时，要相应增加步数，才能比较同一物理时长。

碰撞更新使用指数形式并不意味着整个求解器可任意使用大步长。输运—碰撞分裂仍有误差；热源分为两个半步也不能据此宣称整个算法二阶精度。

<a id="chapter-5"></a>
## 5. 材料数据准备

### 5.1 目录和文件选择

单材料目录例如：

```text
material/Si/
├── POSCAR
└── kappa-fbz.hdf5
```

文件选择优先级为 `kappa-fbz.hdf5`，其次 `kappa.hdf5`；两者都不存在时，只接受目录中唯一且无歧义的其他 `.hdf5`。为了复现，建议一个生产材料目录只放一套明确的输入文件。

POSCAR 提供晶胞体积与倒格矢转换信息，应与 HDF5 所用原胞一致。不可读、不完整或晶格奇异会终止加载，程序不会以任意体积替代。

### 5.2 HDF5 数据约定

其中 `nq` 为完整 q 点数，`nb` 为分支数，`nT` 为温度点数。

| 数据集 | 形状 | 单位与要求 |
|---|---|---|
| `frequency` | `(nq, nb)` | THz；内部转角频率 |
| `qpoint` | `(nq, 3)` | 倒空间分数坐标 |
| `mesh` | `(3)` | 正整数，乘积等于 `nq` |
| `group_velocity` | `(nq, nb, 3)` | Å/ps；不是直接填写 m/s |
| `temperature` | `(nT)` | K，有限、非负、严格递增，至少一个点 |
| `gamma` | `(nT, nq, nb)` | THz 半高半宽线宽，有限、非负 |
| `weight`，可选 | `(nq)` | 若存在，每项必须为 1 |
| `gamma_N`、`gamma_U` | `(nT, nq, nb)` | Callaway 必需，单位与 `gamma` 一致 |
| `gamma_isotope`、`gamma_impurity`、`gamma_defect`，可选 | `(nq, nb)` 或 `(nT, nq, nb)` | 独立附加阻性线宽，与已有项不得重复 |

当前求解器要求完整、等权 q 网格。把不可约网格文件改名成 `kappa-fbz.hdf5` 不会完成展开。Callaway 对网格、第一布里渊区代表点以及 N/U 一致性还有额外检查。

寿命与线宽的关系为 `τ = 1 / (4πγ)`，γ 用 THz 时 τ 用 ps。`gamma = 0` 表示无该通道散射、寿命无穷大；不能用零值掩盖尚未计算的模式。

### 5.3 随附材料

| 材料 | q 网格 | 温度范围 | N/U 数据与用途 |
|---|---|---|---|
| Si | 31×31×31 | 0–1000 K，10 K 间隔 | 教程总线宽，无分离 N/U |
| Ge | 31×31×31 | 0–1000 K，10 K 间隔 | 教程总线宽，无分离 N/U |
| 3C-SiC | 31×31×31 | 100–1000 K，50 K 间隔 | 专用 NEP 势数据，含 N/U 和独立同位素项 |

材料来源和导出检查还应阅读 `material/README.md` 以及 `material/SiC/README.md`、`provenance.json`、`validation.json`。公开材料集合不包含 Graphene。

RTA 对表外温度的线宽查询使用端点值；这只是程序行为，不是可靠的物理外推。线性 Callaway 拒绝表外参考温度；非线性 Callaway 拒绝表外 N/R 查询温度。热力学能量—温度表与线宽表是不同的表，不能因为温度重建成功就认定散射数据覆盖充分。

### 5.4 本地检查材料

下面命令只读取材料文件：

```bash
python - <<'PY'
from pathlib import Path
import h5py
p = Path('material/Si')
preferred = [p / 'kappa-fbz.hdf5', p / 'kappa.hdf5']
path = next((f for f in preferred if f.is_file()), None)
if path is None:
    candidates = list(p.glob('*.hdf5'))
    if len(candidates) != 1:
        raise SystemExit('HDF5 文件缺失或有歧义，请检查材料目录')
    path = candidates[0]
print('材料文件:', path)
with h5py.File(path, 'r') as f:
    for key in ['mesh', 'frequency', 'qpoint', 'group_velocity',
                'temperature', 'gamma', 'gamma_N', 'gamma_U', 'gamma_isotope']:
        print(key, f[key].shape if key in f else '未提供')
    print('温度/K:', f['temperature'][:])
PY
```

**截图 S16：材料检查。** 保留文件名、形状和温度范围。Si 中 `gamma_N`/`gamma_U` 显示“未提供”是预期信息，不能因此直接运行真实 Si Callaway。

### 5.5 加入新材料

1. 用一致的结构和计算设置生成频率、群速度、线宽及需要的 N/U 数据。
2. 检查单位、温度、分支排序和 POSCAR 原胞是否匹配。
3. 若数据是不可约网格，先使用项目转换工具并核查转换结果；查看 `python material/expand_hdf5_to_fbz.py --help` 获取本地脚本的准确参数。
4. 将明确的一套数据放入独立材料目录，在输入中修改材料路径。
5. 先进行等温平衡与小温差测试，再做 q 网格、载体数和时间步收敛。
6. 保存数据来源、结构、势函数或电子结构设置、导出命令、文件校验值与已知限制。

Callaway 要求同一套计算得到的 `gamma_N + gamma_U = gamma`，当前检查的相对容差为 `1e-7`、绝对容差为 `1e-14 THz`。附加同位素、杂质、缺陷项在 RTA 中加到总线宽，在 Callaway 中加到阻性 R 通道；不能把已含这些贡献的总线宽再加一次。

不要在科学计算中启用 `PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1` 来绕过材料读取失败。合成回退仅供测试，不是缺失真实材料数据时的替代方案。

<a id="chapter-6"></a>
## 6. 几何、网格与边界条件

### 6.1 长方体与表面网格

`model = "box"` 时显式提供 `sizes = [Lx,Ly,Lz]`，单位 nm。STL/OBJ 路径则表示从表面网格建立几何；输入坐标同样按 nm 解释。网格加载过程中会平移包围盒，实际坐标以加载摘要和 `grid_centers.csv` 为准。使用非零原点的外部模型时，先核对平移后的热源位置和界面位置，不能直接照搬 CAD 显示坐标。

表面模型应封闭、流形、无零面积三角形，尺度正确。默认 `merge_coplanar_facets = true` 将连通共面三角形合并为面片。边界选择器匹配合并后面片的质心，不是任意落入区域的顶点。

### 6.2 网格

`grid_xyz = [nx,ny,nz]` 划分包围盒。box 保留全部单元；一般表面网格仅保留中心在几何内部的单元。普通网格当前按“域总体积/保留单元数”赋予等体积，并非对每个曲面边缘单元精确求交。

对于适用的等高多边形挤出 STL，可设置 `extruded_z = true`。程序验证垂直侧壁、平顶和平底，使用投影三角形和裁切获得单元体积；一般 FinFET 或复杂曲面不能仅为提高精度就任意启用。

网格越细不必然越准确：固定载体数下，每个单元的样本会减少，空单元和热源待沉积能量可能增加。细化网格通常需要同时提高采样规模。

### 6.3 边界类型

| 代码 | 含义 | `boundary_values` | 检查要点 |
|---|---|---|---|
| `T` | 热浴 | 温度，K，非负 | 热浴放在哪个物理表面 |
| `P` | 周期边界 | 必须为 0 | 每个周期区域恰好配对一次 |
| `R` | 粗糙反射 | RMS 粗糙度，nm，非负 | 0 为名义镜面；离散模式能否表示反射方向仍需检查 |

未选中的表面默认是粗糙度为 0 的反射表面。正式计算建议显式标注重要表面，避免依赖默认行为。

### 6.4 边界区域与配对

每个区域为相对于包围盒的六元数组 `[xmin,ymin,zmin,xmax,ymax,zmax]`，0/1 对应该方向最小/最大坐标。它不是 nm 坐标。

```toml
# 选择 x 最小端面附近的面片质心
boundary_position = [[-0.01, 0, 0, 0.01, 1, 1]]
boundary_conditions = ["T"]
boundary_values = [301]
```

上述片段仅演示一个热浴，不能独立构成双热浴热导率计算。三组数组必须等长，每个声明区域至少匹配一个面片。区域重叠时后声明的条件优先。

`periodic_pair` 填写的是区域坐标数组，连续两个区域构成一对，不是边界编号列表。每个周期区域仅出现一次。程序还检查对面法向、面积、边界边长谱和面片数量，防止把几何不匹配的表面简单连起来。

### 6.5 几何核查流程

先用短算例检查 `summary.txt` 的面片、边界类型、配对数量和包围盒，再看网格中心与三维图。热源位置要与实际保留单元重叠。若粒子频繁恢复、边界追踪失败或一步超过 64 次边界碰撞，应检查网格与时间步，不能仅提高载体数。

**截图 S17：几何与边界摘要。** 标注模型名称、边界方向、网格数和周期配对；FinFET 案例还应保存三维几何及热源位置图，见 S23。

<a id="chapter-7"></a>
## 7. 输入文件参数详解

### 7.1 通用规则

输入文件必须使用 `.toml` 后缀和分节格式。当前解析器拒绝未知节、未知键、重复键和无效组合。字符串加引号，布尔值写 `true`/`false`，数值不加引号。建议使用示例中的简单标量和数组形式，不依赖完整 TOML 标准中的高级语法。

本章表格中的“默认”来自当前配置实现。“应显式填写”表示即使结构体存在内部初值，也不建议依赖它建立真实模型。未启用的可选功能通常省略整节。参数片段应合并到同名节内，不要在完整输入末尾重复添加已有键。

### 7.2 `[geometry]`

| 参数 | 类型、单位 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `model` | 字符串 | `"box"` | `"box"` 或 STL/OBJ 文件路径；不支持旧 `"cylinder"` |
| `sizes` | 3 个正数，nm | box 应显式填写 | x/y/z 长度，例如 `[100,100,100]` |
| `merge_coplanar_facets` | 布尔 | `true` | 合并连通共面三角形；影响边界匹配对象 |
| `extruded_z` | 布尔 | `false` | 等高挤出几何的专用验证和精确裁切体积 |

### 7.3 `[simulation]`

| 参数 | 类型、单位 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `particle_count` | 正数，建议整数 | `10000` | 初始计算载体规模，例如 `100000` |
| `time_step` | 正数，ps | `1.0` | 时间步；正式计算需收敛测试 |
| `iterations` | 正整数 | `10000` | 最大步数 |
| `random_seed` | 无符号 64 位整数 | `12345` | 随机流种子；复现还需固定程序和线程配置 |
| `grid_xyz` | 3 个正整数 | 应显式填写 | `[nx,ny,nz]`，内部初始值为 `[1,1,1]` |
| `initial_temperature` | 非负数 K，或 `"linear"` | `300` | 均匀初温或在有效冷热浴之间插值 |
| `background_temperature` | 有限非负数，K | `300` | 固定偏差背景；Callaway/梯度驱动要求大于 0 |
| `lifetime_temperature` | 非负数 K，或 `"local"` | `"local"` | 普通 RTA/非线性 Callaway 的散射率查询；线性 Callaway 固定在背景温度 |
| `temperature_lookup_dt` | 正数，K | `0.1` | 能量—温度查找表温度间隔，不是时间步 |
| `compute_kappa` | 布尔 | `false` | 输出热导率估计；体热源开启时禁用 |
| `transport_axis` | `"x"`/`"y"`/`"z"` | 未指定则按最大网格维度推断 | 建议总是填写；梯度驱动必须填写 |
| `write_cell_heat_flux` | 布尔 | `false` | 有限域输出逐单元热流；梯度驱动自动输出 |
| `convergence_write_interval` | 正整数，步 | `10` | 诊断输出间隔，不是自动收敛判据 |
| `profile_timers` | 布尔 | `false` | 输出计算阶段耗时 |
| `progress_temperature_summary_only` | 布尔 | `false` | 终端只显示 Tmin/Tavg/Tmax，适合大网格 |
| `conservative_boundaries` | 布尔 | `false` | RTA 选用守恒反射和固定入射模态热浴通量，便于与 Callaway 比较；不会把 Bose 分布线性化 |

### 7.4 `[scattering]`

| 参数 | 类型 | 默认 | 含义与约束 |
|---|---|---|---|
| `model` | 字符串 | `"rta"` | 公开可用值为 `"rta"`、`"callaway"` |
| `callaway_formulation` | 字符串 | `"linearized"` | `"nonlinear"` 仅用于有限域 Callaway，不能结合梯度驱动 |
| `matrix_files` | 路径数组 | 空 | 内部保留键，公开 RTA/Callaway 工作流不要填写；不能借此启用全矩阵求解 |

### 7.5 `[driving]`、`[sampling]`、`[convergence]`

| 完整参数名 | 类型、单位 | 默认 | 含义与约束 |
|---|---|---|---|
| `driving.temperature_gradient` | 有限非零数，K/m | 不启用 | 正负号表示沿显式输运轴的温度梯度方向 |
| `driving.warm_start_ps` | 非负数，ps | `0` | 梯度驱动的体 RTA 初始响应猜测，不是重启文件或真实已运行时间 |
| `sampling.resample_interval` | 非负整数，步 | `0` | 公开现有载体 RTA/Callaway 必须保持 0 |
| `sampling.per_mode_sign` | 整数 1–1024 | `2` | 保留的重采样参数；当前公开路径不启用重采样，通常省略 |
| `convergence.enabled` | 布尔 | `false` | 梯度驱动的经验平稳性提前停止 |
| `convergence.min_time_ps` | 非负数，ps | `2000` | 此时间之前的样本不参与窗口判定 |
| `convergence.window_ps` | 正数，ps | `500` | 每个非重叠统计窗口长度 |
| `convergence.windows` | 整数 4–100 | `4` | 连续参与检查的完整窗口数 |
| `convergence.relative_tolerance` | 正数 | `0.005` | 相对于全局平均热导率的容差 |
| `convergence.absolute_tolerance` | 非负数，W/(m·K) | `0.01` | 绝对容差 |

梯度驱动同时要求：单材料 box；驱动方向的实际面片成对周期；所有边界均无热浴；`initial_temperature = lifetime_temperature = background_temperature > 0`；无体热源、初始模式激发和非线性 Callaway。

提前停止还要求 `compute_kappa = true`，每个窗口至少包含 4 个输出间隔。对全局与所有单元热导率，最近若干窗口均值的极差均需不超过“绝对容差 + 相对容差 × 全局均值绝对值”；全局平均值也须超过绝对容差。最大步数仍有效，必须查看实际停止原因。该检查不提供独立样本置信区间，也不替代时间步、q 网格和载体数收敛。

### 7.6 `[boundary]`

| 参数 | 类型 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `boundary_position` | N×6 数组 | 无显式区域 | 相对包围盒选择器，每项至少匹配一个面片 |
| `boundary_conditions` | N 个字符串 | 与区域配套 | `T`、`P`、`R` |
| `boundary_values` | N 个数 | 与区域配套 | 依次为温度 K、周期占位 0、粗糙度 nm |
| `periodic_pair` | 偶数个六元区域数组 | 无周期时可省略 | 每连续两个区域成对；不是索引数组 |

### 7.7 `[materials]`

| 参数 | 类型、单位 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `folders` | 路径数组 | 多材料时必填 | 至少 2 个，从低坐标到高坐标排序 |
| `interface_axis` | `"x"`/`"y"`/`"z"` | `"x"` | 分层方向 |
| `interface_positions` | 递增数列，nm | 多材料时必填 | 数量为材料数减 1；位于域内、落在网格边界、每层有有效单元 |
| `interface_model` | 字符串 | `"diffuse_mismatch"` | `dmm`/`diffuse_mismatch`、`amm`/`acoustic_mismatch`、`mmm`/`mixed_mismatch` |
| `interface_frequency_bin_thz` | 正数，THz | `0.1` | DMM 的频带采样宽度；AMM 中用于候选搜索索引 |
| `interface_amm_fraction` | 数值 [0,1] | 内部初值 0；MMM 必须显式填写 | MMM 每次事件选择声学通道的概率；AMM 对应 1；DMM 不接受非零值 |
| `interface_parallel_bin_inv_a` | 正数，Å⁻¹ | `0.05` | 声学通道的平行动量搜索索引间隔 |
| `mass_densities_kg_m3` | 正数数组，kg/m³ | 非零声学比例必填 | 每种材料一个密度，顺序与 `folders` 相同 |

只有设置 `folders` 后才能使用其他材料分层键。单材料使用 `[io].material_folder`。热浴面片须垂直于分层轴；分层轴方向不能用周期边界把两种端部材料直接连接。

### 7.8 `[heat_source]`：空间和功率

| 参数 | 类型、单位 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `enabled` | 布尔 | `false`，部分源参数可触发推断 | 建议始终显式填写 |
| `profile` | 字符串 | `"uniform"` | `"uniform"` 或 `"gaussian"`，描述空间分布 |
| `min`、`max` | 各 3 个数，nm | uniform 必需 | 均匀源区域；Gaussian 中可共同提供作为裁切盒 |
| `center` | 3 个数，nm | 单中心 Gaussian 必需 | 绝对空间坐标 |
| `centers` | N×3 数组，nm | 多中心可选 | 与 `center` 二选一，便于避免歧义 |
| `sigma` | 3 个数，nm | Gaussian 必需 | 正值为标准差；非正值表示该轴均匀 |
| `half_width` | 3 个数，nm | 可省略 | 对非正 sigma 的轴，正值可限制中心两侧均匀区半宽；省略则无该局部限制 |
| `power_density` | 正数，W/m³ | 功率密度模式填写 | 均匀区数值或 Gaussian 峰值 |
| `power_densities` | N 个正数，W/m³ | 多中心可选 | 每个中心的峰值；数量与中心数一致；仍保留正 `power_density` 作为配置基础值 |
| `total_power` | 正数，W | 总功率模式填写 | 与 `power_density`/`power_densities` 互斥，多中心共享这一个总功率 |

当前空间 Gaussian 在每个正 sigma 轴的 **±2σ** 外截断，再结合裁切盒、几何内部有效单元进行选取。它不是在无穷空间积分的完整 Gaussian。使用 `total_power` 时会在最终选定网格上重新归一化；使用峰值功率密度时，总功率取决于网格、裁切和几何。

### 7.9 `[heat_source]`：时间包络与频谱

| 参数 | 类型、单位 | 默认/必需 | 含义与约束 |
|---|---|---|---|
| `time_profile` | 字符串 | `"constant"` | `"constant"` 或 `"square"` |
| `start_time` | 数值，ps | `0` | 包络起始时间，通常取非负值 |
| `end_time` | 数值，ps | `-1` | 负值为不设终点；实际结束时间应晚于起始时间 |
| `period` | 正数，ps | square 必需 | 方波周期 |
| `on_duration` | 数值，ps | `-1` | 非负时显式给出开启时长，不大于周期 |
| `duty_cycle` | 数值 [0,1] | `0.5` | 未指定非负 `on_duration` 时使用的占空比 |
| `amplitude` | 非负数 | `1` | 开启状态乘数，同时作用于总功率/功率密度 |
| `spectrum` | 字符串 | `"thermal"` | `thermal`、`weighted`、`gaussian`，与空间 `profile` 独立 |
| `branches` | 不重复的非负整数数组 | 所有分支 | **从 0 开始**，按 HDF5 顺序；不自动判别声学/光学分支 |
| `frequency_min`、`frequency_max` | 非负数，THz | 对应方向不限 | f=ω/(2π)，两者存在时上限不小于下限 |
| `branch_weights` | 非负数数组 | weighted 必需 | 每个材料的每个分支一个相对能量权重，所选模式至少有正权重 |
| `frequency_center` | 非负数，THz | Gaussian 频谱必需 | 模态能量权重的中心频率 |
| `frequency_sigma` | 正数，THz | Gaussian 频谱必需 | 频率标准差，与空间 `sigma` 不同 |

`branch_weights` 只与 `spectrum = "weighted"` 配合；频谱 Gaussian 的中心和宽度只与 `spectrum = "gaussian"` 配合。频率筛选和分支筛选取交集。选不到材料模式会报错；材料存在该模式但某单元暂时没有对应载体时会形成待沉积能量。

### 7.10 `[mode_excitation]`

| 参数 | 类型、单位 | 默认 | 含义与约束 |
|---|---|---|---|
| `enabled` | 布尔 | `false` | 开启一次性的初始模态占据修改 |
| `frequency_min`、`frequency_max` | 非负数，THz | 对应方向不限 | 目标频窗；有上下限时保持有序 |
| `occupation_multiplier` | 正数 | `1` | 初始占据倍数，大于 1 表示过占据 |

该功能用于初始瞬态激发，不是持续热源，不改变模式采样归一化，也不能与周期温度梯度驱动组合。

### 7.11 `[io]` 与路径规则

| 参数 | 类型 | 要求 | 例子 |
|---|---|---|---|
| `material_folder` | 路径字符串 | 单材料时应填写；目录含 POSCAR/HDF5 | `"../../material/Si"` |
| `output_folder` | 路径字符串 | 应填写非空结果基名 | `"../manual_runs/cross_smoke"` |

**输入材料与几何路径：**相对路径优先尝试当前工作目录，其次尝试输入文件所在目录。**输出路径：**相对路径直接以输入文件所在目录为基准。排查路径歧义时使用绝对路径，且不要在不同位置放置同名但内容不同的材料目录。

移动输入文件时，必须重新检查相对路径。程序自动为输出基名增加 `_0`、`_1` 等后缀；已有目录通常不会被覆盖，但本手册不建议多个进程同时使用同一基名创建目录。

<a id="chapter-8"></a>
## 8. 热源与材料界面

### 8.1 总功率和功率密度

对均匀空间源，离散总功率等于所有选中单元的 `功率密度 × 单元体积` 之和，nm³ 转 m³ 的系数是 `10^-27`。固定峰值的 Gaussian 源不能直接用整个包围盒体积乘峰值计算总功率。

若已知进入晶格的总功率，优先使用 `total_power`。程序对选定的离散空间分布归一化，使基准总功率等于输入值，再乘时间包络与 `amplitude`。这提高功率控制的明确性，但不消除空间分布的网格误差。

均匀局域源片段：

```toml
[heat_source]
enabled = true
profile = "uniform"
min = [40, 0, 0]
max = [60, 100, 100]
total_power = 1e-6
spectrum = "thermal"
time_profile = "constant"
start_time = 0
end_time = 50
```

在 0–50 ps 加热、幅值为 1 时，请求的能量为 `1e-6 W × 50e-12 s = 5e-17 J ≈ 312.08 eV`。这是检查单位的解析值；实际沉积能量还需减去待沉积预算。

### 8.2 空间 Gaussian 与频谱 Gaussian

空间 Gaussian 决定在哪里加热，例如 `center = [50,50,50]`、`sigma = [15,0,0]` 表示沿 x 变化、y/z 均匀。频谱 Gaussian 决定符合筛选条件的现有载体如何分配能量，例如中心 10 THz、宽度 2 THz。

```toml
[heat_source]
enabled = true
profile = "gaussian"
center = [50, 50, 50]
sigma = [15, 0, 0]
total_power = 1e-9
spectrum = "gaussian"
frequency_center = 10
frequency_sigma = 2
start_time = 0
end_time = 5
```

频谱权重按每个采样载体计算，重复采到同一模式也重复参与权重。均匀模式采样极限下，频率直方图趋向态密度乘 Gaussian，而不是“每 THz 总能量严格为 Gaussian”。同理，`branch_weights = [1,1,1,0,0,0]` 不表示前三个分支各获得三分之一总功率；分支总量还取决于可用载体数。

### 8.3 方波与多中心源

方波片段可替换已有源的时间设置：

```toml
time_profile = "square"
start_time = 0
end_time = 100
period = 20
on_duration = 5
amplitude = 1
```

每个完整周期开启 5 ps、关闭 15 ps，长时间平均请求功率为峰值基准功率的 1/4。脉冲边缘需要时间分辨率检查。源结束后待沉积预算仍可能在后续有合格载体到达时释放，不能把这种延迟全部解释成材料固有弛豫。

多中心 Gaussian 可以使用 `centers = [[30,50,50],[70,50,50]]` 与共同 `sigma`。若使用 `total_power`，所有中心共享该功率；若使用 `power_densities`，分别指定各中心峰值，并核对离散积分后的总功率。

### 8.4 能量账本

每个单元满足：**累计请求能量 = 累计沉积能量 + 当前待沉积能量**，允许浮点残差。待沉积能量仍保留在原单元，不会转移到其他单元，也不计入已重建的声子温度。

检查 `heat_source_ledger.txt` 和 `hs_pending_total_ev`。若最终待沉积比例显著，应提高空间/模态采样、检查窄频窗与网格，不能把全部请求能量当作已进入晶格。改善该问题后，仍要独立检查空间形状和统计误差。

**截图 S18：热源定义。** 保留空间区域、功率、单位、时间与频谱设置。

**截图 S19：热源账本。** 保留列名和请求/沉积/待沉积量，注明统计时刻。非零待沉积量应在审阅记录中说明。

### 8.5 DMM 双层界面

```toml
[materials]
folders = ["../../material/Si", "../../material/Ge"]
interface_axis = "x"
interface_positions = [50]
interface_model = "dmm"
interface_frequency_bin_thz = 0.1
```

在长 100 nm、沿 x 划分 10 格的几何中，x=50 nm 正好位于单元边界。DMM 按共享频率带和法向相空间通量选择反射/透射模式。频带有限宽度是物理离散近似，必须与 q 网格一起检查；守恒偏差能量并不自动保证任意非参考温度下都严格平衡。

### 8.6 AMM 与 MMM

MMM 片段：

```toml
interface_model = "mmm"
interface_amm_fraction = 0.5
mass_densities_kg_m3 = [2329, 5323]
interface_frequency_bin_thz = 0.1
interface_parallel_bin_inv_a = 0.05
```

该密度数组沿用 Si/Ge 示例值，正式研究须与选用结构和数据保持一致。`p = interface_amm_fraction` 表示事件级混合：概率 p 使用声学通道，概率 1−p 使用 DMM。`p=0` 是 DMM 端点，`p=1` 是声学端点；混合的是散射事件核，不是把两个界面热阻直接线性插值。

当前声学核标识为 `reciprocal_linearized_qcell_overlap_v2`，采用有限 q 单元重叠与局域色散线性化。三个最低分支允许声学透射，光学模式反射，不包含完整偏振转换或完整各向异性弹性边界求解。旧点分箱声学核的 p>0 结果不能直接作为当前实现的基准。

DMM 中频带宽度改变离散采样物理；AMM 中频率与平行动量分箱主要用于搜索加速，应检查结果对搜索索引设置的独立性，并另外验证 q 网格和非参考温度平衡偏差。

### 8.7 界面热阻如何提取

在稳态近似一维输运下，分别在界面两侧选取适合的内部单元拟合温度趋势，外推至同一界面位置求温度跳跃 ΔT，并用实际穿过界面的热流密度 q 计算面积归一化热阻 `R_K = ΔT/q`，单位 m²·K/W。

说明拟合范围、温度跳跃符号、热流方向、统计窗口以及是否包含接触跳跃。仅有整个器件的 `kappa_eff` 不能直接得出界面热阻；远离一维稳态条件时也不能随意使用一个全局热流代替界面局部热流。

<a id="chapter-9"></a>
## 9. 典型算例教程

### 9.1 算例总览与统一操作

所有配套输入都在 `docs/user_manual_examples`。它们完整保留所需的几何、边界和材料配置，不需要把正文的局部片段拼成输入。每次运行都记录 `Results folder:`，再设置 `RESULT_DIR` 绘图。不要修改输入后继续把旧输出当作新结果。

| 输入文件 | 目标 | 最大模拟时间 | 核心输出 |
|---|---|---|---|
| `01_cross_smoke.toml` | 首次短流程检查 | 20 ps | 基础文件与一维曲线 |
| `02_cross_reference.toml` | 跨面 Si 薄膜 | 1000 ps | 温度、热流、两种热导率 |
| `03_inplane_finite.toml` | 有限长度面内薄膜 | 1000 ps | 粗糙表面和端部热浴的联合效应 |
| `04_inplane_gradient.toml` | 无限重复面内薄膜 | 2000 ps | 梯度热导率和局部热流 |
| `05_bilayer_dmm.toml` | Si/Ge DMM 双层 | 1000 ps | 层间温度跳跃 |
| `06_bilayer_mmm.toml` | Si/Ge MMM 双层 | 1000 ps | 混合界面事件与温度分布 |
| `07_total_power.toml` | 定总功率局域脉冲 | 100 ps | 温升与能量账本 |
| `08_gaussian_spectrum.toml` | Gaussian 频谱加热 | 10 ps | 频谱设置及待沉积预算 |
| `09_finfet.toml` | FinFET 局域热源 | 100 ps | 三维温度与切片 |
| `10_callaway_linearized.toml` | 线性 Callaway 合成流程 | 2 ps | 模型加载与短时输出 |
| `11_gradient_stationarity.toml` | 梯度平稳性筛查 | 以输入的步数×步长为上限 | 实际终止时间和窗口诊断 |
| `12_square_source.toml` | 方波热源 | 50 ps | 脉冲输入与温度响应 |
| `13_callaway_nonlinear.toml` | 非线性 Callaway 合成流程 | 2 ps | 固定 N/R 率下的非线性更新 |

这些时长是输入设置，不表示已经到达稳态。AMM/MMM 初始化可能需要建立较大的通道表；FinFET 的材料和几何初始化也可能占用明显时间。

### 9.2 跨面 Si 薄膜

目标、完整短例和操作见第 3 章。参考算例保持 100 nm 输运长度，x 热浴、y/z 周期。正式研究可扫描 x 长度、温差、时间步和采样量，但应一次只改变一个主要因素。

判断顺序：温度剖面是否稳定 → 热流方向是否合理 → 粒子与能量收支是否异常 → 统计窗口是否合适 → 参数收敛。`kappa_int` 与 `kappa_eff` 不一致可能包含接触温跳或非扩散效应，不应强行调整参数使它们相等。

### 9.3 有限长度面内薄膜

本例长 1000 nm、宽 100 nm、厚 100 nm；x 两端 301/299 K 热浴，y 周期，z 两面粗糙度为 1 nm。

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/03_inplane_finite.toml
```

复制实际结果目录后：

```bash
RESULT_DIR="docs/manual_runs/inplane_finite_0"
python tools/plot_convergence.py "$RESULT_DIR" --max-temp-lines 10 --show-legend
```

检查输入中 x 长度和 z 粗糙度，再看温度、热流和热导率。与跨面例比较时，几何、输运长度及边界均改变，不能把全部差异归因于单一表面粗糙度。若研究粗糙度效应，应固定其他条件只改变 z 面粗糙度；研究长度效应则固定厚度与粗糙度扫描 x。

**截图 S20：有限面内例。** 保存几何/边界设置和曲线，图注明长度、厚度、粗糙度与统计时段。

### 9.4 无限重复薄膜的梯度驱动

本例 x/y 周期，z 厚度 100 nm、表面粗糙度 1 nm；沿 x 施加 `G=100000 K/m`。温度和散射率参考都为 300 K，网格为 `[1,1,20]`，用于分辨厚度方向热流。

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/04_inplane_gradient.toml
```

读取实际结果路径后：

```bash
RESULT_DIR="docs/manual_runs/inplane_gradient_0"
python tools/plot_convergence.py "$RESULT_DIR" --max-temp-lines 20
head -n 6 "$RESULT_DIR/grid_heat_flux.csv"
```

正梯度下平均 `qx` 应为负，`kappa_eff = -〈qx〉/G` 应为正。`kappa_int` 为 0 且不适用，不代表材料不导热。厚度方向的 `qx(z)` 应从逐单元热流文件提取，并在已验证稳定的窗口平均。

**截图 S21：梯度例。** 保留 `[driving]`、三种温度设置、周期边界、`kappa_eff` 曲线与局部热流表头。注明温度曲线不含施加的仿射斜坡。

进一步检查：G 改为负值时热流反向；G 的小幅缩放不应改变线性响应热导率；改变 x 周期长度不应改变均匀薄膜的收敛结果。每次独立运行并记录目录。

`11_gradient_stationarity.toml` 开启平稳性停止，保留 `warm_start_ps = 0`，即默认仍从平衡开始；改为正值才开启热启动。热启动改变初始条件，不能将其曲线作为真实开机瞬态；达到平稳性阈值也不能跳过第 11 章的数值收敛检查。

### 9.5 Si/Ge 双层：DMM 与 MMM

两例均沿 x 在 50 nm 处分层，x 总长 100 nm，10 格；热浴位于两端，横向周期。依次运行：

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/05_bilayer_dmm.toml
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/06_bilayer_mmm.toml
```

每条命令完成并记录目录后再运行下一条。MMM 使用 p=0.5 和示例密度；这是教学参数，未针对具体 Si/Ge 界面标定。

使用各自结果目录绘制时间曲线，然后按第 10 章方法得到统计窗口内的温度剖面。检查界面位置、两侧温降和稳态热流，不直接拿时间最后一行给出界面热阻。

**截图 S22：双层例。** 保存层顺序/界面设置、界面诊断和带界面位置标注的温度图。比较时统一温度范围、横轴尺度和统计规则。若要检查纯 AMM，可复制 MMM 输入，改为 `interface_model = "amm"`、移除混合比例键并保留正确密度，再独立命名输出。

### 9.6 总功率与频谱热源

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/07_total_power.toml
```

本例外边界默认反射，0–50 ps 加热，随后继续运行到 100 ps。持续加热阶段没有热浴排热，关闭后能量可在内部重新分布，但全域平均温度不会必然回到 300 K。按第 8 章的解析能量和账本检查单位与沉积完成程度。

频谱演示：

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/08_gaussian_spectrum.toml
```

本例使用空间 Gaussian 与 10±2 THz 权重，0–5 ps 请求加热、运行至 10 ps。不要将设置的频谱中心直接等同于实际采样能量直方图的峰值。若缩小频窗，须检查每格可用载体以及 pending 能量。

截图沿用 S18、S19，再保存温度曲线。方波例可运行 `12_square_source.toml`，其空间网格较细，更要关注每格载体不足和时间采样。

### 9.7 FinFET 局域加热

配套输入读取 `model/finfet.stl`，基于原始教程的 200,000 个载体、0.05 ps 时间步、2000 步和 `[22,36,30]` 网格。底部热浴为 300 K，侧向周期及粗糙区域由输入选择器确定。Gaussian 源叠加空间裁切盒，峰值为 `1e20 W/m³`。

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/09_finfet.toml
```

运行结束，替换实际结果路径：

```bash
RESULT_DIR="docs/manual_runs/finfet_0"
python tools/plot_temperature_3d.py \
  --input docs/user_manual_examples/09_finfet.toml \
  --results "$RESULT_DIR" \
  --tail 50 \
  --x-slice-rel 0.5 \
  --y-slice-rel 0.5
ls "$RESULT_DIR/plots_3d"
```

`--tail 50` 表示对最后 50 条输出记录平均，不是 50 ps。应根据实际输出间隔换算时间范围，并确认这些记录是否属于所需阶段。截面参数为相对包围盒位置。

**截图 S23：FinFET。** 至少保存四类内容：①输入中的模型/边界/热源；②摘要的域尺寸与网格；③三维温度图；④穿过热点的切片及完整色条。初始中部切片未穿过热点时，可改变切片位置，并记录实际参数。

该例是演示局域晶格加热的模型，不包含自动生成的完整晶体管栅极、氧化层和电子耗散求解。温度峰值随热源分布、网格、边界和统计噪声变化，不能仅凭器件外形相似就与实际器件温升定量对比。

### 9.8 线性与非线性 Callaway 合成演示

先生成独立的教学材料目录：

```bash
python tools/create_callaway_demo.py --out docs/manual_runs/callaway_demo_material
```

此目录若已含材料文件，生成脚本会拒绝覆盖。首次成功后直接复用；需要重新生成时使用新目录，并同步修改输入路径。

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/10_callaway_linearized.toml
```

该材料只有 300 K 的人工 N/U 数据，结果只验证流程，不能称为 Si 或任何真实材料的预测。默认线性 Callaway 以背景温度固定 N/R 率和线性化系数。

非线性演示：

```bash
OMP_NUM_THREADS=4 ./build-manual/PhonoMC docs/user_manual_examples/13_callaway_nonlinear.toml
```

该文件显式设定 `lifetime_temperature = 300`。因为合成材料只有一行 300 K 数据，若改为 `"local"` 而局域温度偏离 300 K，会触发表外温度拒绝，不是应该忽略的故障。真实温度相关非线性计算应提供覆盖所有查询温度的 N/U 数据。

**截图 S24：Callaway。** 保存材料生成输出中的 SYNTHETIC/合成属性、模型配置、运行结束与诊断。图注明确写“合成材料流程演示”。

### 9.9 创建自己的输入

选择与目标最接近的配套例，复制到新的文件名，先改输出基名，再依次修改材料、几何、边界、驱动和采样参数。第一轮保持短时长验证配置，随后逐项做收敛扫描。复制到另一个目录时重新检查所有相对路径。

不要把不同物理模型的输入片段直接叠加。例如梯度驱动不能同时保留热浴和体热源；同一个 TOML 节中的键不能重复；Si/Ge 教程材料不能仅把 `rta` 改为 `callaway` 就认为 N/U 数据已经存在。

<a id="chapter-10"></a>
## 10. 输出文件与后处理

### 10.1 文件清单

| 文件/目录 | 内容 | 产生条件 |
|---|---|---|
| `summary.txt` | 实际输入、几何、材料、边界、运行与诊断摘要 | 正常运行流程；异常中断时可能不完整 |
| `grid_centers.csv` | 索引、中心位置和体积 | 几何网格建立后 |
| `convergence.txt` | 时间、温度、热流、热导率、粒子和能量诊断 | 按输出间隔写入 |
| `grid_heat_flux.csv` | 逐单元温度及 x/y/z 热流密度 | 梯度例自动，有限域可显式开启 |
| `heat_source_ledger.txt` | 最新逐单元热源请求、沉积、待沉积账本 | 开启体热源 |
| `plots_1d/` | 一维时间曲线 PNG | 运行一维绘图脚本后 |
| `plots_3d/` | 三维温度与截面 PNG | 运行三维绘图脚本后 |

结果目录自动递增。截图、分析笔记和输入副本必须记录实际后缀，不使用“最后一个目录”代替精确路径。

### 10.2 `grid_centers.csv`

表头为：

```text
index,x_nm,y_nm,z_nm,volume_nm3
```

索引从 0 开始，`convergence.txt` 的 `T_1` 对应索引 0，`T_2` 对应索引 1。非规则几何中有效单元数可能小于 `nx×ny×nz`，不要按完整笛卡尔数组强行 reshape。

### 10.3 `convergence.txt`

第一行以 `#` 开始，列之间为空白。应按列名读取，避免把可变数量温度列后的第几个数硬编码为热流。

| 列名 | 含义、单位 |
|---|---|
| `timestep`、`time_ps` | 已完成步数、模拟时间 ps |
| `T_1 … T_n` | 网格重建温度，K |
| `heatflux` | 输运轴上的体积平均热流密度，W/m² |
| `kappa_int` | 由单元中心温度拟合梯度得到的热导率，W/(m·K) |
| `kappa_eff` | 热浴端点或施加梯度定义的热导率，W/(m·K) |
| `particle_count` | 当前载体数 |
| `absorbed`、`injected`、`recovered`、`net` | 所记录单步的吸收、注入、恢复、净注入计数 |
| `reservoir_absorbed_energy_ev`、`reservoir_injected_energy_ev` | 单步热浴吸收/注入的偏差能量，eV |
| `hs_occupation_updates` | 单步热源占据更新次数，不是新建载体数 |
| `hs_injected_energy_ev` | 单步实际沉积到声子的热源能量，eV |
| `lifetime_energy_residual_ev` | 单步碰撞守恒求解残差，eV |
| `energy_balance_residual_ev` | 计入已记录交换后的单步能量闭合残差，eV |
| `total_thermal_energy_ev` | 域内热能，去除零点能，eV |
| `hs_prescribed_total_ev` | 累计请求的热源能量，eV |
| `hs_pending_total_ev` | 当前待沉积热源能量，eV |
| `hs_ledger_residual_ev` | 请求−沉积−待沉积的账本残差，eV |
| `gradient_occupation_updates`、`gradient_energy_ev` | 梯度例追加的占据更新次数和单步净能量增量 |

**重要：**输出间隔大于 1 步时，单步计数/能量列只记录输出时刻那一步，不是这两个输出时刻之间的累计量。不能把稀疏记录中的 `hs_injected_energy_ev` 直接求和当作整段总能量。累计值应使用摘要或累计账本字段。

### 10.4 两种热导率

热浴驱动中，`kappa_eff` 使用热浴温差和完整域长；`kappa_int` 使用网格中心温度的最小二乘梯度。边界跳跃明显时两者可能不同。均匀温度或极小梯度也会使热导率比值不稳定，不能只检查一个数是否有限。

梯度驱动中，`kappa_eff = -〈q_a〉/G`，`kappa_int` 不适用。局域加热中一般热导率估计被禁用，重点是温升、热流和功率平衡。整个双层器件的有效热导率也不等于任一层的本征热导率。

### 10.5 读取与统计示例

以下片段在已经设置正确 `RESULT_DIR` 的终端执行。时间窗 `500–1000 ps` 仅用于参考跨面例的演示，必须先看曲线确认适用；短流程例不能使用这个窗口。

```bash
python - "$RESULT_DIR" <<'PY'
import sys
from pathlib import Path
import numpy as np
p = Path(sys.argv[1]) / 'convergence.txt'
with p.open() as f:
    names = f.readline().lstrip('#').split()
a = np.loadtxt(p, comments='#', ndmin=2)
col = {name: a[:, i] for i, name in enumerate(names)}
t = col['time_ps']
mask = (t >= 500) & (t <= 1000)
if mask.sum() < 2:
    raise SystemExit('该时间窗内不足两条记录，请选择实际可用的稳定窗口')
print('样本窗口/ps:', t[mask][0], t[mask][-1], '记录数:', mask.sum())
for key in ['heatflux', 'kappa_eff', 'kappa_int']:
    v = col[key][mask]
    print(key, '记录均值=', v.mean(), '记录标准差=', v.std(ddof=1))
print('末次 pending/eV:', col['hs_pending_total_ev'][-1])
PY
```

等间隔样本的记录均值适合快速检查；若有不等间隔记录，应按实际时间做加权积分。这里的标准差不是均值的不确定度，时间相关性意味着不能直接用“标准差/所有记录数平方根”声称置信区间。

生成一维温度剖面时，可在 Python 中按温度列顺序取窗口平均，然后与 `grid_centers.csv` 的相同索引合并。对于 `[nx,1,1]` 的 box，按 `x_nm` 排序画图即可；多维模型需先指定截面或做体积加权分层平均。

### 10.6 逐单元热流与平均

`grid_heat_flux.csv` 的列为：

```text
timestep,time_ps,cell_index,temperature_K,qx_W_m2,qy_W_m2,qz_W_m2
```

通过 `cell_index` 与网格文件连接。不同体积单元的全域平均应使用 `Σ(q_i V_i)/ΣV_i`。普通网格虽可能等体积，使用体积权重仍便于与等高挤出精确裁切网格保持一致。

`heat_source_ledger.txt` 中坐标字段是 `x_A,y_A,z_A`，单位 **Å**；`grid_centers.csv` 的坐标则是 **nm**。绘图前必须按列名换算，不能仅因都表示位置就直接叠加。

### 10.7 图与截图的用途

绘图脚本生成的原始 PNG 用于报告排版；操作截图用于证明命令、路径和参数。两者都应保存。三维图使用相同色标才能比较不同参数下的温度，不能靠自动色标下的“更红”判断温升更高。

截图图注建议包含：材料、模型、几何尺寸、边界、热源/梯度、载体数、时间步、统计窗口、结果目录。图注不要写“已收敛”，除非附有实际收敛证据。

<a id="chapter-11"></a>
## 11. 收敛性与结果验证

### 11.1 三层检查

**配置正确：**材料与几何加载正常，单位正确，边界匹配，热源选中有效单元。**数值稳定与守恒：**无持续异常粒子恢复，能量残差与数值精度一致，温度/热流无无法解释的漂移。**目标量收敛：**载体、时间步、网格和统计时长变化后目标量在所需精度内稳定。

三层均满足后，还要核查物理近似是否适合研究对象，尤其是 Callaway、界面声学近似和材料数据本身。

### 11.2 推荐扫描设计

以下是扫描设计示例，不是普遍适用的最低精度或阈值。先用基准例估算耗时，再选择实际可承担的规模。

| 检查维度 | 示例变化 | 保持不变的量 | 看什么 |
|---|---|---|---|
| 载体数 | N、2N、4N | 几何、时间步、物理时长 | 均值变化、噪声、每格采样和 pending |
| 时间步 | Δt、Δt/2、Δt/4 | 总物理时间与统计窗口 | 分裂误差与边界追踪 |
| 空间网格 | nx 或 nz 逐步加密 | 物理尺寸；适当增加载体 | 温度梯度、热点、局部热流 |
| 总时长 | t、2t | 其他设置 | 长时间尾部漂移 |
| 随机种子 | 至少若干独立种子 | 其他设置与分析规则 | 独立运行之间的离散程度 |
| 材料 q 网格 | 更细且一致生成的数据 | 物理材料和温度 | 模态积分、镜面伙伴、界面通道误差 |
| DMM 频带 | 例如 0.2、0.1、0.05 THz | 材料数据和其他设置 | 频带离散误差及可用出射通道 |
| 梯度 | 改符号、减幅值 | 几何与参考温度 | 线性响应一致性 |

细化空间网格但保持较少载体，可能导致采样问题比原先更严重。对每个单元收敛困难的热源或 Callaway 算例，应同步关注约束求解和模态覆盖。

### 11.3 稳态与统计窗口

先看温度、热流和目标量随时间的变化，舍弃初始化阶段。用多个不重叠窗口检查均值是否仍系统漂移。延长计算检查较慢弛豫过程；不要仅截取看起来最平的一小段。

时间序列通常相关。可以采用足够长的分块均值或独立种子重复运行评估不确定度。报告所用方法、样本数和窗口，不把单条曲线的抖动范围直接写成材料性质误差。

### 11.4 能量与粒子检查

碰撞残差和单步能量闭合应接近相应计算尺度下的浮点误差，观察绝对值及相对于实际能量交换的比例，不使用缺乏量纲说明的统一阈值。有限精度输出也会影响从文本重新计算的残差。

有热源时核对请求/沉积/待沉积总账；有热浴时核对累计吸收与注入。关闭热源的封闭体系不应持续无故增能。开放域载体数波动可以正常，但持续单向漂移、频繁恢复或无法解释的能量变化需要排查。

多材料界面还应做等温平衡测试，尤其选择不同于偏差参考温度的等温值，以暴露有限频带或 q 单元近似的平衡偏差。只在零偏差初态测试可能看不出这些误差。

### 11.5 截图与记录

**截图 S25：收敛比较。** 在统一坐标、色标和统计窗口定义下展示不同 N、Δt、网格、时长或种子的结果。表格至少列出输入名、结果目录、参数、窗口、目标量和待沉积比例。

**截图 S26：最终审阅记录。** 展示可复现配置、验证结论、尚未解决的问题。允许结论为“流程已通过，但载体数收敛未完成”，不能为了完成截图流程而给出没有证据的物理结论。

建议以 `docs/manual_review_notes.md` 保存个人记录，按“截图编号—操作—观察—问题—结论”填写；附录给出可复制的模板。

<a id="chapter-12"></a>
## 12. 性能与故障排查

### 12.1 提高效率

使用 Release 构建；先完成短例再增加规模。确认 OpenMP 在编译时启用后，通过 `OMP_NUM_THREADS=1/2/4/8` 等进行短例计时比较。线程越多不一定越快，材料表、内存带宽、单元统计和几何工作量都会限制扩展性。

在 `[simulation]` 中设 `profile_timers = true`，从摘要观察初始化、粒子推进、碰撞、温度更新、热源和统计等阶段占比。初始化很慢与每个时间步很慢属于不同问题。

粗略扩大 N 或步数可能近似增大计算量，但复杂几何和声学界面通道构建并不保证线性扩展。批量运行时应让“同时运行进程数 × 每进程线程数”与机器实际能力相符，并分别命名输出。

### 12.2 常见问题

| 现象/错误 | 优先检查 | 处理方式 |
|---|---|---|
| `cmake` 或编译器找不到 | 系统依赖与 PATH | 安装/加载相应工具，重新开终端检查 |
| HDF5 C++ 接口缺失 | 是否安装开发库及 C++ 组件 | 更正依赖来源，重新配置新构建目录 |
| OpenMP 未找到 | 运行库、编译器匹配 | 安装运行库；macOS 按第 2 章显式指定；或先用串行版 |
| `BUILD_TESTING` 报缺少测试目录 | 是否在公开目录错误启用本地测试 | 用户运行使用 `-DBUILD_TESTING=OFF` |
| 未知 TOML 键/节 | 拼写、所在节、旧版本输入 | 对照第 7 章；不要复制历史别名或 `full_matrix` 模式 |
| 数字类型错误 | 是否写成 `"300"` 或非法数组 | 使用不带引号数值、合法布尔和整数网格 |
| 边界数组长度不一致 | 三个数组的项数 | 逐项对齐，周期值也保留 0 |
| 区域没有匹配面片 | 相对坐标、合并后质心 | 看摘要与几何，适当扩大边界薄层选择器 |
| 周期配对失败 | 方向、面积、边数、配对重复 | 修正几何与区域，不跳过校验 |
| POSCAR 加载失败 | 路径、原胞、体积 | 使用有效且与 HDF5 匹配的 POSCAR |
| HDF5 形状/权重失败 | 完整网格、数据集、温度维 | 正确展开/导出，不通过改文件名规避 |
| `qdelaunay` 缺失 | Qhull 命令行工具 | 安装并确保 PATH 可见；box 可先验证基础流程 |
| Callaway 缺 N/U | 材料是否提供 `gamma_N`/`gamma_U` | 换用验证过的数据；Si/Ge 教程总线宽不能自动拆分 |
| Callaway 温度越界 | 参考/局域查询与 N/R 温区 | 提供覆盖温区的数据；单温合成演示固定查询温度 |
| 源未选中网格中心 | 源宽度、位置、加载坐标、网格 | 看网格中心并调整分辨率/源定义，不只增大功率 |
| pending 能量很大 | 空单元、窄频窗、分支样本不足 | 增加载体，合理设置网格与目标模式，复查最终账本 |
| 温度持续上升 | 是否为无热浴持续加热 | 先判断物理设置；需要散热稳态时定义实际热沉 |
| 温度或热流噪声大 | 每格样本、统计时长、种子 | 增加采样与合理平均，不用图像平滑代替计算收敛 |
| 粒子频繁恢复/轨迹失败 | 网格封闭性、细长结构、时间步 | 修复几何并减小时间步；适用时使用挤出几何选项 |
| 热导率为零或异常 | 热源禁用、梯度例 `kappa_int`、温差太小 | 先核对输出定义，再判断是否数值错误 |
| 找不到绘图结果 | 实际目录后缀、Python 环境 | 使用打印出的路径，检查脚本错误和依赖 |
| 改种子或线程后结果不同 | 随机流与并行配置 | 以统计一致性比较；逐位复现需固定可执行文件及线程等配置 |

### 12.3 怎样提交一份可诊断的本地问题记录

保存完整输入、终端运行命令、第一条关键错误及其上下文、材料来源、构建方式和实际结果目录。截图要覆盖错误文本，不只截弹窗标题。可以附 `summary.txt` 和相关表头；长文件以文件形式保存，截图作为操作定位证据。

如果本次运行异常中断，不把部分输出当作完整稳态结果。无需重复启动耗时算例来“试试看”，先定位配置、数据或几何错误。

### 12.4 本地测试与文档

普通用户只需构建求解器。维护者若确实拥有本地 `tests/`，可另建目录启用 `BUILD_TESTING=ON`；测试通过证明相关数值契约，不直接证明材料预测精度。不要为用户安装流程依赖不随公开源码提供的测试文件。

原有英文文档使用 Sphinx；本中文手册当前为独立 Markdown 审阅文件，不会自动加入原有 Sphinx 导航。可以直接在支持 Markdown 的编辑器中阅读。是否纳入正式文档站点、发布或上传，属于后续独立安排。

<a id="appendix"></a>
## 附录

### A. 单位速查

| 物理量 | 输入/输出约定 | 常用换算 |
|---|---|---|
| 几何、STL/OBJ、热源位置、粗糙度 | nm | 1 nm = 10 Å = 10⁻⁹ m |
| `grid_centers.csv` 坐标/体积 | nm / nm³ | 1 nm³ = 10⁻²⁷ m³ |
| 热源逐格账本坐标 | Å | 10 Å = 1 nm |
| 时间 | ps | 1 ps = 10⁻¹² s |
| 温度 | K | 温差 1 K 与 1 °C 数值相同，绝对温度不可混用 |
| 总功率 | W | 1 nW = 10⁻⁹ W，1 μW = 10⁻⁶ W |
| 功率密度 | W/m³ | 不要与表面热流密度混淆 |
| 热流密度 | W/m² | 总热流功率为面积积分，不直接等同该列数值 |
| 热导率 | W/(m·K) | 使用一致的 SI 梯度 |
| 界面面积归一化热阻 | m²·K/W | ΔT/热流密度 |
| 温度梯度 | K/m | 1 K/nm = 10⁹ K/m |
| 频率 f、线宽 γ | THz | f=ω/(2π)，不要把角频率直接填入频窗 |
| HDF5 群速度 | Å/ps | 1 Å/ps = 100 m/s |
| 能量账本 | eV | 1 eV = 1.602176634×10⁻¹⁹ J |
| 质量密度 | kg/m³ | 1 g/cm³ = 1000 kg/m³ |

### B. 中英文术语

| 中文 | 英文/缩写 | 本手册含义 |
|---|---|---|
| 声子模式 | phonon mode | q 点与分支的组合 |
| 计算载体 | computational carrier | 具有统计权重的模拟样本 |
| 偏差能量 | deviational energy | 相对于固定平衡参考的能量差，可有正负 |
| 弛豫时间近似 | RTA | 单通道近似碰撞模型 |
| 正常散射 | normal process, N | 保持晶体动量约束的通道 |
| 倒逆散射 | Umklapp process, U | 阻性通道的一部分 |
| 阻性散射 | resistive process, R | U 与独立同位素/杂质/缺陷项的合并 |
| 漫散射失配模型 | DMM | 频带与法向相空间通量加权界面模型 |
| 声学失配近似 | AMM | 当前实现中的标量色散声学通道 |
| 混合失配模型 | MMM | 声学与 DMM 事件核概率混合 |
| 热浴 | thermal reservoir | 吸收与注入声子载体的温度边界 |
| 平稳性 | stationarity | 统计窗口观测量无显著持续漂移 |
| 收敛性 | convergence | 相关数值参数进一步改善后目标量稳定 |
| 待沉积能量 | pending energy | 已请求但尚未分配给合格载体的源预算 |

### C. 主要公式与解释

1. Bose 占据：`n_BE(T,ω)=1/[exp(ℏω/k_B T)−1]`。
2. 单模式偏差能量：`ΔE=ℏω[n−n_BE(T_ref,ω)]`；求宏观能量时还需正确统计权重与晶体归一化。
3. 总 RTA 线宽：`γ_total = γ + γ_isotope + γ_impurity + γ_defect`，仅添加存在且独立的项。
4. 寿命：`τ=1/(4πγ_total)`；有限步松弛因子 `1−exp(−4πγ_total Δt)`。
5. Fourier 定义：`q_a = −κ ∂T/∂a`；梯度驱动 `κ_eff=−〈q_a〉/G`。
6. 恒定功率请求能量：`E=PΔt`，时间换算为秒后得到 J。
7. 离散体积平均：`〈A〉=Σ_i A_i V_i / Σ_i V_i`。
8. 界面面积归一化热阻：`R_K=ΔT_interface/q_normal`，需给出跳跃定义与统计条件。

这些公式用于理解与检查结果；实际碰撞步骤还包含在现有载体样本上求解的守恒约束，不能直接用某一个公式代替完整算法。

### D. 可复现记录模板

```text
案例名称：
审阅日期与手册日期：
源码快照/本地版本标记（无需执行 Git）：
可执行文件路径与校验值：
编译器、构建类型、OpenMP 状态：
运行命令与线程数：
完整输入文件及校验值：
材料目录、数据来源与校验值：
材料温区、q 网格、N/U/附加散射项：
几何尺寸、网格、边界、热源/梯度：
载体数、随机种子、时间步、最大步数：
实际最终时间与停止原因：
实际输出目录：
统计窗口、平均方法、独立种子/分块方法：
目标量与不确定度：
能量、粒子、pending 检查：
载体数/时间步/网格/时长/q 网格收敛证据：
截图编号与原始 PNG 路径：
已知限制与待解决问题：
```

可用 Python 标准库计算某个具体文件的 SHA-256，无需版本控制命令：

```bash
python - <<'PY'
from pathlib import Path
import hashlib
p = Path('docs/user_manual_examples/01_cross_smoke.toml')
print(p, hashlib.sha256(p.read_bytes()).hexdigest())
PY
```

### E. 本手册依据与进一步阅读

本地输入实现优先用于确定当前可接受的键、默认值和组合限制。旧文档局部描述若与后续新增功能冲突，不应将旧段落单独作为当前功能范围。

- [项目概览](../README.md)：功能、构建、示例和公开范围。
- [配置定义](../include/SimulationConfig.h)、[配置解析与校验](../src/SimulationConfig.cpp)：参数表的主要依据。
- [英文输入参考](source/input_files.rst)、[基本方法](source/methodology.rst)：输入及方法解释。
- [Callaway](source/callaway.rst)、[梯度驱动](source/temperature_gradient.rst)：线性/非线性选项与限制。
- [材料说明](source/materials.rst)、[材料目录说明](../material/README.md)：数据格式和来源。
- [几何与边界](source/geometry_boundaries.rst)、[结果分析](source/results_analysis.rst)：网格和输出说明。
- [输出格式实现](../src/solver/ResultWriter.cpp)、[网格输出实现](../src/SimulationDomain.cpp)：列名、单位和索引。
- [热源实现](../src/solver/PrescribedHeatSource.cpp)：空间截断、离散功率归一化及账本。
- [参考文献目录](source/reference.rst)：研究论文与理论出处。引用论文时还应核对原文及所用模型，不把本手册当作独立物理验证。

### F. 审阅版状态

本手册与配套输入用于本地审阅。输入校验范围与文档检查结果见 [检查记录](user_manual_examples/VALIDATION.md)。未展示的实际操作截图需由读者逐步补充；配套教学参数未被宣称为已收敛物理结果。编写工作不包含提交、推送、发布或上传操作。
