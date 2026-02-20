# EPMC (Energy Phonon Monte Carlo) 仿真框架

EPMC 是一个基于 C++ 开发的、利用声子蒙特卡洛（Phonon Monte Carlo）方法模拟半导体材料（如 Si, Ge, SiN）中声子输运过程的仿真框架。该工具支持复杂几何结构，并能通过 HDF5 文件加载真实的声子色散和散射数据，用于预测微纳尺度的热流分布和等效热导率。

## 核心功能

* **蒙特卡洛粒子追踪**：模拟海量声子粒子的运动、碰撞、边界反射及透射过程。
* **复杂几何支持**：集成 `Qhull` 库实现四面体网格划分，支持快速的点定位与边界判定。
* **真实声子色散**：支持从 HDF5 文件加载完整的声子模式、群速度和寿命数据。
* **高性能并行**：基于 **OpenMP** 的多线程并行架构，适配多核 CPU 和超算节点。
* **统计输出**：实时计算温度梯度、热流以及基于线性拟合的有效热导率。

---

## 环境要求

在 Linux (Ubuntu/Debian) 环境下，请确保安装以下依赖：

| 依赖项 | 推荐版本 | 安装命令 (Ubuntu) |
| :--- | :--- | :--- |
| **编译器** | GCC 9+ / Clang 10+ | `sudo apt install build-essential` |
| **CMake** | 3.20+ | `sudo apt install cmake` |
| **HDF5** | 1.10+ | `sudo apt install libhdf5-dev` |
| **Qhull** | 2020+ | `sudo apt install qhull-bin` |
| **OpenMP** | - | `sudo apt install libomp-dev` |

---

## 编译指南

执行以下步骤进行本地编译：

1. **创建构建目录**：
   ```bash
   mkdir build && cd build
   ```

2. **配置 CMake**：
   ```bash
   cmake ..
   ```

3. **编译项目**：
   ```bash
   make -j
   ```

编译完成后，可执行文件 `epmc` 将生成在 `build` 文件夹中。

## 使用方法

1. **运行模拟**

   在项目根目录下运行，并指定配置文件（TOML 格式）：
   ```bash
   ./build/epmc input_cross_100nm.toml
   ```

2. **手动指定并行线程数**
   
   如果需要使用 64 核，可以运行：
   ```bash
   OMP_NUM_THREADS=64 ./build/epmc input_cross_100nm.toml
   ```

## 输出文件说明

模拟开始后，结果将保存在 results/ 目录下以项目名命名的文件夹中：

  `convergence.txt`: 核心数据文件，包含步数、各子体积温度、实时热流及热导率预测值。

  `geometry_summary.txt`: 输出当前几何体的体积、表面积及网格划分统计。

  `subvolume_centers.csv`: 用于后处理，记录各测温区间的中心位置。
