# MagPackage: 高性能三维低频磁场有限元计算包

**MagPackage** 是一个基于 MATLAB 开发的高性能三维电磁场有限元仿真框架。它专为低频磁场问题设计，采用 **Nedelec (Edge) 单元** 求解磁矢位 $\mathbf{A}$ 公式，并结合 **C++/MEX** 加速内核与 **MUMPS** 直接求解器，实现了接近原生 C++ 代码的计算效率，同时保留了 MATLAB 的易用性。

## 🌟 核心特性 (Features)

### 1. 多物理场求解器

- **静磁场 (Magnetostatic)**: 线性/非线性各向同性材料求解。
- **频域涡流场 (Frequency Domain / Eddy Current)**: 求解时谐磁场，支持趋肤效应分析。
- **瞬态非线性场 (Transient Nonlinear)**: 基于 Implicit Euler 和 Newton-Raphson 的时域步进求解，支持非线性 B-H 曲线。
- **场路耦合 (Field-Circuit Coupled)**: 支持电压驱动线圈，求解包含外电路（电阻、电感）的耦合系统。
- **谐波平衡法 (HBFEM)**: 针对周期性非线性问题的高效频域求解器，支持 DC 偏置及高次谐波分析。
- **标量场 (Scalar Potential)**: 求解电流传导、静电场或热传导问题 ($V$-formulation)。

### 2. 高性能计算核心

- **混合内核架构**:
  - **C++/MEX 加速**: 核心组装函数（刚度、质量、雅可比矩阵）使用 C++ OpenMP 并行实现，速度提升 10-50 倍。
  - **MATLAB 并行回退**: 自动检测环境，若未编译 MEX 则自动回退到 MATLAB `parfor` 实现，保证兼容性。
- **高级线性求解器**: 集成 **MUMPS** (Multifrontal Massively Parallel Sparse direct Solver) 接口，支持复用符号分析 (`ReuseAnalysis`)，大幅加速非线性迭代和瞬态计算。

### 3. 先进功能

- **材料库**: 支持 B-H 样条插值、导数平滑处理及深度饱和区的物理修正防御机制。
- **线圈建模**: 包含 `CoilGeometryUtils`，支持任意形状线圈（跑道型、圆角矩形）的几何特征自动识别与源项加载。
- **后处理**: 内置损耗计算器（欧姆损耗、Steinmetz 铁损）、电感参数提取及场可视化工具。

## 📂 目录结构 (Directory Structure)

```
magpackage/
├── src/                # 源代码核心
│   ├── Assembler.m     # 总装配器 (核心调度类)
│   ├── Mesh.m          # 网格对象
│   ├── DofHandler.m    # 自由度管理
│   ├── kernels/        # MATLAB版参考内核
│   ├── mex/            # C++ 高性能内核源码 (.cpp)
│   ├── mumps/          # MUMPS 求解器接口
│   ├── solver/         # 各类物理场求解器 (Frequency, Transient, HBFEM...)
│   └── ...
├── tests/              # 单元测试与基准案例
├── tutorials/          # 教程 (如 TEAM7 Benchmark)
├── install.m           # 安装脚本 (添加路径)
├── make.m              # 编译脚本 (构建 MEX)
├── package_project.m   # 自动打包工具
└── run_all_tests.m     # 自动化测试运行器
```

## 🛠️ 安装与配置 (Installation)

### 1. 环境要求

- MATLAB R2020b 或更高版本。
- C++ 编译器 (用于编译 MEX 文件)。
  - Windows: Visual Studio 2019/2022 (推荐) 或 MinGW-w64。
  - Linux/Mac: GCC 或 Clang。

### 2. 安装步骤

1. **下载源码** 并解压到任意目录。

2. 启动 MATLAB，切换到 `magpackage` 根目录。

3. 运行安装脚本以配置路径：

   ```
   install
   ```

4. (关键) 编译高性能内核：

   运行 make.m 自动编译 C++ MEX 文件。

   ```
   make
   ```

   *成功编译后，系统将自动使用 `src/mex/` 下的高速内核。*

## 🚀 快速开始 (Quick Start)

以下代码展示了如何求解一个简单的 TEAM 7 问题（线圈+铝板涡流）：

```
% 1. 加载网格
mesh = Mesh.load('tutorials/TEAM7/Team7.mphtxt');

% 2. 定义物理空间 (Nedelec Edge Elements)
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

% 3. 定义材料与物理参数
MatLibData(1) = MaterialLib.createLinear(1.0); % 空气
MatLibData(2) = MaterialLib.createLinear(1.0); % 铝板
SigmaMap(2) = 3.526e7; % 铝电导率

% 4. 配置求解器 (50Hz 频域)
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, 50);

% 5. 求解
% (SourceMap 定义略，详见 tutorials/TEAM7/run_team7_v3.m)
A_sol = solver.solve(space, MatLibData, SigmaMap, SourceMap, fixedDofs);

% 6. 后处理与绘图
viz = Visualizer(PostProcessor(assembler));
viz.plotFieldMagnitude(A_sol, space, 'B');
```

更多完整案例请参考 `tutorials/` 和 `tests/` 目录。

## ✅ 自动化测试 (Testing)

项目包含一套完整的自动化测试套件，覆盖了从底层内核到高层求解器的所有功能。

运行所有测试：

```
run_all_tests
```

脚本将显示图形化进度条，并生成通过率报告：

```
############################################################
                     TEST SUMMARY                           
############################################################
Total Time: 25.42 seconds
Total Tests: 23
  [PASS]: 23
  [WARN]: 0
  [FAIL]: 0
############################################################
```

## 📚 文档 (Documentation)

详细的设计文档和技术说明位于 `docs/` 目录下：

- [HBFEM 求解器设计文档](https://www.google.com/search?q=docs/20251127-1712-HBFEM.md)
- [源项处理模块文档](https://www.google.com/search?q=docs/20251128-0445-源项处理模块技术文档%20(Source%20Processing%20Module%20Documentation).md)
- [TEAM7 基准测试报告](https://www.google.com/search?q=docs/20251128-1152-TEAM7.md)

## 📦 打包与分发

如果您需要将项目打包分发给其他用户，可以使用内置的打包工具：

```
package_project
```