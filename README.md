# MagPackage: 高性能三维低频磁场有限元仿真代码包

**MagPackage** 是一套基于 MATLAB 开发的高性能三维有限元（FEM）代码库，专用于低频电磁场仿真。该项目采用面向对象编程（OOP）架构，利用 C++ MEX 混合编程加速组装，实现了全流程并行化计算，并集成 MUMPS 高性能稀疏矩阵直接求解器，实现了在 MATLAB 环境下对复杂三维电磁问题的高效求解。

## 🌟 核心特性 (Key Features)

- **多物理场求解器**:
  - **静磁场/稳态场**: 支持非线性材料的标量/矢量位求解。
  - **频域 (Frequency Domain)**: 解决涡流问题，支持趋肤效应分析。
  - **时域瞬态 (Transient)**: 基于时间步进法的动态磁场仿真。
  - **谐波平衡法 (HBFEM)**: 针对周期性非线性问题的高效频域求解算法。
  - **场路耦合 (Field-Circuit Coupling)**: 支持电压源/电流源激励下的绞线圈与有限元区域的强耦合仿真。
- **高性能计算**:
  - **C++ MEX 加速**: 刚度矩阵、质量矩阵及雅可比矩阵的组装核心采用 C++ 编写，通过 MEX 接口调用，显著提升计算速度。
  - **MUMPS 求解器集成**: 内置 MUMPS (MUltifrontal Massively Parallel sparse direct Solver) 接口，高效处理大规模稀疏线性方程组。
- **先进的单元与材料库**:
  - 支持四面体网格（Tetrahedral Mesh）。
  - **基函数**: 采用 Nedelec 棱边元（一阶）处理矢量场，Lagrange 节点元处理标量场。
  - **材料模型**: 内置丰富的 B-H 曲线库，支持非线性铁磁材料建模。
- **网格接口**:
  - 支持 **Gmsh** (`.msh`) 网格文件导入。
  - 支持 **COMSOL** (`.mphtxt`) 网格文件导入。
- **完备的后处理**:
  - 损耗计算（铁耗、铜耗）。
  - 电感、电阻等电路参数提取。
  - 场量可视化与能量分析。

## 🛠️ 安装与环境配置 (Installation)

### 环境要求

- **MATLAB**: 推荐 R2020b 或更高版本。
- **C++ 编译器**: 需安装受 MATLAB 支持的 C++ 编译器（如 MinGW-w64, Visual Studio 等），用于编译 MEX 文件。

### 安装步骤

1. 下载代码库

   将本代码包下载至本地目录。

2. 配置编译器

   在 MATLAB 命令行窗口中运行以下命令以设置 C++ 编译器：

   ```
   mex -setup cpp
   ```

3. 编译内核

   运行项目根目录下的 make.m 脚本，自动编译所有 C++ MEX 加速内核：

   ```
   make
   ```

   *注意：确保编译过程中无报错信息。*

4. 安装路径

   运行 install.m 脚本将项目路径添加到 MATLAB 搜索路径中：

   ```
   install
   ```

5. 验证安装

   运行测试套件以确保所有模块正常工作：

   ```
   run_all_tests
   ```

## 🚀 快速开始 (Quick Start)

### 运行 TEAM7 基准测试

本项目包含了经典的 TEAM Problem 7（非对称导体板涡流问题）基准测试案例。

```
% 进入教程目录
cd tutorials/TEAM7

% 运行仿真脚本
run_team7
```

该脚本将执行以下操作：

1. 加载 `Team7.mphtxt` 网格文件。
2. 配置线圈激励和铝板材料参数。
3. 构建频域耦合求解器。
4. 求解并进行后处理（计算阻抗、损耗等）。

### 基础代码示例

以下是一个构建简单频域求解流程的代码片段：

```
% 1. 读取网格
mesh = Mesh('tests/example_iron_coil.msh');

% 2. 配置仿真参数
femConfig = FemConfig();
femConfig.Frequency = 50; % 50 Hz
femConfig.SolverType = 'FrequencyCoupled';

% 3. 定义材料与物理属性
% (此处需根据具体 Mesh 的 ID 设置材料和边界条件)
% ...

% 4. 实例化求解器并求解
solver = FrequencyCoupledSolver(mesh, femConfig);
solver.solve();

% 5. 后处理与可视化
viz = Visualizer(mesh, solver.Solution);
viz.plotField('B', 'mag'); % 绘制磁通密度模值
```

## 🧪 自动化测试 (Automated Testing)

为了保证代码的稳定性和正确性，MagPackage 内置了一套自动化测试框架。

- 运行所有测试:

  在 MATLAB 命令行中运行 run_all_tests.m 脚本：

  ```
  run_all_tests
  ```

  该脚本会自动扫描 `tests/` 目录下的所有测试文件（如 `test_*.m`），依次执行并报告通过/失败状态。

- 测试覆盖范围:

  测试用例涵盖了基础单元测试（矩阵组装验证）、物理场求解测试（静场、频域、瞬态）以及复杂耦合系统测试（HBFEM、场路耦合）。

## 📦 打包与分发 (Packaging & Distribution)

MagPackage 提供了便捷的脚本用于项目打包，方便分发给其他用户或在不同机器上部署。

- 执行打包:

  运行根目录下的 package_project.m：

  ```
  package_project
  ```

  该脚本会执行以下操作：

  1. 清理临时的编译文件（`.mex`, `.o`, `.obj`）。
  2. 重新编译所有内核以确保兼容性。
  3. 将源代码、文档及必要的网格文件打包成 `.zip` 归档文件或 MATLAB 工具箱安装包 (`.mltbx`)。

- 查看项目树:

  在打包前，可运行 print_project_tree.m 查看当前项目的文件结构，确认没有遗漏重要文件。

## 📂 项目结构 (Project Structure)

```
MagPackage/
├── bin/                  # MEX编译二进制文件
├── docs/                 # 说明文档目录
├── src/                  # 源代码核心目录
│   ├── Assembler.m       # 矩阵组装器类
│   ├── Mesh.m            # 网格处理类
│   ├── solver/           # 各类求解器核心逻辑
│   ├── mex/              # C++ 源代码 (高性能组装内核)
│   ├── kernels/          # MATLAB版本的组装内核 (参考/调试用)
│   ├── mumps/            # MUMPS 线性求解器接口
│   ├── mesh_io/          # 网格读取接口 (Gmsh, COMSOL)
│   ├── post/             # 后处理 (可视化, 参数提取, 损耗计算)
│   └── utils/            # 通用工具函数
├── tests/                # 单元测试与功能测试用例
├── tutorials/            # 教程与基准案例 (如 TEAM7)
├── make.m                # MEX 编译脚本
├── install.m             # 安装/路径配置脚本
├── package_project.m     # 项目打包脚本
├── run_all_tests.m       # 自动化测试入口
└── README.md             # 项目说明文档
```

## 🏗️ 代码包架构 (Software Architecture)

MagPackage 采用分层设计的面向对象架构，旨在平衡 MATLAB 的易用性与 C++ 的计算效率。

### 1. 面向对象设计 (OOP Design)

代码逻辑被封装在几个核心类中，通过组合与继承实现灵活的物理场定义：

- **数据层 (Data Layer)**:
  - `Mesh`: 管理网格拓扑（节点、坐标、单元连接关系）。
  - `MaterialLib`: 封装材料属性（B-H 曲线、导电率），提供统一查询接口。
- **离散层 (Discretization Layer)**:
  - `FunctionSpace`: 定义有限元基函数（Nedelec 棱边元或 Lagrange 节点元）。
  - `DofHandler`: 负责局部自由度到全局自由度的映射，处理悬挂节点和编号优化。
- **组装层 (Assembly Layer)**:
  - `Assembler`: 核心引擎，连接网格、材料与离散空间。它根据 `FemConfig` 调用底层内核。
  - **Kernels**: 具体的积分计算逻辑（如刚度矩阵 $K$、质量矩阵 $M$）被隔离在独立的内核函数中。
- **求解层 (Solver Layer)**:
  - `LinearSolver` / `NonlinearSolver`: 定义牛顿迭代、残差计算等通用流程。
  - `TransientSolver` / `FrequencySolver`: 继承自基础求解器，实现特定物理方程的时间步进或频域迭代。
- **物理组件 (Physics Components)**:
  - `Winding`: 专门处理绞线圈的几何与电路属性（匝数、截面、电阻），实现场路耦合接口。

### 2. 并行与高性能策略 (Parallelism & Performance)

为了克服 MATLAB 在循环处理上的性能瓶颈，MagPackage 采用了混合编程模型：

- **MATLAB 代码层面的切片化与并行 (MATLAB-level Slicing & Parallelization)**:
  - **数据切片 (Data Slicing)**: 在数据预处理阶段，利用 MATLAB 强大的矩阵索引功能，将大规模网格数据进行“切片”处理。通过将节点坐标、单元索引等数据以大块（Block）形式批量提取并传递给底层，避免了逐个单元访问的循环开销，极大提高了数据吞吐量。
  - **任务并行 (Task Parallelism)**: 代码架构设计兼容 MATLAB Parallel Computing Toolbox。支持在上层应用逻辑（如多频率扫描、敏感度分析或独立时间步的后处理）中使用 `parfor` 进行粗粒度的并行加速。
- **C++ MEX 组装加速**:
  - 有限元最耗时的步骤是刚度矩阵的元素级组装（Element-wise Assembly）。本项目将这一过程下沉至 C++ 层（`src/mex/`）。
  - 通过 MATLAB 的 MEX 接口，将大数组指针直接传递给 C++，避免内存拷贝。
  - C++ 内核设计为无状态函数，天然支持 **OpenMP 多线程并行**，可同时计算多个单元的局部矩阵。
- **高性能直接求解器 (MUMPS)**:
  - 对于组装好的稀疏线性方程组 $Ax=b$，调用 **MUMPS** 求解器。
  - MUMPS 利用多核 CPU 进行并行因子分解（Factorization）和回代（Substitution），在处理数百万自由度的大规模问题时，速度远超 MATLAB 内置的 `mldivide` (`\`)。
- **矢量化与稀疏矩阵优化**:
  - 在 MATLAB 层面，尽量使用向量化操作处理后处理和数据预处理。
  - 全程采用 CSR/CSC 稀疏存储格式，最小化内存占用。

## 🧠 求解器与算法说明 (Solvers & Algorithms)

MagPackage 在 `src/solver/` 目录下提供了多种针对不同物理场景的求解器。以下是各求解器的详细算法说明：

### 1. ScalarSolver (标量场求解器)

- **适用场景**: 静电场、直流传导场、标量磁位静磁场。

- 物理方程: 求解广义拉普拉斯/泊松方程：

  $$ -\nabla \cdot (\sigma \nabla u) = f $$

  其中 $u$ 为标量势（如电位或磁位），$\sigma$ 为电导率或导磁率。

- **算法**: 采用 **Lagrange P1 (节点元)** 进行离散，组装对称正定矩阵并求解。

### 2. FrequencySolver (频域涡流求解器)

- **适用场景**: 交流激励下的线性/非线性涡流问题（固定频率）。

- 物理方程: $A-\phi$ 方法或纯 $A$ 方法：

  $$ \nabla \times (\nu \nabla \times \mathbf{A}) + j\omega\sigma \mathbf{A} = \mathbf{J}_s $$

- **算法**:

  - **基函数**: 磁矢位 $\mathbf{A}$ 采用 **Nedelec Edge Elements (棱边元)**，自然满足散度条件。
  - **非线性处理**: 对于非线性材料，采用 Newton-Raphson 迭代法，更新雅可比矩阵。

### 3. FrequencyCoupledSolver (频域场路耦合求解器)

- **适用场景**: 需要考虑线圈外电路（电压源激励）的涡流问题。

- 算法: 在 FrequencySolver 的基础上，引入电路方程作为约束：

  $$ U = I R_{dc} + j\omega \Psi(\mathbf{A}) $$

  构建 扩展的系统矩阵（Stranded Coil Formulation），同时求解磁矢位 $\mathbf{A}$ 和线圈电流 $I$。

### 4. TransientSolver (时域瞬态求解器)

- **适用场景**: 任意波形激励、瞬态过程、非线性材料动态响应。

- 物理方程:

  $$ \nabla \times (\nu(\mathbf{B}) \nabla \times \mathbf{A}) + \sigma \frac{\partial \mathbf{A}}{\partial t} = \mathbf{J}_s(t) $$

- **算法**:

  - **时间离散**: 采用 **Backward Euler (后向欧拉)** 或 **Crank-Nicolson** 差分格式。
  - **非线性迭代**: 每个时间步内进行 Newton-Raphson 迭代，直到残差收敛。

### 5. TransientCoupledSolver (时域场路耦合求解器)

- **适用场景**: 瞬态电压源激励、外部电路包含电容/电感等动态元件。

- 算法: 将有限元方程与电路微分方程联立：

  $$ V(t) = i(t)R + L_{ext}\frac{di}{dt} + \frac{d\Psi}{dt} $$

  形成全局耦合的大型系统矩阵，在每个时间步同步求解场量和电路变量。

### 6. HBFEMSolver (谐波平衡求解器)

- **适用场景**: 周期性激励下的稳态非线性磁场（如变压器空载、铁芯饱和谐波分析）。相比时域方法，能直接获得稳态解，无需经历漫长的暂态过程。
- **算法**: **Harmonic Balance FEM (HBFEM)**。
  - 将 $\mathbf{A}$ 和磁阻率 $\nu(B)$ 展开为傅里叶级数。
  - 在频域内构建耦合矩阵，通过 DFT/IDFT 在频域和时域间转换以计算非线性材料属性。
  - 同时求解基波及各次谐波分量。

### 7. HBFEMCoupledSolver (谐波平衡场路耦合求解器)

- **适用场景**: 电压源激励下的周期性非线性问题（输入电压含谐波或负载非线性）。
- **算法**: 将 HBFEM 与外电路的频域方程（阻抗矩阵）耦合，实现包含高次谐波的场路系统直接求解。

## ⚠️ 注意事项

- **MUMPS 依赖**: 本项目依赖 MUMPS 库进行大规模矩阵求解。请确保您的系统环境支持相关 MEX 文件的调用。
- **网格 ID**: 设置材料属性和边界条件时，请务必检查网格文件（`.msh` 或 `.mphtxt`）中定义的 Physical Domain ID 和 Boundary ID。

## 🤝 贡献与反馈

如果您在使用过程中遇到问题或有改进建议，欢迎提交 Issue 或 Pull Request。

## 📄 License

[在此处添加您的许可证，例如 MIT, GPL 等]