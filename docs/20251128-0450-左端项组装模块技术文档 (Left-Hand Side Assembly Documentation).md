# 左端项组装模块技术文档 (Left-Hand Side Assembly Documentation)

[toc]

## 1. 概述 (Overview)

本模块负责有限元求解器中**左端项 (Left-Hand Side, LHS)** 矩阵的组装。这包括刚度矩阵、质量矩阵、耦合矩阵以及边界条件的施加。这些矩阵共同构成了离散化后的线性代数方程组 $\mathbf{A}\mathbf{x} = \mathbf{b}$ 中的系数矩阵 $\mathbf{A}$。

本模块采用了**高性能计算 (HPC)** 策略，特别是利用 MATLAB 的 `parfor` 进行并行组装，并结合了稀疏矩阵技术以优化内存和速度。

## 2. 理论基础 (Theoretical Basis)

### 2.1 有限元离散化 (FEM Discretization)

在低频电磁场中，我们通常求解关于磁矢位 $\mathbf{A}$ 和电标量位 $V$ 的方程。将连续场变量用基函数展开：

$$\mathbf{A} \approx \sum_{j=1}^{N_e} a_j \mathbf{N}_j, \quad V \approx \sum_{k=1}^{N_n} v_k \phi_k$$

其中 $\mathbf{N}_j$ 是棱边元 (Edge Element) 基函数，$\phi_k$ 是节点元 (Nodal Element) 基函数。

### 2.2 矩阵定义

根据伽辽金法 (Galerkin Method)，我们可以得到以下几个核心矩阵：

1. 磁场刚度矩阵 (Magnetic Stiffness Matrix, $\mathbf{K}$):

   来源于 $\nabla \times (\nu \nabla \times \mathbf{A})$ 项。

   $$K_{ij} = \int_\Omega (\nabla \times \mathbf{N}_i) \cdot \nu \cdot (\nabla \times \mathbf{N}_j) \, d\Omega$$

2. 质量矩阵 (Mass Matrix, $\mathbf{M}$):

   来源于 $\sigma \frac{\partial \mathbf{A}}{\partial t}$ 或 $j\omega\sigma \mathbf{A}$ 项。

   $$M_{ij} = \int_\Omega \mathbf{N}_i \cdot \sigma \cdot \mathbf{N}_j \, d\Omega$$

3. 标量拉普拉斯矩阵 (Scalar Laplacian Matrix, $\mathbf{G}$):

   来源于 $\nabla \cdot (\sigma \nabla V)$ 项。

   $$G_{ij} = \int_\Omega \nabla \phi_i \cdot \sigma \cdot \nabla \phi_j \, d\Omega$$

4. 耦合矩阵 (Coupling Matrix, $\mathbf{C}$):

   来源于 $\sigma \nabla V$ 对 $\mathbf{A}$ 方程的贡献，或 $\nabla \cdot (\sigma \mathbf{A})$ 对 $V$ 方程的贡献。

   $$C_{ij} = \int_\Omega \mathbf{N}_i \cdot \sigma \cdot \nabla \phi_j \, d\Omega$$

## 3. 模块实现详解 (Implementation Details)

### 3.1 磁场刚度矩阵组装 (`assemble_magnetic_stiffness.m`)

**功能**: 计算 $\mathbf{K}$ 矩阵。支持各向同性及非线性磁阻率 $\nu$。

**算法流程**:

1. **输入解析**: 接收模型数据 `Model` 和可选的磁阻率向量 `Nu_Override` (用于非线性迭代更新)。如果未提供 `Nu`，则从材料库计算初始线性值。
2. **数据分块 (Pre-chunking)**:
   - 为了最大化 `parfor` 效率，将所有单元分成多个块（Chunk）。
   - 每个块包含独立的拓扑 (`T`, `T2E`, `Signs`) 和物理参数 (`Nu`) 数据。
3. **并行内核 (Parallel Kernel)**:
   - 在每个 Worker 中，调用符号积分生成的内核函数 `Ke_curl_curl` 计算单元矩阵 (6x6)。
   - 应用棱边方向符号修正: $K_{local}^{corr} = K_{local} \cdot (\text{sign}_i \cdot \text{sign}_j)$。
   - 将局部索引映射到全局索引，构建三元组 `(rows, cols, vals)`。
4. **全局组装**: 将所有块的三元组拼接，利用 `sparse` 函数一次性生成稀疏矩阵。

**特点**:

- **无竞争并行**: 利用分块策略避免了并行循环中的数据竞争和通信开销。
- **非线性支持**: 接口设计允许外部传入更新后的 $\nu$，无缝对接牛顿-拉夫逊或定点迭代算法。

### 3.2 质量矩阵组装 (`assemble_mass_matrix.m`)

**功能**: 计算 $\mathbf{M}$ 矩阵。支持区域筛选。

**核心优化**:

- **Early Filtering (早期筛选)**: 对于涡流问题，只需在导体区域积分。代码首先检查材料电导率 `Sigma`，仅将 `Sigma > 0` 的单元放入并行队列。这在大量空气单元存在时能显著加速（速度提升可达 10 倍以上）。
- **内核**: 调用 `Me_edge_edge` 计算 $\int \mathbf{N}_i \cdot \mathbf{N}_j$。

### 3.3 非线性雅可比矩阵组装 (`assemble_nonlinear_jacobian.m`)

**功能**: 为牛顿-拉夫逊迭代计算切线刚度矩阵 $\mathbf{K}_{tan}$ 和残差向量 $\mathbf{R}$。

理论:

$$\mathbf{K}_{tan} = \frac{\partial \mathbf{R}}{\partial \mathbf{a}} = \int_\Omega (\nabla \times \mathbf{N})^T \cdot \frac{\partial \mathbf{H}}{\partial \mathbf{B}} \cdot (\nabla \times \mathbf{N}) \, d\Omega$$

其中微分磁阻率张量 $\frac{\partial \mathbf{H}}{\partial \mathbf{B}} = \nu \mathbf{I} + 2 \frac{\partial \nu}{\partial B^2} \mathbf{B}\mathbf{B}^T$。

**实现**:

- **解向量映射**: 使用 `prepare_element_A` 将全局解 $\mathbf{x}$ 快速映射回单元局部解 $\mathbf{a}_e$。
- **B场计算**: 在积分点计算磁通密度 $\mathbf{B}$ 及其模平方 $B^2$。
- **材料微分**: 调用 `eval_material_nu` 同时返回 $\nu$ 和导数 $\frac{d\nu}{dB^2}$。
- **张量构建**: 动态构建各向异性张量，并调用 `Ke_aniso` 内核进行积分。
- **同步输出**: 同时返回刚度矩阵和残差向量，避免重复计算几何量。

### 3.4 绕组耦合向量 (`assemble_winding_vector.m`)

功能: 将复杂的线圈几何映射到有限元网格上，生成源项向量 $\mathbf{W}$。

$$\Psi = \int \mathbf{A} \cdot \mathbf{J}_{coil} \, dV \approx \mathbf{W}^T \mathbf{A}$$

**算法**:

1. **射线追踪**: 对线圈的每条线段，在网格中进行射线追踪或密集采样。
2. **点定位**: 使用 MATLAB `triangulation` 对象的 `pointLocation` 快速查找采样点所在的四面体单元。
3. **形状函数插值**: 在单元内计算矢量基函数 $\mathbf{N}$ 的值。
4. **投影**: 计算 $\mathbf{W}_i = \sum \mathbf{N}_i(\mathbf{x}_k) \cdot d\mathbf{l}_k$。

**特点**:

- **网格无关性**: 线圈几何独立于网格，无需在其路径上强制生成节点。
- **高精度**: 采用自适应步长积分，保证了源项的平滑性。

### 3.5 边界条件施加 (`apply_gauge_strategy.m`, `apply_dirichlet_bc.m`)

**功能**: 处理规范条件（Gauge Condition）和狄利克雷边界条件。

**策略**:

1. **Lagrange Multiplier (拉格朗日乘子法)**:

   - 引入额外的非物理变量（类似节点电压），强制满足 $\nabla \cdot \sigma \mathbf{A} = 0$。

   - 构建鞍点系统：

     $$\begin{bmatrix} \mathbf{K} & \mathbf{C} \\ \mathbf{C}^T & 0 \end{bmatrix} \begin{bmatrix} \mathbf{A} \\ \lambda \end{bmatrix} = \begin{bmatrix} \mathbf{b} \\ 0 \end{bmatrix}$$

   - 自动识别边界节点并固定 $\lambda=0$，保证系统非奇异。

2. **Penalty Method (罚函数法)**:

   - 添加正则化项 $\alpha \mathbf{M}$ 到刚度矩阵，其中 $\alpha$ 是一个小量（如 $10^{-6} \times \text{scale}$）。
   - 不增加系统维度，但在低频下可能引入误差。

**Dirichlet 实现**:

- 采用**置 1 置 0 法**：将边界自由度对应的行/列置零，对角线置 1，RHS 设为定值。保持了矩阵及其对称性，有利于直接求解器。

### 3.6 辅助工具

- `project_source_A_on_edges.m`: 直接计算 Biot-Savart 源场并投影到棱边，用于处理等效磁化电流源项 $\mathbf{J}_{mag} = \nabla \times \mathbf{M}$。
- `assemble_coupling_matrix.m`: 组装混合元耦合矩阵 $\mathbf{C}$，连接棱边元和节点元。
- `assemble_rhs_reduced.m`: 组装由磁化强度差异 $(\nu_0 - \nu)\mathbf{B}_s$ 产生的等效载荷向量。

## 4. 总结 (Summary)

本模块构建了一套完整的、高度优化的有限元组装流水线。其核心优势在于：

1. **模块化设计**: 不同的物理项（刚度、质量、耦合）被拆分为独立函数，便于组合求解各类方程（静场、涡流、瞬态、谐波）。
2. **高性能**: 全面采用向量化运算和并行处理，能高效处理百万级自由度问题。
3. **非线性完备**: 内置了对 Jacobian 和 Residual 的精确计算支持，为非线性求解器的稳健收敛奠定了基础。