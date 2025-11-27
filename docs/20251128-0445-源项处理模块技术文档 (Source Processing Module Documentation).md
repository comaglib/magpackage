# 源项处理模块技术文档 (Source Processing Module Documentation)

[toc]

## 1. 概述 (Overview)

本模块负责电磁场有限元分析 (FEM) 中的**源项处理**。其核心任务是将复杂的线圈几何离散化，并基于**毕奥-萨伐尔定律 (Biot-Savart Law)** 计算源电流在任意空间点产生的**磁矢位 (Magnetic Vector Potential,** $\mathbf{A}$**)** 和 **磁感应强度 (Magnetic Flux Density,** $\mathbf{B}$**)**。

这些计算结果通常作为 FEM 求解器的右端项 (RHS)，用于处理开路问题、线圈建模或作为初始场猜测。

## 2. 理论基础 (Theoretical Basis)

### 2.1 毕奥-萨伐尔定律 (Biot-Savart Law)

对于载流线圈，其产生的磁场可以通过对电流路径进行线积分得到。

**磁感应强度** $\mathbf{B}$ **(T):**

$$\mathbf{B}(\mathbf{r}) = \frac{\mu_0 I}{4\pi} \oint_C \frac{d\mathbf{l}' \times (\mathbf{r} - \mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|^3}$$

**磁矢位** $\mathbf{A}$ **(Wb/m):**

$$\mathbf{A}(\mathbf{r}) = \frac{\mu_0 I}{4\pi} \oint_C \frac{d\mathbf{l}'}{|\mathbf{r} - \mathbf{r}'|}$$

其中：

- $\mu_0 = 4\pi \times 10^{-7}$ H/m (真空磁导率)
- $I$: 电流 (A)
- $\mathbf{r}$: 目标场点坐标
- $\mathbf{r}'$: 源点（线圈上的点）坐标
- $d\mathbf{l}'$: 线圈上的微元矢量

### 2.2 线段积分公式 (Analytical Integration for Segments)

为了数值计算，我们将线圈近似为一系列直线段。对于一段从 $\mathbf{P}_1$ 到 $\mathbf{P}_2$ 的直载流导线，在目标点 $\mathbf{r}$ 处的场有解析解。

#### 磁感应强度 $\mathbf{B}$ 的解析解:

$$\mathbf{B} = \frac{\mu_0 I}{4\pi d} (\cos\theta_1 - \cos\theta_2) \cdot \mathbf{e}_{\phi}$$

或者向量形式（代码中采用的形式，更适合计算机）：

$$\mathbf{B} = \frac{\mu_0 I}{4\pi} \frac{(\mathbf{R}_1 \times \mathbf{R}_2)}{|\mathbf{R}_1 \times \mathbf{R}_2|^2} \left( \frac{\mathbf{R}_1 \cdot \mathbf{u}}{R_1} - \frac{\mathbf{R}_2 \cdot \mathbf{u}}{R_2} \right)$$

其中 $\mathbf{R}_1 = \mathbf{r} - \mathbf{P}_1$, $\mathbf{R}_2 = \mathbf{r} - \mathbf{P}_2$, $\mathbf{u} = \mathbf{P}_2 - \mathbf{P}_1$。

#### 磁矢位 $\mathbf{A}$ 的解析解:

$$\mathbf{A} = \frac{\mu_0 I}{4\pi} \ln \left( \frac{R_1 + R_2 + L}{R_1 + R_2 - L} \right) \cdot \hat{\mathbf{u}}$$

其中 $L = |\mathbf{P}_2 - \mathbf{P}_1|$ 是线段长度，$\hat{\mathbf{u}}$ 是线段单位方向向量。此公式精确且避免了数值积分的奇异性。

## 3. 模块功能与实现 (Module Implementation)

### 3.1 线圈几何生成 (`create_racetrack_coil.m`)

**功能**: 生成跑道型 (Racetrack) 线圈的离散线段数据。跑道型线圈由两条直线段和两个半圆弧组成，广泛应用于电机和变压器。

**算法流程**:

1. **参数解析**: 读取线圈中心、法向量、长宽、半径、匝数和电流。
2. **局部坐标建模**:
   - 在局部 XY 平面上生成点云。
   - 直线段 1: `[-L/2, -R]` 到 `[L/2, -R]`
   - 圆弧 1: 半径 R，角度 -90° 到 90°
   - 直线段 2: `[L/2, R]` 到 `[-L/2, R]`
   - 圆弧 2: 半径 R，角度 90° 到 270°
3. **坐标变换 (三维旋转)**:
   - 构建局部坐标系 $(u, v, w)$，其中 $w$ 轴对齐线圈法向量 `Params.Normal`。
   - 使用旋转矩阵 `RotMat` 将局部点转换到全局坐标系。
   - 平移到 `Params.Center`。
4. **数据封装**:
   - `P1`: 线段起点矩阵 (3 x N)
   - `P2`: 线段终点矩阵 (3 x N)
   - `I`: 每段电流 (考虑匝数 N_turns，通常 `I_seg = I_input`)。*注意：某些求解器会在外部处理匝数乘积。*

**特点**:

- **支持任意空间姿态**: 通过法向量定义方向。
- **离散度可控**: `N_seg` 控制圆弧的光滑度。

### 3.2 磁矢位计算 (`compute_biot_savart_A.m`)

**功能**: 并行计算目标点集上的 $\mathbf{A}$ 场。

**实现细节**:

1. **预处理**:
   - 剔除长度极短的无效线段 (`L_sq > 1e-24`)，防止除零错误。
   - 支持标量或向量形式的输入电流 `I`。
   - 计算线段单位方向向量 `UnitDir`。
2. **HPC 优化 (数据分块)**:
   - 为了充分利用多核 CPU 并行计算 (`parfor`)，代码将大量的目标点 (`TargetPoints`) 分割成多个 **块 (Chunks)**（默认 10000 点/块）。
   - 使用 `parallel.pool.Constant` 广播线圈几何数据 (`S1`, `S2`, `I`, `UnitDir`)，避免在每次循环中重复传输大数据，显著降低通信开销。
3. **核心计算 (Kernel)**:
   - 对每个块内的点，遍历所有线段。
   - 应用解析公式 $\mathbf{A} \propto \ln(\dots) \hat{\mathbf{u}}$。
   - **数值稳定性**: 在对数的分母中加入 `1e-12` (`denom(denom < 1e-12) = 1e-12`)，防止当目标点刚好在线段延长线上时出现奇异性 (NaN/Inf)。

**特点**:

- **解析解**: 精度高，无积分误差。
- **并行加速**: 适合大规模场点计算 (如密集网格)。
- **鲁棒性**: 处理了分母为零的边界情况。

### 3.3 磁感应强度计算 (`compute_biot_savart_B.m`)

**功能**: 并行计算目标点集上的 $\mathbf{B}$ 场。

**实现细节**:

1. **预处理与分块**: 与 `compute_biot_savart_A` 相同的逻辑（无效线段剔除、电流兼容、`chunkSize=10000`、`parfor` 并行）。
2. **广播变量**: 使用 `parallel.pool.Constant` 优化内存传输。
3. **核心计算 (Kernel)**:
   - 计算向量 $\mathbf{R}_1, \mathbf{R}_2$。
   - 计算叉积 $\mathbf{C} = \mathbf{R}_1 \times \mathbf{R}_2$ 及其模平方 `cross_sq`。
   - **奇异点保护**: `cross_sq(cross_sq < 1e-20) = 1e-20`。当目标点位于导线上时，叉积为零，会导致 $B \to \infty$。此保护防止程序崩溃，虽然物理上该点场强仍无定义（或需考虑导线半径）。
   - 应用向量形式的毕奥-萨伐尔积分公式。
4. **结果聚合**: 将各块计算结果拼接为完整的 `B_field`。

**特点**:

- **高性能**: 针对 MATLAB 的向量化运算进行了优化。
- **通用性**: 适用于任意形状的电流回路（只要离散为线段）。

### 3.4 串行版 B 场计算 (`compute_biot_savart_B_serial.m`)

**功能**: `compute_biot_savart_B` 的串行版本。

**设计意图**:

- **避免嵌套并行**: 当主求解器（如 `solve_hbfem_voltage`）已经在最外层使用了 `parfor`（例如并行更新材料属性）时，内部调用的函数不能再次开启 `parfor`，否则会导致并行池冲突或效率下降。
- **轻量级**: 去除了分块和 `parallel.pool.Constant` 的开销，适合在已有并行环境的内部被调用，或者处理少量点。

**实现**:

- 逻辑与并行版完全一致，只是将 `parfor` 替换为普通的 `for` 循环，并移除了分块逻辑。

## 4. 总结 (Summary)

| **文件名**                       | **核心功能**      | **关键技术**                       | **适用场景**               |
| -------------------------------- | ----------------- | ---------------------------------- | -------------------------- |
| `create_racetrack_coil.m`        | 几何建模          | 局部坐标系 + 旋转矩阵              | 电机/变压器线圈建模        |
| `compute_biot_savart_A.m`        | 计算 $\mathbf{A}$ | 解析积分 + `parfor` + 数据分块     | 求解器 RHS 组装、边界条件  |
| `compute_biot_savart_B.m`        | 计算 $\mathbf{B}$ | 向量化叉积 + `parfor` + 奇异点保护 | 后处理、场分布查看         |
| `compute_biot_savart_B_serial.m` | 计算 $\mathbf{B}$ | 串行计算                           | **被**并行循环调用的子任务 |

这套源项处理代码构建了一个**从几何到物理场**的完整流水