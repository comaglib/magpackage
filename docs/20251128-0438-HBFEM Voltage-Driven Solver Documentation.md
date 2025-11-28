# **HBFEM Voltage-Driven Solver Documentation**

[toc]

## **1\. 概述 (Overview)**

solve_hbfem_voltage 是一个基于 **谐波平衡有限元法 (Harmonic Balance FEM, HBFEM)** 的高性能电磁场求解器。它专门用于解决包含 **非线性铁磁材料** 和 **外部电路耦合** 的低频电磁问题（如变压器、电抗器、电机等）。

### **核心能力**

* **强非线性处理**: 能够处理深度饱和的 B-H 曲线，准确捕捉磁滞回线（通过复数磁阻率扩展）和饱和效应。  
* **场路耦合 (Field-Circuit Coupling)**: 直接求解电压方程 $V = RI + \frac{d\Psi}{dt}$，支持电压源驱动，自动计算非线性励磁电流。  
* **谐波分析**: 直接在频域求解稳态响应，一次性获得基波及高次谐波分量，避免了时域法漫长的暂态过程。

## **2\. 理论基础 (Theoretical Basis)**

### **2.1 控制方程 (Governing Equations)**

从麦克斯韦方程组出发，引入磁矢位 $\boldsymbol{A}$ (Magnetic Vector Potential)，忽略位移电流（低频近似）：

$$\nabla \times \left( \nu(\boldsymbol{B}) \nabla \times \boldsymbol{A} \right) = \boldsymbol{J}_s - \sigma \frac{\partial \boldsymbol{A}}{\partial t} - \sigma \nabla V$$  
其中：

* $\nu(\boldsymbol{B})$: 非线性磁阻率 (Reluctivity, $1/\mu$)。  
* $\boldsymbol{J}_s$: 源电流密度。  
* $\sigma$: 电导率。

对于电压驱动的线圈，源电流密度 $\boldsymbol{J}_s$ 不是已知的，而是由电路方程决定：

$$U(t) = R \cdot i(t) + L_{leak} \frac{di}{dt} + \frac{d\Psi}{dt}$$  
其中 $\Psi$ 是线圈磁链，由场变量 $\boldsymbol{A}$ 决定：$\Psi = \oint_{coil} \boldsymbol{A} \cdot d\boldsymbol{l}$。

### **2.2 谐波平衡法 (Harmonic Balance Method)**

由于 $\nu(\boldsymbol{B})$ 是随时间变化的非线性参数，传统的频域法无法直接应用。HBFEM 将所有变量展开为傅里叶级数：

$$\boldsymbol{A}(t) = \sum_{k} \boldsymbol{A}_k e^{j \omega_k t}, \quad \nu(t) = \sum_{m} \nu_m e^{j \omega_m t}$$  
然而，直接求解所有谐波的耦合方程组会导致矩阵极其庞大（Dense Matrix）。本求解器采用 **分解式 (Decomposed) / 定点迭代 (Fixed-Point)** 策略：

将非线性磁阻率分解为 直流分量 $\nu_{dc}$ 和 时变扰动分量 $\nu_{pert}(t)$：

$$\nu(t) = \nu_{dc} + (\nu(t) - \nu_{dc})$$  
代入控制方程并移项，得到第 $k$ 次谐波的线性化方程：

$$\nabla \times (\nu_{dc} \nabla \times \boldsymbol{A}_k) + j\omega_k \sigma \boldsymbol{A}_k = \boldsymbol{J}_k - \nabla \times \boldsymbol{M}_{resid, k}$$  
其中 $\boldsymbol{M}_{resid}$ (Residual Magnetization) 是非线性残差项，视为等效磁化电流源：

$$\boldsymbol{M}_{resid}(t) = (\nu(t) - \nu_{dc}) \cdot \boldsymbol{B}(t)$$

### **2.3 物理修正源项 (Corrected Magnetization Source)**

为了精确处理开路磁芯等强退磁问题，我们引入真空磁阻率 $\nu_0$ 作为参考，将源项重写为：

$$\nabla \times (\nu_{dc} \boldsymbol{B}_{total}) = \nabla \times (\nu_0 \boldsymbol{B}_s)$$  
这等价于在右端项引入一个修正的磁化源：

$$\boldsymbol{M}_{linear} = (\nu_0 - \nu_{dc}) \boldsymbol{B}_s$$  
这保证了当材料退化为空气时，方程精确回归到 Biot-Savart 定律。

## **3\. 算法流程 (Algorithm Workflow)**

求解器采用 **连续负载推进 (Continuous Source Stepping)** 结合 **对偶右端项 (Dual-RHS)** 的迭代策略。

### **Step 1: 初始化 (Initialization)**

1. **频率设置**: 确定基频和谐波次数 (如 1, 3, 5, 7, 9)。  
2. **FFT 参数**: 计算满足采样定理的时间步数 $N_t$。  
3. **几何预计算**:  
   * 计算线圈的 **Winding Vector** ($\boldsymbol{W}$)，用于计算磁链 $\Psi = \boldsymbol{W}^T \boldsymbol{A}$。  
   * 计算单位电流产生的 **源磁场** $\boldsymbol{B}_s$ (Biot-Savart)。  
4. **初值猜测**:  
   * $\nu_{dc}$ 初始化为材料的线性磁导率对应的磁阻率。  
   * 电流 $I$ 根据线性阻抗估算。

### **Step 2: 连续负载推进 (Continuous Ramping Loop)**

为了解决强非线性带来的收敛困难，电压负载从 0% 线性增加到 100%。

**For** iter = 1 : MaxIter

1. **更新负载**: $V_{target} = V_{rated} \times \min(\frac{iter}{RampSteps}, 1.0)$。  
2. **组装 LHS**: 使用当前的 $\nu_{dc}$ 组装刚度矩阵 $\mathbf{K}$。  
3. **材料属性更新 (Parallel)**:  
   * 时域重构: $\boldsymbol{A}(t) = \text{IFFT}(\boldsymbol{A}_k)$。  
   * 计算 B 场: $\boldsymbol{B}(t) = \nabla \times \boldsymbol{A}(t) + \boldsymbol{B}_s \cdot i(t)$。  
   * 查表: 根据 $B-H$ 曲线更新 $\nu(t)$。  
   * 频域变换: 计算 $\nu_{dc}$ 和非线性残差 $\boldsymbol{M}_{resid}$ 的谐波分量。  
4. **组装 RHS**:  
   * 线性源项: $\boldsymbol{F}_{lin} = \int \nabla \times \boldsymbol{N} \cdot (\nu_0 - \nu_{dc})\boldsymbol{B}_s dV$。  
   * 非线性残差: $\boldsymbol{F}_{res} = -\int \nabla \times \boldsymbol{N} \cdot \boldsymbol{M}_{resid} dV$。  
5. **求解线性方程 (Dual-RHS)**:  
   * 分别求解残差响应 $\boldsymbol{A}_{res}$ 和单位电流响应 $\boldsymbol{A}_{unit}$。  
   * 利用叠加原理建立电路方程：  
     $$V = Z_{leak} I + j\omega (\Psi_{res} + I \cdot \Psi_{unit})$$  
   * 解出电流 $I$，合成最终磁矢位 $\boldsymbol{A} = \boldsymbol{A}_{res} + I \cdot \boldsymbol{A}_{unit}$。  
6. **收敛检查**: 计算相对误差。如果负载为 100% 且误差小于容差，则退出。  
7. **松弛更新**: $X_{new} = (1-\alpha)X_{old} + \alpha X_{new}$。

## **4\. 关键技术实现 (Implementation Details)**

### **4.1 稀疏矩阵掩码 (Sparse Masking for BC)**

apply_dirichlet_bc_multi 函数经过深度优化。传统的 K(:, dofs) = 0 操作在 MATLAB 中非常缓慢。本求解器构建对角掩码矩阵 $\mathbf{D}$ (边界处为0，内部为1)：

$K_new = D * K * D + I_fixed;$

这种纯矩阵乘法操作利用了稀疏矩阵的存储特性，速度提升显著。

### **4.2 电流阻尼与限幅 (Current Damping & Limiter)**

在电压驱动模式下，电感 $L$ 随饱和急剧下降，导致电流 $I = V / (j\omega L)$ 可能出现数值爆炸。

* **Soft-Start Limiter**: 强制限制单步迭代中电流幅值的增长倍数（如最大 2.5 倍），防止物理量逃逸。  
* **Low-Pass Filter**: $I_{new} = 0.7 I_{calc} + 0.3 I_{old}$，引入惯性，抑制高频震荡。

### **4.3 物理一致性修正**

代码中严格区分了 **源场 (Source Field)** 和 **反应场 (Reaction Field)**。

* **源场** $\boldsymbol{B}_s$: 由线圈电流在真空中产生，通过 Biot-Savart 积分预计算。  
* **反应场** $\boldsymbol{B}_r$: 由铁磁材料磁化产生，通过 FEM 求解。  
* 总场: $\boldsymbol{B}_{total} = \boldsymbol{B}_s + \boldsymbol{B}_r$。  
  这种分离技术避免了在源区域进行奇异性积分，大幅提高了精度。

## **5\. 使用指南 (Usage Guide)**

### **输入参数结构体**

``` matlab
Model.Mesh          % 网格对象 (P, T, Edges)  
Model.Materials     % 材料库 (含 B-H 曲线 spline)  
Coil.P1, Coil.P2    % 线圈几何线段  
Coil.Turns          % 匝数  
Circuit.Voltage     % 电压向量 [V1, V3, V5...]  
HBFEMParams.Harmonics % 需要求解的谐波阶次 [1, 3, 5...]
```

### **典型调用示例**

``` matlab
% 定义 50Hz, 300V 电压源，求解 1,3,5 次谐波  
Circuit.Voltage = [300; 0; 0];  
Params.Frequency = 50;  
Params.Harmonics = [1, 3, 5];

[Sol, Info] = solve_hbfem_voltage(Model, Coil, Circuit, Params);

% 后处理  
I_fundamental = abs(Sol.Current(1));  
I_3rd_harmonic = abs(Sol.Current(2));
```

## **6\. 故障排查 (Troubleshooting)**

| 现象 | 可能原因 | 解决方案 |
| :---- | :---- | :---- |
| **迭代不收敛 (Oscillation)** | 电压过高导致深度饱和 | 减小 Circuit.R (阻尼作用) 或增加 HBFEMParams.TimeSteps。 |
| **电流非常小 (uA 级)** | 磁路未闭合 (Open Circuit) | 检查网格连接性，确保铁芯与空气节点重合。 |
| **无三次谐波** | 材料工作在通过线性区 | 增加电压 Circuit.Voltage 或检查 B-H 曲线膝点。 |
| **计算速度慢** | 网格过密或谐波过多 | 减少网格数量；减少谐波阶次。 |

