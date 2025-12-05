% run_ecore_hbfem.m - E-Core Transformer HBFEM Benchmark
% 
% 描述:
%   该脚本使用谐波平衡有限元法 (HBFEM) 求解 E 型变压器的稳态响应。
%   考虑了铁芯的非线性 B-H 特性和外部电路耦合。
%
% 主要修正记录 [Fix Log]:
%   1. 动态正则化: 自动计算刚度与质量矩阵范数比，设定最优 epsilon 以防止矩阵奇异。
%   2. 求解器配置: 显式设置 MUMPS 为非对称模式，以处理复数耦合矩阵。
%   3. 结果验证: 增加了对计算结果是否爆炸或为零的自动检查。

clear; clc; tic;
fprintf('=========================================================\n');
fprintf('   E-Core HBFEM Analysis (Robust Regularization)         \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
% 加载 COMSOL 格式网格文件，单位转换为米
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges(); % 生成棱边拓扑数据 (用于 Nedelec 单元)

% 定义原始网格中的区域 ID
tags_primary   = [5, 6, 8, 9];   % 原边线圈部分
tags_secondary = [3, 4, 7, 10];  % 副边线圈部分
tag_core       = 2;              % 铁芯
tag_air        = 1;              % 空气

% 定义新的统一 ID，便于后续管理
NEW_TAG_PRIM = 50;
NEW_TAG_SEC  = 60;

% 更新网格中的 RegionTags
mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料与物理属性配置 ---
% 初始化材料库
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 2.1 铁芯材料 (非线性)
% 定义 B-H 曲线数据点 (T, A/m)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 2.2 线性材料 (空气与线圈)
% 这里假设空气和铜线的相对磁导率均为 1
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% 2.3 有限元空间设置
% 使用一阶 Nedelec 棱边元求解磁矢位 A
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

%% --- 3. 线圈几何与方向场设置 ---
% 3.1 自动检测圆形线圈几何参数 (圆心, 半径, 截面积)
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);

% 3.2 计算线圈切向方向场
% 生成每个单元内的电流方向向量，确保电流沿圆环流动
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

% 3.3 创建 Winding 对象
% 参数: 名称, 区域ID, 匝数(300), 电阻(100Ω), 截面积, 全局方向(被场覆盖)
winding = Winding('Primary', NEW_TAG_PRIM, 300, 100, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

%% --- 4. HBFEM 求解器配置 ---
assembler = Assembler(mesh, dofHandler);

% 4.1 定义谐波集合
% 包含基波及主要奇次谐波 [1, 3, 5, 7]
harmonics = [1, 3, 5, 7]; 
baseFreq = 50;            
aft = AFT(harmonics, [], baseFreq); % 初始化时频转换工具

% 4.2 电路参数
R_circuit = 100; % 外电路电阻 Ohm

% 4.3 创建耦合求解器
solver = HBFEMCoupledSolver(assembler, aft, winding, R_circuit);
solver.Tolerance = 1e-3; % 相对收敛容差
solver.MaxIter = 50;     % 最大迭代次数

% [关键配置] 线性求解器设置
% HBFEM 耦合矩阵是非对称且复数的，必须设置 Symmetry = 0
solver.LinearSolver.MumpsSymmetry = 0; 
solver.LinearSolver.MumpsICNTL.i14 = 300; % 增加内存分配预估

% 4.4 定义电压激励
% V(t) = 76 * sin(wt) -> 相量 V = -76j (在基频处)
idx_fund = find(harmonics == 1);
V_harmonics = zeros(length(harmonics), 1);
V_harmonics(idx_fund) = -76j; 

% 4.5 边界条件
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

% 4.6 动态正则化 (Robust Regularization)
% 原因: 本算例中电导率 sigma=0，静态/低频 Curl-Curl 算子存在零空间，需添加微小质量项
fprintf('[Auto-Scale] Estimating Stiffness Norm...\n');

% (1) 估算刚度矩阵 K 的范数 (使用空气磁阻率作为参考)
nu_vec = ones(mesh.NumElements, 1) * (1/(4*pi*1e-7)); 
K_sample = assembler.assembleStiffness(space_A, nu_vec);
K_norm = norm(K_sample, 1);
fprintf('   -> Stiffness Norm ~ %.2e\n', K_norm);

% (2) 估算质量矩阵 M 的范数
M_reg = assembler.assembleMass(space_A);
M_norm = norm(M_reg, 1);
fprintf('   -> Mass Norm ~ %.2e\n', M_norm);

% (3) 计算 epsilon
% 目标: Regularization项 (epsilon*M) 约为 Stiffness项 的 1e-7 倍
% 既保证矩阵非奇异，又不影响物理精度
target_ratio = 1e-7; 
epsilon = (K_norm * target_ratio) / M_norm;
fprintf('   -> Calculated Robust Epsilon = %.4e\n', epsilon);

% 将正则化矩阵传递给求解器
solver.MatrixK_Add = epsilon * M_reg;

%% --- 5. 执行求解 ---
fprintf('[Step 5] Solving HBFEM...\n');
try
    % 调用求解器主函数
    [X_sol, I_sol, info] = solver.solve(space_A, matLib, V_harmonics, fixedDofs_A);
    
    % 结果合理性检查
    I_fund_mag = abs(I_sol(idx_fund));
    if I_fund_mag > 1e10 || isnan(I_fund_mag)
        error('Solution Exploded (I = %.2e). Matrix is still singular.', I_fund_mag);
    elseif I_fund_mag < 1e-6
        warning('Solution is Zero. Check Excitation or Impedance.');
    end
    
catch ME
    fprintf('Error during solve: %s\n', ME.message);
    return;
end

%% --- 6. 后处理与可视化 ---
fprintf('[Step 6] Post-Processing...\n');

% 6.1 重构时域电流波形
% 公式: I(t) = Re( sum( I_k * exp(j * w_k * t) ) )
num_time_pts = 100;
t_cycle = linspace(0, 2/baseFreq, num_time_pts); % 绘制 2 个周期
I_time = zeros(size(t_cycle));

for k = 1:length(harmonics)
    h = harmonics(k);
    omega = 2 * pi * baseFreq * h;
    phasor = I_sol(k);
    I_time = I_time + real(phasor * exp(1j * omega * t_cycle));
end

% 6.2 打印谐波分量数据
fprintf('\n--- Current Harmonics (Amp) ---\n');
for k = 1:length(harmonics)
    fprintf('   Order %d:  %.4f A  (Phase: %.1f deg)\n', ...
        harmonics(k), abs(I_sol(k)), angle(I_sol(k)) * 180/pi);
end

% 6.3 绘图
figure('Name', 'E-Core HBFEM Analysis', 'Position', [100, 100, 1000, 500]);

% 子图1: 稳态电流波形
subplot(1, 2, 1);
plot(t_cycle*1000, I_time, 'r-', 'LineWidth', 2); grid on;
xlabel('Time (ms)'); ylabel('Current (A)');
title('Steady-State Current Waveform');
legend('Total Current (HBFEM)');

% 子图2: 峰值时刻磁通密度 B 分布
% 找到电流峰值时刻
[~, max_idx] = max(abs(I_time));
t_peak = t_cycle(max_idx);

% 重构该时刻的磁矢位 A 场
A_peak = zeros(dofHandler.NumGlobalDofs, 1);
for k = 1:length(harmonics)
    h = harmonics(k);
    phasor = exp(1j * 2*pi*baseFreq*h * t_peak);
    A_peak = A_peak + real(X_sol(:, k) * phasor);
end

% 计算并映射 B 场
post = PostProcessor(assembler);
B_elems = post.computeElementB(A_peak, 'Nedelec_P1');
B_mag = sqrt(sum(abs(B_elems).^2, 1)); % 计算模值
B_nodes = post.mapElementsToNodes(B_mag); % 平滑到节点

subplot(1, 2, 2);
viz = Visualizer(post);
B_core = post.combineRegionData(B_nodes, tag_core); % 仅提取铁芯数据
viz.plotFieldOnSurface(tag_core, B_core(:), 'FaceAlpha', 1.0, 'EdgeColor', 'none');

title(sprintf('|B| at Peak Current (t=%.2f ms)', t_peak*1000));
view(3); axis equal;