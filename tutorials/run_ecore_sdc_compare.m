% run_ecore_sdc_compare.m
% =========================================================================
% E-Core 变压器瞬态分析对比测试脚本
% 对比对象：
%   1. 标准 BDF1 求解器 (TransientCoupledSolver) - 作为基准 (Reference)
%   2. 高阶 SDC 求解器 (SDCSolver)               - 待验证的高阶算法
%
% 测试目标：
%   验证 SDC 算法在粗时间步长下，能否达到与细步长 BDF1 相当甚至更高的精度，
%   特别是针对变压器涌流这种强非线性刚性问题。
% =========================================================================

clear; clc;
fprintf('=========================================================\n');
fprintf('   Comparison: Standard BDF1 vs High-Order SDC           \n');
fprintf('=========================================================\n');

%% --- 1. 模型共享设置 (Mesh & Physics) ---
% 目的：确保两个求解器使用完全相同的网格、材料和电路参数，保证对比公平性。

fprintf('[Setup] Loading Mesh and Materials...\n');

% 1.1 加载网格
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found: %s', meshFile); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

% 1.2 区域重映射 (Re-tagging)
% 原网格标签分散，需统一为：空气(1), 铁芯(2), 原边线圈(50), 副边线圈(60)
tags_primary = [5, 6, 8, 9]; 
tags_secondary = [3, 4, 7, 10];
tag_core = 2; 
tag_air = 1;
NEW_TAG_PRIM = 50; 
NEW_TAG_SEC = 60;

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

% 1.3 材料定义 (MaterialLib v4.0 PCHIP Smoothness)
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% B-H 曲线数据 (DW465-50 硅钢片)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];

% 创建非线性材料 (自动应用 PCHIP 平滑和真空修正)
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 线性材料 (空气和线圈)
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air) = matLinear; 
matLib(NEW_TAG_PRIM) = matLinear; 
matLib(NEW_TAG_SEC) = matLinear;

% 电导率设置 (本例忽略涡流，设为 0)
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

% 1.4 有限元空间与自由度分发
space_A = FunctionSpace('Nedelec', 1); % 磁矢位 A (边缘元)
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

space_P = FunctionSpace('Lagrange', 1); % 标量电位 P (节点元，用于规范)
dofHandler.distributeDofs(space_P);

% 1.5 线圈几何参数自动提取
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

winding = Winding('Primary', NEW_TAG_PRIM, 3000, 10, area_S, [0,0,0]);
winding.setDirectionField(dir_map);

% 1.6 外电路参数
circuit = struct();
circuit.R = 10; % 电阻 (Ohm)
circuit.L = 0;  % 漏电感 (H)
% 电压源: 500V, 50Hz 正弦波
circuit.V_source_func = @(t) 500 * sin(2 * pi * 50 * t);

% 1.7 边界条件与组装器
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);
assembler = Assembler(mesh, dofHandler);

% 探针位置 (铁芯中心)
probePoint = [0, 0, 0]; 

%% --- 2. 运行 BDF1 (Reference) ---
fprintf('\n[Run 1] Standard BDF1 Solver (Fine Step)...\n');

dt_bdf = 4e-4;           % 细步长 (0.4 ms)
timeSim = 0.01;          % 仿真时长 (0.01 s, 半个周期)
timeSteps_bdf = repmat(dt_bdf, round(timeSim/dt_bdf), 1);

% tic;
% solver_bdf = TransientCoupledSolver(assembler);
% solver_bdf.Tolerance = 1e-3; % BDF1 收敛容差
% [~, info_bdf] = solver_bdf.solve(space_A, space_P, matLib, sigmaMap, ...
%                                  circuit, winding, timeSteps_bdf, ...
%                                  fixedDofs_A, fixedDofs_P, [], [], probePoint);
% time_bdf = toc;
% fprintf('-> BDF1 Completed in %.2f seconds.\n', time_bdf);

%% --- 3. 运行 SDC (High-Order) ---
fprintf('\n[Run 2] High-Order SDC Solver (Coarse Step)...\n');

% SDC 策略: 
% 使用粗步长 (2ms, 即 BDF1 的 5 倍)，通过高阶多项式 (P=4) 来弥补精度。
dt_slab = 2e-3;          
timeSteps_sdc = repmat(dt_slab, round(timeSim/dt_slab), 1);

solver_sdc = SDCSolver(assembler);
solver_sdc.PolyOrder = 4;     % 多项式阶数 P=4 (5个节点)
solver_sdc.MaxSDCIters = 20;  % SDC 最大修正次数
solver_sdc.SDCTolerance = 1e-4; % SDC 收敛容差 (相对/绝对混合判据)

tic;
[~, info_sdc] = solver_sdc.solve(space_A, space_P, matLib, sigmaMap, ...
                                 circuit, winding, timeSteps_sdc, ...
                                 fixedDofs_A, fixedDofs_P, [], [], probePoint);
time_sdc = toc;
fprintf('-> SDC Completed in %.2f seconds.\n', time_sdc);

%% --- 4. 结果对比绘图与误差分析 ---
fprintf('\n[Post] Plotting comparison...\n');

% 提取数据
% t_bdf = cumsum(timeSteps_bdf);
% I_bdf = info_bdf.CurrentHistory;
% B_bdf = info_bdf.ProbeB_History;

t_sdc = cumsum(timeSteps_sdc);
I_sdc = info_sdc.CurrentHistory;
B_sdc = info_sdc.ProbeB_History;

% 绘图
figure('Name', 'BDF1 vs SDC Comparison', 'Position', [200, 200, 1000, 600]);

% 子图 1: 电流对比
subplot(2, 1, 1);
% plot(t_bdf*1e3, I_bdf, 'k-', 'LineWidth', 1.5, 'DisplayName', sprintf('BDF1 (dt=%.1fms)', dt_bdf*1e3));
% hold on;
plot(t_sdc*1e3, I_sdc, 'r-o', 'LineWidth', 1.2, 'MarkerSize', 5, ...
    'DisplayName', sprintf('SDC (dt=%.1fms, P=%d)', dt_slab*1e3, solver_sdc.PolyOrder));
grid on;
% legend('Location', 'best');
xlabel('Time (ms)'); ylabel('Current (A)');
title('Comparision of Inrush Current');

% 子图 2: 磁密对比
subplot(2, 1, 2);
% plot(t_bdf*1e3, B_bdf, 'k-', 'LineWidth', 1.5, 'DisplayName', 'BDF1');
% hold on;
plot(t_sdc*1e3, B_sdc, 'b-d', 'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', 'SDC');
grid on;
% legend('Location', 'best');
xlabel('Time (ms)'); ylabel('|B| at Origin (T)');
title('Comparision of B-Field at Core Center');

% 误差分析 (在共有时间点计算相对误差)
% BDF1 步长 0.4ms, SDC 步长 2.0ms -> 每 5 个 BDF1 步对应 1 个 SDC 步
% step_ratio = round(dt_slab / dt_bdf);
% common_indices = step_ratio : step_ratio : length(I_bdf);

% 截断以匹配较短的向量 (如果仿真没跑完)
% num_common = min(length(common_indices), length(I_sdc));
% if num_common > 0
%     common_indices = common_indices(1:num_common);
%     I_bdf_sub = I_bdf(common_indices);
%     I_sdc_sub = I_sdc(1:num_common);
% 
%     rel_err = norm(I_bdf_sub - I_sdc_sub) / norm(I_bdf_sub);
% 
%     fprintf('\n--- Accuracy Report ---\n');
%     fprintf('Step Ratio (BDF/SDC): 1 : %d\n', step_ratio);
%     fprintf('Relative Error (Current) at Slab Nodes: %.2e%%\n', rel_err * 100);
% else
%     fprintf('\n--- Accuracy Report ---\n');
%     fprintf('Warning: No common time points found for error calculation.\n');
% end

fprintf('\nDone.\n');