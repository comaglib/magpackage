% run_ecore_transient.m - E-Core Transformer Transient Benchmark
% 
% 描述:
%   模拟 E 型变压器在正弦电压源激励下的瞬态响应（场路耦合）

clear; clc; tic;

fprintf('=========================================================\n');
fprintf('   E-Core Transient Analysis (No Eddy, No V-Dofs)        \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
meshFile = 'data/meshData/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end

mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

% --- 区域 ID 合并 ---
tags_primary   = [5, 6, 8, 9];
tags_secondary = [3, 4, 7, 10];
tag_core       = 2;
tag_air        = 1;

NEW_TAG_PRIM = 50;
NEW_TAG_SEC  = 60;

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料与物理属性配置 ---

% 2.1 铁芯非线性 B-H
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];

matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 2.2 线性材料 (空气/线圈)
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% 2.3 电导率 (全局设为 0，不分析涡流)
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags)
    sigmaMap(all_tags(i)) = 0.0;
end

% 2.4 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

space_P = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_P);

%% --- 3. 线圈与电路耦合设置 ---
fprintf('[Step 3] Setup Field-Circuit Coupling...\n');

% 3.1 自动计算线圈几何信息
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);

% 3.2 计算切向电流方向场 (顺时针)
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

fprintf('   -> Coil Config: R=%.4fm, Clockwise Direction Set.\n', radius);

% 3.3 创建 Winding 对象
winding = Winding('Primary', NEW_TAG_PRIM, 300, 100, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

% 3.4 外电路定义
circuit = struct();
circuit.R = 100; 
circuit.L = 0;
circuit.V_source_func = @(t) 76 * sin(2 * pi * 50 * t);

%% --- 4. 求解器配置 ---
assembler = Assembler(mesh, dofHandler);
solver = TransientCoupledSolver(assembler);
solver.Tolerance = 1e-2;
solver.RelTolerance = 1e-2;
solver.LinearSolver.MumpsSymmetry = 0;
solver.LinearSolver.MumpsICNTL.i14 = 300;

dt = 5e-4;
timeSteps = repmat(dt, round(0.04/dt), 1);

fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);

plotFunc = @(t, I, t_vec, I_vec) plot(t_vec, I_vec, 'r-o', 'LineWidth', 1.5); grid on;

%% --- 5. 求解 ---
fprintf('[Step 5] Solving...\n');

[sols, info] = solver.solve(space_A, space_P, matLib, sigmaMap, ...
                            circuit, winding, timeSteps, ...
                            fixedDofs_A, fixedDofs_P, [], plotFunc);

%% --- 6. 后处理 ---
fprintf('[Step 6] Post-Processing...\n');
post = PostProcessor(assembler);
viz  = Visualizer(post);

% 提取电流
num_steps = length(sols);
time_vec = cumsum(timeSteps);
current_vec = info.CurrentHistory;

% 绘制电流曲线
figure('Name', 'E-Core Transient Analysis', 'Position', [100, 100, 1000, 500]);
subplot(1, 2, 1);
plot([0;time_vec], [0;current_vec], 'b-o', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('Current (A)');
title('Primary Coil Current (No Eddy)');

% 绘制铁芯 B 场
target_t = 0.006;
[~, step_idx] = min(abs(time_vec - target_t));
A_sol = sols(1:dofHandler.NumGlobalDofs);

subplot(1, 2, 2);
B_elems = post.computeElementB(A_sol, 'Nedelec_P1');
B_mag = sqrt(sum(abs(B_elems).^2, 1));
B_nodes = post.mapElementsToNodes(B_mag);
B_core = post.combineRegionData(B_nodes, tag_core);

viz.plotFieldOnSurface(tag_core, B_core(:), 'FaceAlpha', 1.0, 'EdgeColor', 'none');
title(sprintf('|B| at t=%.3fs', time_vec(step_idx)));
view(3); axis equal;

toc;