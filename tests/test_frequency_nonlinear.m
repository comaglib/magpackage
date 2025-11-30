% test_v3_frequency_nonlinear.m
% 验证非线性频域求解器 (Coil + Nonlinear Core with Eddy Currents)
%
% 模型说明:
%   - Region 1 (Core): 铁芯，非线性，导电 (sigma > 0)，无源。
%   - Region 2 (Coil): 空气/线圈区域，线性，不导电 (sigma = 0)，有电流源。
%   - 物理现象: 线圈电流产生磁场 -> 磁场进入铁芯 -> 铁芯产生感应涡流 (Eddy Current)。

clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Nonlinear Eddy Current (Coil + Core)            \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 (双区域) ---
% Core: [-0.05, 0.05] x [-0.05, 0.05] x [-0.2, 0.2]
% Coil: [-0.10, 0.10] x [-0.10, 0.10] x [-0.2, 0.2] (包裹 Core)

fprintf('[Mesh] Generating Composite Mesh...\n');
% 内部铁芯点
[X1, Y1, Z1] = meshgrid(linspace(-0.05, 0.05, 5), linspace(-0.05, 0.05, 5), linspace(-0.2, 0.2, 9));
DT1 = delaunayTriangulation(X1(:), Y1(:), Z1(:));

% 外部线圈点 (挖空? 简单起见，我们在外部堆叠两个块，或者直接定义不重叠的区域)
% 为简化，我们建立两个并列的块 (Side-by-side)，磁通通过空气耦合
% Block 1: Core (0,0,0)
% Block 2: Source (0.15, 0, 0) -> 这不产生包围感应。
% 
% 更好方案: 同心结构。
% 我们手动构建两套不重叠的四面体比较麻烦。
% 简化模型: 上下结构。
% Region 1 (Core): z in [0, 0.1]
% Region 2 (Coil): z in [0.1, 0.2] (Air gap source)
% 这样磁通会穿过界面进入 Core。

[X, Y, Z] = meshgrid(linspace(-0.05, 0.05, 7), linspace(-0.05, 0.05, 7), linspace(0, 0.2, 11));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';

% 定义区域 Tag
% z < 0.1 -> Tag 1 (Core)
% z >= 0.1 -> Tag 2 (Coil Source)
centers = (mesh.P(:, mesh.T(1,:)) + mesh.P(:, mesh.T(2,:)) + mesh.P(:, mesh.T(3,:)) + mesh.P(:, mesh.T(4,:))) / 4;
mask_core = centers(3, :) < 0.1;
mesh.RegionTags = ones(1, size(mesh.T, 2)) * 2; % Default 2
mesh.RegionTags(mask_core) = 1;

mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

fprintf('[Mesh] Elements: %d (Core: %d, Coil: %d)\n', mesh.NumElements, sum(mesh.RegionTags==1), sum(mesh.RegionTags==2));

% --- 2. 材料定义 ---
mu0 = 4*pi*1e-7;

% Tag 1: Core (Nonlinear, Conductive)
matCore.Type = 'Nonlinear';
matCore.nu0 = 1/(2000 * mu0); 
% B-H 曲线
B_knots = [0, 0.5, 1.0, 1.4, 1.6, 2.0, 5.0]; 
H_vals  = [0, 200, 500, 2000, 10000, 50000, 500000]; 
nu_vals = H_vals ./ (B_knots + 1e-9); nu_vals(1) = nu_vals(2); 
matCore.SplineNu = spline(B_knots.^2, nu_vals);
matCore.MaxBSq = max(B_knots)^2;

% Tag 2: Coil (Linear, Non-conductive)
matCoil.Type = 'Linear';
matCoil.Nu_Linear = 1/mu0;

MatLib = containers.Map({1, 2}, {matCore, matCoil});

% --- 3. 物理场配置 ---
% Tag 1: Sigma = 1e6 (Eddy currents allowed)
% Tag 2: Sigma = 0 (No eddy currents)
sigma_core = 5e3; 
SigmaMap = containers.Map({1, 2}, {sigma_core, 0});

% Tag 2: Source Current (Driving Field)
% 在上方区域施加 X 方向电流，产生 Y 方向磁场
J_drive = 1e6; 
freq = 50; 
SourceMap = containers.Map({1, 2}, {[0,0,0], [J_drive, 0, 0]});

% --- 4. 空间与自由度 ---
dofHandler = DofHandler(mesh);

% (A) 磁矢位
space_A = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space_A);
% A 边界: 外表面 PEC (A x n = 0 -> A_tan = 0 for full Dirichlet? Or naturally handled?)
% 为简化，固定外边界 A=0 (限制磁通区域)
mask_bnd_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_A = find(mask_bnd_A);

% (V) 电标位
% 仅在 Tag 1 (Core) 分配 V 自由度
space_V = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_V, 1); 

% V 边界: 仅仅为了定标 (Gauge fixing)。
% 如果整个导体的边界都未连接电势源，V 是悬浮的 (Floating)。
% 必须固定至少一个点或边界，防止矩阵奇异。
% 这里我们将 Core 的底部 (z=0) 接地 V=0。
nodes = mesh.P;
mask_nodes_V = (abs(nodes(3, :) - 0.0) < 1e-6); % z=0 surface
% 获取这些节点对应的 V-DoF
% 注意: findOuterBoundaryDofs 返回的是所有外表面。我们这里需要更精细的选择。
% DofHandler 的 EntityMaps 存储了 NodeID -> GlobalDoF 的映射。
map_V = dofHandler.EntityMaps(space_V.toString());
fixedDofs_V = [];
if ~isempty(map_V)
    valid_nodes = find(mask_nodes_V);
    % 过滤掉未分配 DoF 的节点 (如 Tag 2 的节点)
    valid_dofs = map_V(valid_nodes); 
    fixedDofs_V = valid_dofs(valid_dofs > 0);
end

fprintf('[DoF Info] A: %d, V: %d (Fixed V: %d)\n', ...
    dofHandler.SpaceLocalSizes(space_A.toString()), ...
    dofHandler.SpaceLocalSizes(space_V.toString()), length(fixedDofs_V));

% --- 5. 求解 ---
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, freq);

% 策略配置
solver.Tolerance = 1e-4;
solver.Relaxation = 0.5;

fprintf('\n[Solver] Running Coupled Solve...\n');
try
    [A_sol, V_sol] = solver.solve(space_A, space_V, MatLib, SigmaMap, SourceMap, fixedDofs_A, fixedDofs_V);
catch ME
    fprintf('[FAIL] Solver Crashed: %s\n', ME.message);
    return;
end

% --- 6. 验证 ---
if isempty(A_sol), fprintf('[FAIL] Solution is empty.\n'); return; end

max_A = full(max(abs(A_sol)));
max_V = full(max(abs(V_sol)));
fprintf('\n[Result] Max |A|: %.4e, Max |V|: %.4e\n', max_A, max_V);

% 判定: 
% 1. A 应该非零 (有源)
% 2. V 应该非零 (有感应涡流)
%    由于是对称模型且 J_src 是线性的，感应电流 J_eddy = -sigma(jwA + gradV)。
%    如果几何不对称或边界条件导致电荷积累，V 将显著非零。
if max_A > 1e-9 && max_V > 1e-9
    fprintf('[PASS] Significant eddy currents detected (|V| > 0).\n');
elseif max_A > 1e-9
    fprintf('[WARN] Magnetic field exists, but V is trivial (Symmetric eddy loops?).\n');
    fprintf('       This is physically possible if J_eddy is purely divergence-free without charge accumulation.\n');
else
    fprintf('[FAIL] Solution is trivial.\n');
end