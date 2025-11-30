% test_v3_nonlinear_solve.m
% 验证非线性静磁求解器 (含导电区域 A-V 耦合)
%
% 更新:
%   1. 启用电导率 (Sigma > 0)，模拟导电铁芯。
%   2. 增加 space_V 和 V 场边界条件。
%   3. 验证 A-V 耦合下的非线性收敛性。

clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Nonlinear Magnetostatic Solve (Conductive)      \n');
fprintf('=========================================================\n');

% --- 1. 生成网格 (铁芯线圈) ---
L_iron = 0.05;
[X, Y, Z] = meshgrid(linspace(-L_iron, L_iron, 7), ...
                     linspace(-L_iron, L_iron, 7), ...
                     linspace(0, 0.2, 5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); % 1=Iron
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();
fprintf('[Mesh] Elements: %d, Edges: %d\n', mesh.NumElements, size(mesh.Edges,2));

% --- 2. 材料定义 (非线性) ---
matData = struct();
matData.Type = 'Nonlinear';
matData.nu0 = 1/(1000 * 4*pi*1e-7); 
matData.MaxBSq = 2.5^2;

% B-H 曲线
H_samples = [0, 100, 500, 1000, 5000, 10000, 50000];
B_samples = [0, 0.5, 1.0, 1.3,  1.6,  1.7,   1.8];
nu_samples = H_samples ./ (B_samples + 1e-9); 
nu_samples(1) = nu_samples(2); 

pp = spline(B_samples.^2, nu_samples);
matData.SplineNu = pp;
matData.SplineDNu = fnder(pp, 1); % 必须提供导数

MatLibData = containers.Map({1}, {matData});

% [Change] 设置非零电导率，激活 A-V 形式
SigmaMap = containers.Map({1}, {5.8e7}); 

% --- 3. 求解器配置 ---
dofHandler = DofHandler(mesh);

% (A) 磁矢位
space_A = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space_A);

% (V) 电标位 [新增]
space_V = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_V, 1); % 仅在导电区分配

% [Check] 确认自由度
num_A = dofHandler.SpaceLocalSizes(space_A.toString());
num_V = dofHandler.SpaceLocalSizes(space_V.toString());
fprintf('[DoF] A-DoFs: %d, V-DoFs: %d\n', num_A, num_V);

if num_A == 0, error('Active DoFs is 0.'); end

% A 边界: 外边界 A=0
mask_bnd_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_A = find(mask_bnd_A);

% V 边界: 底部接地 (V=0) [新增]
% 必须固定 V 以确定电势参考点，否则矩阵奇异
nodes = mesh.P;
mask_bottom = abs(nodes(3, :) - 0) < 1e-6;
% 查找这些节点对应的 V 自由度
map_V = dofHandler.EntityMaps(space_V.toString());
fixedDofs_V = [];
if ~isempty(map_V)
    idx_nodes = find(mask_bottom);
    valid_dofs = map_V(idx_nodes);
    fixedDofs_V = valid_dofs(valid_dofs > 0);
end
fprintf('[BC] Fixed A: %d, Fixed V: %d\n', length(fixedDofs_A), length(fixedDofs_V));

assembler = Assembler(mesh, dofHandler);
solver = NonlinearSolver(assembler);

% 求解器设置
solver.AutoRegularization = true;
solver.RegScale = 1e-6; 
solver.Relaxation = 0.8; 

% --- 4. 载荷爬升测试 (Ramping) ---
fprintf('\n[Ramping] Starting Load Ramping (3 steps)...\n');

J_target = 1e7; 
steps = 3;
x_curr = []; 

for i = 1:steps
    factor = i / steps;
    fprintf('  Step %d: Load factor %.2f\n', i, factor);
    
    current_J = [factor * J_target, 0, 0];
    SourceMap = containers.Map({1}, {current_J});
    
    % [Call] 传入 space_V 和 fixedDofs_V
    [A_sol, V_sol, info] = solver.solve(space_A, space_V, MatLibData, SigmaMap, SourceMap, fixedDofs_A, fixedDofs_V, x_curr);
    
    if ~info.Converged
        warning('Solver failed to converge at step %d', i);
        break;
    end
    
    x_curr = info.Solution;
    fprintf('    -> Converged in %d iters. Res: %.2e\n', info.Iterations, info.FinalResidual);
end

% --- 5. 验证 ---
max_A = max(abs(A_sol));
max_V = max(abs(V_sol));
fprintf('\n[Result] Max |A|: %.4e, Max |V|: %.4e\n', max_A, max_V);

if max_A > 1e-9
    fprintf('[PASS] Nonlinear solver produced valid solution.\n');
else
    fprintf('[FAIL] Solution is trivial (Zeros).\n');
end