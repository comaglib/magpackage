% test_v3_transient_eddy.m
% 测试非线性瞬态涡流求解器 (Class API)
% ---------------------------------------------------------
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Transient Nonlinear Eddy Current (v3.1)         \n');
fprintf('=========================================================\n');

% --- 1. 网格与模型 ---
fprintf('[Step 1] Generating Mesh (6x6x6)...\n');
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,1,6));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

% --- 2. 自由度 ---
fprintf('[Step 2] Distributing DoFs...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
fprintf('         Total DoFs: %d\n', dofHandler.NumGlobalDofs);

% --- 3. 材料 ---
fprintf('[Step 3] Defining Nonlinear Material...\n');
mu0 = 4*pi*1e-7;
B_curve = linspace(0, 3.0, 50);
mu_r_curve = 1000 ./ (1 + 5*B_curve.^2) + 1; % 简单的饱和曲线
H_curve = B_curve ./ (mu_r_curve * mu0);
MatLibData(1) = MaterialLib.createNonlinear(B_curve, H_curve);

% --- 4. 求解参数 ---
fprintf('[Step 4] Configuring Transient Solver...\n');
freq = 50;
sigma = 1e4; % S/m
dt = 1 / (20 * freq); % 每个周期 20 步
N_steps = 20;

TimeParams.dt = dt;
TimeParams.N_steps = N_steps;
TimeParams.Sigma = sigma;

% 源项: 正弦电流 Jz
J_mag = 1e6;
TimeParams.SourceFunc = @(t) [0, 0, J_mag * sin(2*pi*freq*t)];

fprintf('         Freq: %d Hz, Conductivity: %.2e S/m\n', freq, sigma);
fprintf('         Steps: %d, dt: %.2e s\n', N_steps, dt);

% --- 5. 执行求解 (使用新类接口) ---
fprintf('[Step 5] Running Solver...\n');
assembler = Assembler(mesh, dofHandler);

% [NEW] 实例化求解器对象
solver = TransientNonlinearSolver(assembler);

% [NEW] 调用 solve 方法
[Solution, Info] = solver.solve(space, MatLibData, TimeParams);

% --- 6. 结果检查 ---
fprintf('[Step 6] Analyzing Results...\n');
A_final = Solution.A_history{end};
max_A = max(abs(A_final));
fprintf('   Max A (End): %.4e\n', max_A);

% 简单的物理检查: 能量是否有限且未发散
if max_A > 1e-9 && max_A < 1e5
    fprintf('[PASS] Transient solution is within physical range.\n');
else
    fprintf('[FAIL] Solution diverged or is zero.\n');
end