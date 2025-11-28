% test_v3_transient_circuit.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Coupled Field-Circuit Solver (RL Step Response) \n');
fprintf('=========================================================\n');

% --- 1. 网格与模型 ---
fprintf('[Step 1] Setup Model...\n');
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,2,11));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs==0
    all = dofHandler.DofMaps(space.toString()); dofHandler.NumGlobalDofs = max(all(:));
end

% 线性材料 (空气)
mu0 = 4*pi*1e-7;
MatLibData(1) = MaterialLib.createLinear(1.0); % mu_r = 1

% --- 2. 预计算电感 L (用于验证) ---
fprintf('[Step 2] Calculating Static Inductance...\n');
assembler = Assembler(mesh, dofHandler);
coil = Winding('Coil1', 1, 100, 1.0, 1.0, [0,0,1]); % N=100, R=1
C = assembler.assembleWinding(space, coil);
K = assembler.assembleStiffness(space, 1/mu0);

% 静态求解 I=1
is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);
[K_stat, F_stat] = BoundaryCondition.applyDirichlet(K, C*1.0, is_bnd);
ls = LinearSolver('Auto');
A_stat = ls.solve(K_stat, F_stat);
L_calc = C' * A_stat;
fprintf('  - Calculated Inductance L: %.6e H\n', L_calc);

% --- 3. 瞬态耦合求解 ---
fprintf('[Step 3] Running Transient Coupled Solver...\n');
R_circuit = 1.0; % Ohm
V_target = 10.0; % Volt

% 求解器配置
solver = TransientCoupledSolver(assembler, coil, R_circuit);

% 时间参数
tau = L_calc / R_circuit;
T_end = 3 * tau; % 模拟 3 个时间常数
dt = tau / 20;
N_steps = ceil(T_end / dt);

TimeParams.dt = dt;
TimeParams.N_steps = N_steps;
TimeParams.V_func = @(t) V_target * (t>=0); % 阶跃电压

[Sol, Info] = solver.solve(space, MatLibData, TimeParams);

% --- 4. 验证结果 ---
fprintf('[Step 4] Validating Results...\n');
t_axis = Sol.Time;
I_fem = Sol.I;

% 解析解: I(t) = V/R * (1 - exp(-R/L * t))
I_analytical = (V_target / R_circuit) * (1 - exp(-t_axis / tau));

% 误差分析
err_norm = norm(I_fem - I_analytical) / norm(I_analytical);
fprintf('  - Max Current FEM: %.4f A\n', max(I_fem));
fprintf('  - Max Current Ana: %.4f A\n', max(I_analytical));
fprintf('  - Relative L2 Error: %.2e\n', err_norm);

% 简单的绘图 (可选)
% plot(t_axis, I_fem, 'o', t_axis, I_analytical, '-'); legend('FEM', 'Analytical');

if err_norm < 1e-2
    fprintf('[PASS] Transient response matches analytical RL circuit.\n');
else
    fprintf('[FAIL] Transient response mismatch.\n');
end