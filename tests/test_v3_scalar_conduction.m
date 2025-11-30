% test_v3_scalar_conduction.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Scalar Solver (DC Conduction)                   \n');
fprintf('=========================================================\n');

% --- 1. 网格: 长方体导体 ---
% 尺寸: L=1.0m, Area = 0.1 x 0.1 m^2
L = 1.0; W = 0.1;
[X, Y, Z] = meshgrid(linspace(0, W, 4), linspace(0, W, 4), linspace(0, L, 21));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges(); % Lagrange 其实不需要 Edge，但 Assembler 初始化可能需要

% --- 2. 空间 (Lagrange P1) ---
dofHandler = DofHandler(mesh);
space = FunctionSpace('Lagrange', 1); % 标量节点元
dofHandler.distributeDofs(space);

% --- 3. 物理参数 ---
sigma = 5.8e7; % 铜
SigmaMap = sigma;

% --- 4. 边界条件 (两端电压) ---
% Z=0 面: V=0
% Z=L 面: V=0.1 V
nodes = mesh.P;
tol = 1e-6;
mask_bottom = abs(nodes(3, :) - 0) < tol;
mask_top    = abs(nodes(3, :) - L) < tol;

bcNodes = [find(mask_bottom), find(mask_top)];
bcValues = [zeros(1, sum(mask_bottom)), 0.1 * ones(1, sum(mask_top))];

% --- 5. 求解 ---
assembler = Assembler(mesh, dofHandler);
solver = ScalarSolver(assembler);

V_sol = solver.solve(space, SigmaMap, bcNodes', bcValues');

% --- 6. 结果验证 ---
% 理论电阻 R = L / (sigma * A)
Area = W * W;
R_theo = L / (sigma * Area);
I_theo = 0.1 / R_theo;

fprintf('[Validation] Resistance Calculation\n');
fprintf('  - Geometry: L=%.2f m, Area=%.2e m^2\n', L, Area);
fprintf('  - Theoretical R: %.4e Ohm\n', R_theo);
fprintf('  - Theoretical I: %.4f A\n', I_theo);

% 计算 FEM 总电流
% I = Integral( J . dS ) = Integral( sigma * E . dS )
% 在 Z=L/2 截面计算? 或者利用功率 P = V*I
% P_fem = Integral( sigma * |grad V|^2 ) dV
% I_fem = P_fem / V_diff

lossCalc = LossCalculator(assembler);
% 我们需要给 LossCalculator 增加一个通用的 Grad-Grad 积分接口，或者直接用刚度矩阵
% P = V' * S * V
S = assembler.assembleScalarLaplacian(space, SigmaMap);
P_fem = V_sol' * S * V_sol;

I_fem = P_fem / 0.1;

fprintf('  - FEM Power:     %.4e W\n', P_fem);
fprintf('  - FEM Current:   %.4f A\n', I_fem);

err = abs(I_fem - I_theo) / I_theo;
fprintf('  - Error:         %.2e\n', err);

if err < 1e-3
    fprintf('[PASS] Scalar solver verified.\n');
else
    fprintf('[FAIL] Current mismatch.\n');
end