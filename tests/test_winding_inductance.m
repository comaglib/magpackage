% test_v3_winding_inductance.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Winding Coupling & Inductance Calculation       \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 ---
fprintf('[Step 1] Generating Mesh...\n');
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,2,11));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

% --- 2. 自由度 ---
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs==0
    all = dofHandler.DofMaps(space.toString()); dofHandler.NumGlobalDofs = max(all(:));
end

% --- 3. 定义绕组 ---
% 名称: Coil1, ID: 1, 匝数: 100, R: 1 Ohm, S: 1*1=1m^2, Dir: Z
fprintf('[Step 2] Defining Winding...\n');
coil = Winding('Coil1', 1, 100, 1.0, 1.0, [0, 0, 1]);

% --- 4. 组装 ---
fprintf('[Step 3] Assembling Matrices...\n');
assembler = Assembler(mesh, dofHandler);

% 刚度矩阵 (空气 mu0)
mu0 = 4*pi*1e-7;
K = assembler.assembleStiffness(space, 1/mu0);

% 耦合向量 C
C = assembler.assembleWinding(space, coil);
fprintf('         C vector assembled. nnz = %d\n', nnz(C));

% --- 5. 求解静磁场 ---
fprintf('[Step 4] Solving Magnetostatic Problem...\n');
I_current = 10.0;
F = C * I_current;

% 边界条件 (PEC)
is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);
[K_sys, F_sys] = BoundaryCondition.applyDirichlet(K, F, is_bnd);

% 求解
ls = LinearSolver('Auto');
A_sol = ls.solve(K_sys, F_sys);

% --- 6. 电感计算与验证 ---
fprintf('[Step 5] Calculating Inductance...\n');

% 方法 A: 磁链法 L = Psi / I = (C' * A) / I
Psi = C' * A_sol;
L_flux = Psi / I_current;

% 方法 B: 能量法 L = 2 * W / I^2
W_mag = 0.5 * A_sol' * K * A_sol; % 使用原始 K
L_energy = 2 * W_mag / (I_current^2);

fprintf('  - Current I       : %.2f A\n', I_current);
fprintf('  - Flux Linkage    : %.6e Wb\n', Psi);
fprintf('  - Magnetic Energy : %.6e J\n', W_mag);
fprintf('  -----------------------------------\n');
fprintf('  - Inductance (Flux)  : %.6e H\n', L_flux);
fprintf('  - Inductance (Energy): %.6e H\n', L_energy);

diff_rel = abs(L_flux - L_energy) / L_energy;
fprintf('  - Relative Diff      : %.2e\n', diff_rel);

if diff_rel < 1e-5
    fprintf('[PASS] Inductance calculation methods match.\n');
else
    fprintf('[FAIL] Discrepancy in inductance calculation.\n');
end