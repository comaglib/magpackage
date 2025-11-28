% test_v3_full_system.m
clear; clc;
addpath(genpath('src'));

fprintf('==============================================\n');
fprintf('   MagPackage v3.0 - Full System Test         \n');
fprintf('==============================================\n');

% --- 1. 前处理 (Preprocessing) ---
fprintf('\n[Step 1] Mesh Generation...\n');
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,1,6));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.generateEdges();
mesh.stats();

% --- 2. 自由度分发 (DoF Distribution) ---
fprintf('\n[Step 2] DoF Distribution...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs == 0
    dof_vals = dofHandler.DofMaps(space.toString());
    dofHandler.NumGlobalDofs = max(dof_vals(:));
end

% --- 3. 物理定义 (Physics Definition) ---
fprintf('\n[Step 3] Defining Materials & Sources...\n');
mu0 = 4*pi*1e-7;
nu0 = 1/mu0;
MaterialMap = nu0;
J_mag = 1e6;
SourceMap = [0, 0, J_mag]; 

% --- 4. 矩阵组装 (Assembly) ---
fprintf('\n[Step 4] Assembly (Parallel)...\n');
assembler = Assembler(mesh, dofHandler);

K = assembler.assembleStiffness(space, MaterialMap);
F = assembler.assembleSource(space, SourceMap);

% 组装质量矩阵用于正则化
M = assembler.assembleMass(space);

% --- 5. 边界与正则化 (BCs & Regularization) ---
fprintf('\n[Step 5] Applying BCs and Regularization...\n');
is_bnd_dof = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% [Fix] 必须使用 full() 将稀疏标量转为全标量，否则 mean(sparse) 返回 sparse，导致 fprintf 报错
ref_val = mean(abs(diag(K)));
epsilon = full(ref_val * 1e-6); 

fprintf('[Solver] Adding Mass Regularization (eps = %.2e)...\n', epsilon);
K_reg = K + epsilon * M;

% 施加 Dirichlet BC
[K_sys, F_sys] = BoundaryCondition.applyDirichlet(K_reg, F, is_bnd_dof);

% --- 6. 求解 (Solve) ---
fprintf('\n[Step 6] Solving Linear System...\n');
fprintf('  - System Size: %d x %d\n', size(K_sys,1), size(K_sys,2));

tic;
% warning('off', 'MATLAB:nearlySingularMatrix');
A_sol = K_sys \ F_sys; 
% warning('on', 'MATLAB:nearlySingularMatrix');
solve_time = toc;

fprintf('  - Solved in %.4f seconds.\n', solve_time);
fprintf('  - Solution Norm |A| = %.4e\n', norm(A_sol));

% --- 7. 后处理验证 (Post-Validation) ---
fprintf('\n[Step 7] Physics Validation...\n');

% 验证 1: 磁场能量 W = 0.5 * A' * K_orig * A
W = 0.5 * A_sol' * K * A_sol;
fprintf('  - Magnetic Energy: %.4e J\n', full(W));

% 验证 2: 能量正定性
if W > 1e-12
    fprintf('  [PASS] Energy is positive.\n');
else
    fprintf('  [FAIL] Energy is non-positive or too small! (%.4e)\n', full(W));
end

% 验证 3: 残差检查
free_mask = ~is_bnd_dof;
resid_vec = K_sys * A_sol - F_sys;
residual = norm(resid_vec(free_mask));

fprintf('  - Solver Residual (Free DoFs): %.4e\n', residual);

if residual < 1e-5
    fprintf('  [PASS] Linear system solved accurately.\n');
else
    fprintf('  [WARN] Residual might be too high.\n');
end

fprintf('\n=== Full System Test Complete ===\n');