% test_v3_nonlinear_solve.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Nonlinear Solver with MUMPS Integration (v4.2)  \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 ---
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,1,6));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.generateEdges();

% --- 2. 自由度分发 ---
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs == 0
    all_dofs = dofHandler.DofMaps(space.toString());
    dofHandler.NumGlobalDofs = max(all_dofs(:));
end

fprintf('[Setup] Global DoFs: %d\n', dofHandler.NumGlobalDofs);

% --- 3. 定义绝对稳定的非线性材料 (Rational Model) ---
mu0 = 4*pi*1e-7;
mu_r_init = 2000;
B_sat = 1.5; 

B_arr = linspace(0, 5.0, 100);
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);

MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

% --- 4. 配置求解器 ---
assembler = Assembler(mesh, dofHandler);
solver = NonlinearSolver(assembler);

% 配置线性求解器 (MUMPS)
solver.LinearSolver.Method = 'Auto'; 
solver.LinearSolver.MumpsICNTL.i14 = 50; 

fprintf('[Setup] Linear Solver Strategy: %s\n', solver.LinearSolver.Method);

% --- 5. 正则化设置 (Regularization) ---
M = assembler.assembleMass(space);
temp_MatMap = 1/mu0;
K_lin = assembler.assembleStiffness(space, temp_MatMap);
ref_val = full(mean(abs(diag(K_lin))));

% [修复] 使用新的通用接口 MatrixK_Add 来施加正则化
epsilon = ref_val * 1e-6;
solver.MatrixK_Add = epsilon * M; 

fprintf('[Setup] Regularization Epsilon: %.2e (Applied via MatrixK_Add)\n', epsilon);

% --- 6. 边界条件 ---
is_bnd_dof = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% --- 7. 负载步进求解 ---
Target_J = 5e5; 
NumSteps = 5; 
x_curr = zeros(dofHandler.NumGlobalDofs, 1);

fprintf('\n[Ramping] Starting Load Ramping (%d steps)...\n', NumSteps);

for step = 1:NumSteps
    scale = step / NumSteps;
    J_step = Target_J * scale;
    
    fprintf('\n--- Step %d/%d (J=%.2e) ---\n', step, NumSteps, J_step);
    
    [x_curr, info] = solver.solve(space, MatLibData, [0,0,J_step], is_bnd_dof, x_curr);
    
    if ~info.Converged
        error('Solver failed at step %d.', step);
    end
end

A_sol = x_curr;

% --- 8. 结果验证 ---
fprintf('\n[Analysis] Checking results...\n');
center = [0.5; 0.5; 0.5];
elem_centers = (mesh.P(:, mesh.T(1,:)) + mesh.P(:, mesh.T(2,:)) + ...
                mesh.P(:, mesh.T(3,:)) + mesh.P(:, mesh.T(4,:))) / 4;
dists = sum((elem_centers - center).^2, 1);
[~, center_elem_idx] = min(dists);

all_dofs_map = dofHandler.DofMaps(space.toString());
dofs = all_dofs_map(:, center_elem_idx);
s = double(mesh.T2E_Sign(:, center_elem_idx));
A_local = A_sol(dofs) .* s;

nodes = mesh.T(:, center_elem_idx);
p_elem = mesh.P(:, nodes);

[J_mat, detJ] = compute_jacobian_tet(p_elem);
q_pt = [0.25; 0.25; 0.25];
[~, curl_ref] = nedelec_tet_p1(q_pt);
C_ref = reshape(curl_ref, 6, 3);
C_phy = (C_ref * J_mat') / detJ;
B_avg = C_phy' * A_local;
B_mag = norm(B_avg);

fprintf('\n[Results]\n');
fprintf('  - Final B_mag at center: %.4f T\n', B_mag);

if B_mag > 0.01 && B_mag < 10.0
    fprintf('  [PASS] Solution is physical and converged.\n');
else
    fprintf('  [FAIL] Solution out of physical range (B=%.2f T).\n', B_mag);
end