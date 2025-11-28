% test_v3_hbfem_solve.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Harmonic Balance FEM Solver (HBFEM v3.0)        \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 ---
fprintf('[Step 1] Generating Mesh (5x5x5)...\n');
[X, Y, Z] = meshgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.generateEdges();

% --- 2. 自由度分发 ---
fprintf('[Step 2] Distributing DoFs...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs == 0
    all_dofs = dofHandler.DofMaps(space.toString());
    dofHandler.NumGlobalDofs = max(all_dofs(:));
end
numDofs = dofHandler.NumGlobalDofs;
fprintf('         Base DoFs: %d\n', numDofs);

% --- 3. 配置 AFT ---
fprintf('[Step 3] Configuring AFT Module...\n');
base_freq = 50;
harmonics = [1, 3, 5]; 
aft = AFT(harmonics, [], base_freq); 

fprintf('         Harmonics: %s\n', mat2str(harmonics));
fprintf('         Time Steps: %d\n', aft.NumTimeSteps);

% --- 4. 定义材料 ---
fprintf('[Step 4] Defining Nonlinear Material...\n');
mu0 = 4*pi*1e-7;
mu_r_init = 2000;
B_sat = 1.5; 

B_arr = linspace(0, 5.0, 100);
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);

MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

% --- 5. 求解器配置 ---
assembler = Assembler(mesh, dofHandler);
solver = HBFEMSolver(assembler, aft);

solver.LinearSolver.Method = 'Auto';
solver.LinearSolver.MumpsICNTL.i14 = 60; 
solver.Tolerance = 1e-4;

% --- 6. 负载步进 ---
fprintf('[Step 5] Running HBFEM Solver with Ramping...\n');

Target_J = 2e5; 
NumSteps = 5; 
x_curr = zeros(numDofs, aft.NumHarmonics); 

idx_fund = find(harmonics == 1);
is_bnd_dof = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

for step = 1:NumSteps
    scale = step / NumSteps;
    J_step = Target_J * scale;
    
    fprintf('\n--- Step %d/%d (J_fund=%.2e) ---\n', step, NumSteps, J_step);
    
    SourceMaps = cell(aft.NumHarmonics, 1);
    SourceMaps{idx_fund} = [0, 0, J_step]; 
    
    [x_curr, info] = solver.solve(space, MatLibData, SourceMaps, is_bnd_dof, x_curr);
    
    if ~info.Converged
        error('HBFEM Failed at step %d', step);
    end
end

X_sol = x_curr;

% --- 7. 结果分析 ---
fprintf('\n[Step 7] Analyzing Harmonic Content...\n');

center = [0.5; 0.5; 0.5];
elem_centers = (mesh.P(:, mesh.T(1,:)) + mesh.P(:, mesh.T(2,:)) + ...
                mesh.P(:, mesh.T(3,:)) + mesh.P(:, mesh.T(4,:))) / 4;
dists = sum((elem_centers - center).^2, 1);
[~, center_elem_idx] = min(dists);

all_dofs_map = dofHandler.DofMaps(space.toString());
dofs = all_dofs_map(:, center_elem_idx);
s = double(mesh.T2E_Sign(:, center_elem_idx));

A_elem_harm = X_sol(dofs, :) .* s;

nodes = mesh.T(:, center_elem_idx);
p_elem = mesh.P(:, nodes);
[J_mat, detJ] = compute_jacobian_tet(p_elem);
q_pt = [0.25; 0.25; 0.25];
[~, curl_ref] = nedelec_tet_p1(q_pt);
C_ref = reshape(curl_ref, 6, 3);
C_phy = (C_ref * J_mat') / detJ;

B_harm = C_phy' * A_elem_harm; 
Bz_harm = abs(B_harm(3, :));   

fprintf('\n--- Harmonic Analysis (Center Element) ---\n');
fprintf(' Order |  Bz Amplitude (T)  | Ratio to Fund. (%%) \n');
fprintf('----------------------------------------------\n');

fund_mag = Bz_harm(idx_fund);

for k = 1:aft.NumHarmonics
    h = harmonics(k);
    mag = Bz_harm(k);
    ratio = (mag / fund_mag) * 100;
    
    fprintf('   %d   |     %8.4f       |     %6.2f %% \n', h, mag, ratio);
end

idx_3rd = find(harmonics == 3);
mag_3rd = Bz_harm(idx_3rd);

if mag_3rd > 1e-4 * fund_mag
    fprintf('\n[PASS] Significant 3rd harmonic detected (Nonlinearity verified).\n');
else
    fprintf('\n[WARN] 3rd harmonic is too weak. Is the material saturated?\n');
end

% [Fix] 使用 Converged 标志而不是比较绝对残差
if info.Converged
    fprintf('[PASS] Solver converged successfully.\n');
else
    fprintf('[FAIL] Solver reported non-convergence.\n');
end