% test_v3_hbfem_solve.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Harmonic Balance FEM Solver (HBFEM v3.3 Final)  \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 (Mesh Generation) ---
fprintf('[Step 1] Generating Mesh (5x5x5)...\n');
[X, Y, Z] = meshgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 

% [Critical Fix] 手动更新网格统计信息
% 必须在 generateEdges() 之前设置，因为拓扑构建依赖 NumElements
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);

% (可选) 如果其他模块需要 Faces，可在此处生成，但目前求解器核心不需要
% mesh.Faces = ... 

mesh.generateEdges();

% --- 2. 自由度分发 (DoF Distribution) ---
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

% --- 3. 配置 AFT (谐波设置) ---
fprintf('[Step 3] Configuring AFT Module...\n');
base_freq = 50;
harmonics = [1, 3, 5]; 
aft = AFT(harmonics, [], base_freq); 

fprintf('         Harmonics: %s\n', mat2str(harmonics));
fprintf('         Time Steps: %d\n', aft.NumTimeSteps);

% --- 4. 定义材料 (Stable Rational Model) ---
fprintf('[Step 4] Defining Nonlinear Material...\n');
mu0 = 4*pi*1e-7;
mu_r_init = 2000;
B_sat = 1.5; 

B_arr = linspace(0, 5.0, 100);
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);

MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

% --- 5. 组装器与求解器 ---
assembler = Assembler(mesh, dofHandler);
solver = HBFEMSolver(assembler, aft);

% 使用 MUMPS 加速
solver.LinearSolver.Method = 'Auto';
solver.LinearSolver.MumpsICNTL.i14 = 60; 

% HBFEM 收敛参数
solver.Tolerance = 1e-3;

% --- 6. 负载步进求解 (Load Ramping) ---
fprintf('[Step 5] Running HBFEM Solver with Ramping...\n');

Target_J = 2e5; 
NumSteps = 2; 
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

% --- 7. 后处理 (Post-Processing) ---
fprintf('\n[Step 7] Post-Processing Results...\n');

% 初始化后处理模块
post = PostProcessor(assembler);

% (A) 探测中心点场值
center_point = [0.5, 0.5, 0.5];
[B_center_harm, elem_idx] = post.probeB(X_sol, center_point);
Bz_center = abs(B_center_harm(3, :)); % 取 Z 分量模值

fprintf('\n--- Harmonic Analysis at Center (Element %d) ---\n', elem_idx);
fprintf(' Order |  Bz Amplitude (T)  | Ratio to Fund. (%%) \n');
fprintf('----------------------------------------------\n');

fund_mag = Bz_center(idx_fund);

for k = 1:aft.NumHarmonics
    h = harmonics(k);
    mag = Bz_center(k);
    ratio = (mag / fund_mag) * 100;
    fprintf('   %d   |     %8.4f       |     %6.2f %% \n', h, mag, ratio);
end

% (B) 计算全场 B 并统计最大值
fprintf('\n[Full Field Statistics]\n');
% 这里的 computeElementB 内部 parfor 现在应该能正常工作了，因为 NumElements 已赋值
B_all = post.computeElementB(X_sol); 
B_mag_all = post.computeMagnitude(B_all); % [Ne x K]

max_B_fund = max(B_mag_all(:, idx_fund));
fprintf('  - Max Fundamental B: %.4f T\n', max_B_fund);

% 验证非线性效应
idx_3rd = find(harmonics == 3);
mag_3rd = Bz_center(idx_3rd);

if mag_3rd > 1e-4 * fund_mag
    fprintf('\n[PASS] Significant 3rd harmonic detected (Nonlinearity verified).\n');
else
    fprintf('\n[WARN] 3rd harmonic is too weak.\n');
end

if info.Converged
    fprintf('[PASS] Solver converged successfully.\n');
else
    fprintf('[FAIL] Solver did not reach tolerance.\n');
end