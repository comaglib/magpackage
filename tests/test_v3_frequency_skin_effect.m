% test_v3_frequency_skin_effect.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Frequency Domain Solver (Skin Effect)           \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 ---
fprintf('[Step 1] Generating Mesh...\n');
L = 0.1; % 导体宽度
[X, Y, Z] = meshgrid(linspace(-L/2, L/2, 11), ...
                     linspace(-L/2, L/2, 11), ...
                     linspace(0, 0.5, 5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));

mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 

% [Critical Fix] 显式更新网格统计信息
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);

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
fprintf('         Total DoFs: %d\n', dofHandler.NumGlobalDofs);

% --- 3. 物理参数设置 ---
f = 1000;             
omega = 2 * pi * f;
mu0 = 4*pi*1e-7;
mu_r = 1000;          
sigma = 1e7;          

% 理论趋肤深度
delta = sqrt(2 / (omega * (mu_r * mu0) * sigma));
fprintf('[Physics] Frequency: %.1f Hz\n', f);
fprintf('          Skin Depth: %.4f mm (Mesh width: %.1f mm)\n', delta*1000, L*1000);

NuMap = 1/(mu_r * mu0); 
SigmaMap = sigma;       

J_drive = 1e6;
SourceMap = [0, 0, J_drive]; 

is_bnd_dof = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% --- 4. 求解 ---
fprintf('[Step 3] Solving Frequency Domain System...\n');
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);

% 使用 MUMPS
solver.LinearSolver.Method = 'Auto';

A_phasor = solver.solve(space, NuMap, SigmaMap, SourceMap, is_bnd_dof);

% --- 5. 后处理 ---
fprintf('[Step 4] Analyzing Results (Skin Effect)...\n');
post = PostProcessor(assembler);

% 探测中心点
pt_center = [0, 0, 0.25];
% 注意: probeB 返回的是 B=Curl(A)，但趋肤效应主要体现在电流 J
% J_tot = J_src - j*w*sigma*A
% 我们需要 A 的值。PostProcessor 目前只封装了 B 的探测。
% 我们可以扩展 PostProcessor，或者在这里手动插值 A。
% 为方便起见，我们利用 PostProcessor 的底层 Triangulation 查找单元，然后手动计算 A。

elem_center = pointLocation(post.Triangulation, pt_center);
pt_surface = [L/2*0.9, 0, 0.25];
elem_surface = pointLocation(post.Triangulation, pt_surface);

% 辅助函数: 计算单元电流密度
function J_val = compute_J_total(elem_idx, A_vec, J_src_vec, omega, sigma, mesh, dofHandler, spaceStr)
    dofMap = dofHandler.DofMaps(spaceStr);
    nodes = mesh.T(:, elem_idx);
    p_elem = mesh.P(:, nodes);
    dofs = dofMap(:, elem_idx);
    s = double(mesh.T2E_Sign(:, elem_idx));
    
    A_local = A_vec(dofs) .* s;
    
    % 简化的 P1 Nedelec 插值 (取中心值)
    q_center = [0.25; 0.25; 0.25];
    [val_ref, ~] = nedelec_tet_p1(q_center); 
    [J_mat, ~] = compute_jacobian_tet(p_elem);
    
    % Piola for H(curl): N = J^{-T} * N_ref
    % 注意: 这里我们需要 A 的值 (Value)，不是 Curl。
    % H(curl) 的基函数值变换规则是 covariant: v = J^{-T} * \hat{v}
    N_phy = (inv(J_mat)') * reshape(val_ref, 6, 3)'; 
    
    A_center = N_phy * A_local; % [3x1]
    
    % J = J_src - j*w*sigma*A
    J_val = J_src_vec(:) - 1j * omega * sigma * A_center;
end

J_c = compute_J_total(elem_center, A_phasor, SourceMap, omega, sigma, mesh, dofHandler, space.toString());
J_s = compute_J_total(elem_surface, A_phasor, SourceMap, omega, sigma, mesh, dofHandler, space.toString());

mag_J_c = abs(J_c(3));
mag_J_s = abs(J_s(3));

fprintf('\n--- Current Density Distribution (Jz) ---\n');
fprintf('  Center  |J| : %.2e A/m^2\n', mag_J_c);
fprintf('  Surface |J| : %.2e A/m^2\n', mag_J_s);
fprintf('  Ratio (Surf/Center): %.2f\n', mag_J_s / mag_J_c);

if mag_J_s > mag_J_c
    fprintf('\n[PASS] Skin effect observed: Current concentrates at the surface.\n');
else
    fprintf('\n[FAIL] No skin effect observed. Check parameters.\n');
end