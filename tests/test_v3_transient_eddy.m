% test_v3_transient_eddy.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Transient Nonlinear Eddy Current (v3.0)         \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 (Mesh) ---
fprintf('[Step 1] Generating Mesh (6x6x6)...\n');
% 生成单位立方体网格
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,1,6));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); % 全域为 Region 1
mesh.generateEdges();

% --- 2. 自由度分发 (DoF) ---
fprintf('[Step 2] Distributing DoFs...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs == 0
    all_dofs = dofHandler.DofMaps(space.toString());
    dofHandler.NumGlobalDofs = max(all_dofs(:));
end
fprintf('         Total DoFs: %d\n', dofHandler.NumGlobalDofs);

% --- 3. 定义材料 (Rational Nonlinear Model) ---
fprintf('[Step 3] Defining Nonlinear Material...\n');
mu0 = 4*pi*1e-7;
mu_r_init = 2000; % 初始相对磁导率
B_sat = 1.5;      % 饱和磁感应强度

B_arr = linspace(0, 3.0, 100);
% 使用 Frohlich-Kennelly 模型保证光滑性
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);

MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

% --- 4. 组装器 (Assembler) ---
assembler = Assembler(mesh, dofHandler);

% --- 5. 设置瞬态参数 (TimeParams) ---
fprintf('[Step 4] Configuring Transient Solver...\n');

% 物理参数
f = 50;             % 频率 50 Hz
omega = 2*pi*f;
J_peak = 1e5;       % 峰值电流密度 100 kA/m^2
Sigma_val = 1e4;    % 电导率 S/m (设为中等值以观察明显的相位滞后，又不至于趋肤效应太强导致中心为0)

% 时间步进参数
N_cycles = 1.0;         % 模拟 1 个周期
T_period = 1/f;
dt = T_period / 20;     % 每个周期 20 步
N_steps = ceil(N_cycles * 20);

% 定义源函数 J(t)
% 输入: t (当前时间)
% 输出: [MaxTag x 3] 矩阵。这里只有 Region 1。
% J(t) = [0, 0, J_peak * sin(omega*t)]
TimeParams.SourceFunc = @(t) [0, 0, J_peak * sin(omega * t)];
TimeParams.dt = dt;
TimeParams.N_steps = N_steps;
TimeParams.Sigma = Sigma_val;

fprintf('         Freq: %d Hz, Conductivity: %.2e S/m\n', f, Sigma_val);
fprintf('         Steps: %d, dt: %.2e s\n', N_steps, dt);

% --- 6. 执行求解 ---
fprintf('[Step 5] Running Solver...\n');
[Solution, Info] = solve_transient_nonlinear(assembler, TimeParams, MatLibData, space);

% --- 7. 后处理与验证 ---
fprintf('[Step 6] Analyzing Results...\n');

% 选取中心单元进行观测
center = [0.5; 0.5; 0.5];
elem_centers = (mesh.P(:, mesh.T(1,:)) + mesh.P(:, mesh.T(2,:)) + ...
                mesh.P(:, mesh.T(3,:)) + mesh.P(:, mesh.T(4,:))) / 4;
dists = sum((elem_centers - center).^2, 1);
[~, center_elem_idx] = min(dists);

% 预计算几何因子
nodes = mesh.T(:, center_elem_idx);
p_elem = mesh.P(:, nodes);
[J_mat, detJ] = compute_jacobian_tet(p_elem);
q_pt = [0.25; 0.25; 0.25];
[~, curl_ref] = nedelec_tet_p1(q_pt);
C_ref = reshape(curl_ref, 6, 3);
C_phy = (C_ref * J_mat') / detJ;

% 提取数据历史
all_dofs_map = dofHandler.DofMaps(space.toString());
dofs = all_dofs_map(:, center_elem_idx);
s = double(mesh.T2E_Sign(:, center_elem_idx));

B_history = zeros(N_steps, 1);
J_history = zeros(N_steps, 1);
Time_axis = Solution.Time;

for n = 1:N_steps
    % 提取 B
    A_vec = Solution.A_history{n};
    A_local = A_vec(dofs) .* s;
    B_vec = C_phy' * A_local;
    B_history(n) = B_vec(3); % 取 Z 分量
    
    % 提取源 J
    src_val = TimeParams.SourceFunc(Time_axis(n));
    J_history(n) = src_val(3);
end

% 打印对比表
fprintf('\n--- Time History Data (Center Element) ---\n');
fprintf(' Step | Time (ms) | Source Jz (kA/m2) | Field Bz (T) \n');
fprintf('-----------------------------------------------------\n');
for n = 1:N_steps
    fprintf('  %2d  |   %5.2f   |     %8.2f      |   %8.4f   \n', ...
        n, Time_axis(n)*1000, J_history(n)/1000, B_history(n));
end

% 简单物理验证
% 1. 检查 B 是否随 J 变化 (相关性)
corr_coeff = corrcoef(J_history, B_history);
R_val = corr_coeff(1,2);

% 2. 检查滞后 (对于导电介质，B 应该滞后于 J)
% 找到峰值位置
[~, idx_J_max] = max(J_history);
[~, idx_B_max] = max(B_history);

fprintf('\n[Validation]\n');
fprintf('  - Correlation (J vs B): %.4f (Should be close to 1 but < 1 due to phase shift)\n', R_val);
fprintf('  - Peak J at step: %d\n', idx_J_max);
fprintf('  - Peak B at step: %d\n', idx_B_max);

if idx_B_max >= idx_J_max
    fprintf('  [PASS] Eddy current lag observed (B peak >= J peak).\n');
else
    fprintf('  [WARN] B peaks before J? Check coordinate system or sigma.\n');
end

if max(abs(B_history)) > 0.01
    fprintf('  [PASS] Significant magnetic field generated.\n');
else
    fprintf('  [FAIL] Magnetic field too weak.\n');
end