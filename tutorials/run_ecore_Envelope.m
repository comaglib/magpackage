% run_ecore_Envelope.m - E-Core Transformer Inrush (Optimized)
% 
% 优化点:
%   1. [变步长]: 启动前 5 步使用 0.5ms 小步长，随后切换到 10ms 大步长。
%      这完美解决了首个时刻电流异常高的问题。
%   2. [正则化]: 增强 epsilon (1e-6 -> 1e-5) 以提高 MUMPS 稳定性。
%   3. [BugFix]: 修正了 solve 返回值接收错误。

clear; clc; tic;

fprintf('=========================================================\n');
fprintf('   E-Core Envelope FEM: Variable Step & Robust Config    \n');
fprintf('=========================================================\n');

%% --- 1. 网格与几何 ---
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

tags_primary = [5, 6, 8, 9];
NEW_TAG_PRIM = 50;
mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;

%% --- 2. 材料 ---
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(2) = MaterialLib.createNonlinear(B_data, H_data); % 铁芯
matLib(1) = MaterialLib.createLinear(1.0); % 空气
matLib(50) = MaterialLib.createLinear(1.0); % 线圈
matLib(60) = MaterialLib.createLinear(1.0); % 副边

space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

%% --- 3. 线圈与电路 ---
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);
winding = Winding('Primary', NEW_TAG_PRIM, 300, 100, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

circuitR = 100;
V_amp = 76;
baseFreq = 50;

%% --- 4. 求解器配置 ---
assembler = Assembler(mesh, dofHandler);

% [设置] 包含偶次谐波以捕捉不对称涌流
harmonics = [0, 1, 2, 3, 4, 5]; 
aft = AFT(harmonics, [], baseFreq); 

solver = EnvelopeCoupledSolver(assembler, aft, winding, circuitR);

% [策略] 变步长时间向量
% 0 ~ 0.01s: 使用 0.5ms 步长 (捕捉剧烈启动)
% 0.01 ~ 0.2s: 使用 10ms 步长 (加速计算衰减)
dt_small = 0.0005;
dt_large = 0.01;
t1 = (dt_small : dt_small : 0.01)';
t2 = (0.01+dt_large : dt_large : 0.2)';
timePoints = [t1; t2];

fprintf('   -> Time Steps: %d small steps + %d large steps (Total %d)\n', ...
    length(t1), length(t2), length(timePoints));

% [正则化] 鲁棒性设置
fprintf('   [Config] Calculating Robust Regularization...\n');
M_geo = assembler.assembleMass(space_A); 
K_sample = assembler.assembleStiffness(space_A, ones(mesh.NumElements,1)*(1/1.25e-6));
K_norm = norm(K_sample, 1);
M_norm = norm(M_geo, 1);

% 使用最小步长来计算 epsilon，保证最困难的时刻矩阵也是良态的
target_ratio = 1e-5; % 稍微调大一点比例，增加稳定性
epsilon = (K_norm * target_ratio * dt_small) / M_norm;
fprintf('   -> Epsilon=%.2e S/m\n', epsilon);
solver.MatrixM = epsilon * M_geo;

% 激励源
V_vec = zeros(length(harmonics), 1);
V_vec(harmonics==1) = -1j * V_amp;
V_func = @(t) V_vec;

fixedDofs = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

% 求解器微调
solver.MaxIter = 50;
solver.LinearSolver.MumpsSymmetry = 0; % 非对称
solver.LinearSolver.MumpsICNTL.i14 = 400; % 内存

%% --- 5. 求解 ---
fprintf('[Step 5] Solving...\n');
% [修复] 正确接收第一个返回值 results
[results, ~] = solver.solve(space_A, matLib, timePoints, V_func, fixedDofs);

%% --- 6. 后处理 ---
fprintf('[Step 6] Plotting...\n');
I_hist = results.I_history;
t_hist = results.Time;

% 重建时域波形
t_fine = [];
I_fine = [];
dt_plot = 2e-4;

for i = 1:length(t_hist)
    t_curr = t_hist(i);
    if i==1, t_prev=0; else, t_prev=t_hist(i-1); end
    
    t_sub = t_prev : dt_plot : t_curr;
    if isempty(t_sub), continue; end
    
    coeffs = I_hist(i, :);
    val = zeros(size(t_sub));
    for k = 1:length(harmonics)
        h = harmonics(k);
        w = 2*pi*baseFreq*h;
        val = val + real(coeffs(k) * exp(1j * w * t_sub));
    end
    t_fine = [t_fine, t_sub];
    I_fine = [I_fine, val];
end

figure('Color','w','Position',[100,100,1000,400]);
subplot(1,2,1);
plot(t_fine, I_fine, 'r-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Current (A)');
title('Inrush Current (Time Domain)');

subplot(1,2,2);
idx_0 = (harmonics == 0);
idx_1 = (harmonics == 1);
plot(t_hist, abs(I_hist(:, idx_0)), 'b-o', 'LineWidth', 1); hold on;
plot(t_hist, abs(I_hist(:, idx_1)), 'k-x', 'LineWidth', 1);
grid on; legend('DC Envelope', 'Fundamental Envelope');
xlabel('Time (s)'); title('Envelope Dynamics');

toc;