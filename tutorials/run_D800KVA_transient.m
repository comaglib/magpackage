% run_D800KVA_compare.m
% =========================================================================
% D800KVA 变压器瞬态分析 (双绕组并行 BDF2 版本)
% 
% 物理工况：
%   - 高压侧 (Primary) 双线圈并联励磁
%   - 电压 49497V, 50Hz, 初始相位 0
%   - 低压侧开路
%   - 铁芯非线性 B-H 曲线
% =========================================================================

clear; clc;
fprintf('=========================================================\n');
fprintf('   D800KVA Transient Analysis (Multi-Winding BDF2)       \n');
fprintf('=========================================================\n');

%% --- 1. 模型设置 (Mesh & Physics) ---
fprintf('[Setup] Loading Mesh and Materials...\n');

% 1.1 加载网格
meshFile = 'data/D800KVA.mphtxt'; 
if ~exist(meshFile, 'file')
    error('Mesh file %s not found.', meshFile);
end
mesh = Mesh.load(meshFile, 'm');
mesh.generateEdges();

% 1.2 区域 ID 定义
tags_p1 = [5, 6, 10, 13];    % Primary 1
tags_p2 = [14, 15, 19, 22];  % Primary 2
tags_s1 = [7, 8, 11, 12];    % Secondary 1
tags_s2 = [16, 17, 20, 21];  % Secondary 2
tag_core = [2, 3, 4, 9, 18, 23];
tag_air = 1;

% 区域标签重映射 (保持独立的线圈 ID)
NEW_TAG_PRIM1 = 50; 
NEW_TAG_PRIM2 = 51; 
NEW_TAG_SEC1 = 60;
NEW_TAG_SEC2 = 61;

mesh.RegionTags(ismember(mesh.RegionTags, tags_p1)) = NEW_TAG_PRIM1;
mesh.RegionTags(ismember(mesh.RegionTags, tags_p2)) = NEW_TAG_PRIM2;
mesh.RegionTags(ismember(mesh.RegionTags, tags_s1)) = NEW_TAG_SEC1;
mesh.RegionTags(ismember(mesh.RegionTags, tags_s2)) = NEW_TAG_SEC2;

% 1.3 材料定义
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% B-H 曲线数据
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, ...
          1, 1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118163704, 2.23];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, ...
          52.1, 56.77, 61.67, 67.19, 74.03, 83.86, 100.17, 137.5, 295.39, 1014.01, ...
          2500, 12500, 112500, 212500];

% 创建非线性材料 (铁芯)
matNonlinear = MaterialLib.createNonlinear(B_data, H_data);
for t = tag_core
    matLib(t) = matNonlinear;
end

% 线性材料 (空气和线圈)
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air) = matLinear;
for t = [NEW_TAG_PRIM1, NEW_TAG_PRIM2, NEW_TAG_SEC1, NEW_TAG_SEC2]
    matLib(t) = matLinear; 
end

% 电导率设置 (无涡流)
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

% 1.4 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

space_P = FunctionSpace('Lagrange', 1); 
dofHandler.distributeDofs(space_P);

% 1.5 线圈绕组配置 (使用双绕组独立定义)
N_turns_single = 2933;
R_single = 23.42;
V_amp_single = 49497;

fprintf('[Setup] Computing coil directions...\n');
% 计算 Primary 1 方向
[c1, ~, S1, ax1] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM1);
dir_map1 = CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM1, c1, ax1);

% 计算 Primary 2 方向
[c2, ~, S2, ax2] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM2);
dir_map2 = CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM2, c2, ax2);

% Primary 2 物理反接处理
dir_map2 = -dir_map2; 

% 创建两个独立的 Winding 对象
w1 = Winding('Primary_1', NEW_TAG_PRIM1, N_turns_single, R_single, S1, [0,0,0]);
w1.setDirectionField(dir_map1);

w2 = Winding('Primary_2', NEW_TAG_PRIM2, N_turns_single, R_single, S2, [0,0,0]);
w2.setDirectionField(dir_map2);

windings = [w1, w2]; % 构造绕组数组

% 1.6 外电路参数 (对应两个绕组)
% V_source = V_amp * sin(wt) (并行励磁)
V_func = @(t) V_amp_single * sin(2 * pi * 50 * t);

% 电路 1
c1 = struct();
c1.R = R_single; c1.L = 0;
c1.V_source_func = V_func;

% 电路 2
c2 = struct();
c2.R = R_single; c2.L = 0;
c2.V_source_func = V_func;

circuits = [c1, c2]; % 构造电路数组

% 1.7 边界条件与组装器
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);
assembler = Assembler(mesh, dofHandler);

% 1.8 探测点
probePoint = [-0.25, 0, 0]; 

% 1.9 绘图回调 (修改为调用本地实时绘图函数)
plotFunc = @(t, I_vals, t_vec, I_hist, B_curr, B_hist) ...
    monitor_callback(t, I_vals, t_vec, I_hist, B_curr, B_hist);

%% --- 2. 运行 BDF2 求解器 ---
fprintf('\n[Run] Starting TransientBDF2Solver (Multi-Winding)...\n');
dt_bdf = 5e-4;       % 步长
timeSim = 0.02;      % 仿真时长 20ms
timeSteps_bdf = repmat(dt_bdf, round(timeSim/dt_bdf), 1);

solver_bdf = TransientBDF2Solver(assembler);
solver_bdf.Tolerance = 1e-3; 
solver_bdf.RelTolerance = 1e-3;

tic;
% 传入绕组数组和电路数组
[~, info_bdf] = solver_bdf.solve(space_A, space_P, matLib, sigmaMap, ...
                                 circuits, windings, timeSteps_bdf, ...
                                 fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);
time_bdf = toc;
fprintf('\n-> BDF2 Completed in %.2f seconds.\n', time_bdf);

%% --- 3. 结果对比绘图 ---
fprintf('\n[Post] Plotting results...\n');

t_vec = cumsum(timeSteps_bdf);
I_history = info_bdf.CurrentHistory; % [Step x 2] 矩阵

I1 = I_history(:, 1);
I2 = I_history(:, 2);
I_total = I1 + I2; % 计算总电流

B_hist = info_bdf.ProbeB_History;

figure('Name', 'D800KVA Multi-Winding Result', 'Position', [200, 200, 1000, 600]);

% 子图 1: 电流
subplot(2, 1, 1);
plot(t_vec*1e3, I_total, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Total Current');
hold on;
plot(t_vec*1e3, I1, 'b--', 'LineWidth', 1.0, 'DisplayName', 'Branch 1');
plot(t_vec*1e3, I2, 'g--', 'LineWidth', 1.0, 'DisplayName', 'Branch 2');
grid on;
legend('Location', 'best');
xlabel('Time (ms)'); ylabel('Current (A)');
title('D800KVA Inrush Current (Parallel Excitation)');

% 子图 2: 磁密
subplot(2, 1, 2);
plot(t_vec*1e3, B_hist, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)'); ylabel('|B| (T)');
title(sprintf('Magnetic Flux Density at [%.2f, %.2f, %.2f]', probePoint));

fprintf('\nDone.\n');

%% --- 辅助局部函数 ---
function monitor_callback(t, I_vals, t_vec, I_hist, B_curr, B_hist)
    % I_vals: 当前时刻各绕组电流向量
    % I_hist: 历史电流矩阵 [Step x NumWindings]
    
    persistent hFig;
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Real-time Monitor', 'Position', [100, 100, 800, 600]);
    end
    set(0, 'CurrentFigure', hFig);

    % --- 子图 1: 多绕组电流波形 ---
    subplot(2, 1, 1);
    
    % 绘制各个分支电流
    plot(t_vec, I_hist(:, 1), 'b-', 'LineWidth', 1.0, 'DisplayName', 'I_1'); 
    hold on;
    if size(I_hist, 2) >= 2
        plot(t_vec, I_hist(:, 2), 'g--', 'LineWidth', 1.0, 'DisplayName', 'I_2');
    end
    
    % 绘制总电流 (如果有多于1个绕组)
    if size(I_hist, 2) > 1
        I_total = sum(I_hist, 2);
        plot(t_vec, I_total, 'r:', 'LineWidth', 1.5, 'DisplayName', 'I_{total}');
    end
    hold off;
    
    grid on;
    legend('Location', 'best');
    ylabel('Current (A)', 'FontSize', 10);
    
    % 标题显示当前数值
    if length(I_vals) >= 2
        title(sprintf('t=%.4fs | I1=%.2fA | I2=%.2fA', t, I_vals(1), I_vals(2)), 'FontSize', 11);
    else
        title(sprintf('t=%.4fs | I=%.2fA', t, I_vals(1)), 'FontSize', 11);
    end
    
    xlim([0, max(t_vec(end), 1e-6)]); % 动态调整轴范围
    
    % --- 子图 2: 磁密波形 ---
    subplot(2, 1, 2);
    if isempty(B_hist)
        B_plot = 0; t_plot = 0;
    else
        B_plot = B_hist; t_plot = t_vec;
    end
    plot(t_plot, B_plot, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 3);
    grid on;
    ylabel('|B| (T)', 'FontSize', 10);
    xlabel('Time (s)', 'FontSize', 10);
    title(sprintf('Probe |B|: %.4f T', B_curr), 'FontSize', 11);
    xlim([0, max(t_vec(end), 1e-6)]);
    
    drawnow; 
end