% run_ecore_soft_switching_audit.m
% 
% 描述:
%   基于 BDF2BreakerSolver 的 "软关断" 能量审计验证脚本。
%   
%   测试目标:
%   验证在物理过程平滑（非狄拉克脉冲）的情况下，算法是否能实现完美的能量守恒。
%   如果此案例通过，说明算法本身没有问题，之前的误差纯粹源于"硬开关"的奇点性质。
%   
%   设置:
%   1. Circuit: t > 0.01s 后，电压线性斜坡下降 (Soft Ramp Down)。
%   2. Solver:  禁用强制断路逻辑 (BreakerTime = Inf)，允许自然衰减。
%   3. Sigma:   Effective Sigma (5e5 S/m) 用于产生物理阻尼。

clear; clc; tic;
fprintf('=========================================================\n');
fprintf('   E-Core Soft Switching & Energy Audit Benchmark        \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
fprintf('[Step 1] Loading Mesh and Geometry...\n');
meshFile = 'data/Ecore.mphtxt'; 
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

% 区域重映射
tags_primary   = [5, 6, 8, 9];
tags_secondary = [3, 4, 7, 10];
tag_core       = 2;
tag_air        = 1;
NEW_TAG_PRIM = 50; 
NEW_TAG_SEC  = 60; 

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料配置 (Effective Sigma Mode) ---
fprintf('[Step 2] Defining Materials and Physics...\n');
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 2.1 铁芯 (非线性)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];

mat_iron = MaterialLib.createNonlinear(B_data, H_data);
% 注入 Bertotti 参数 (仅用于对比，主要依靠 Solver 损耗)
if ismethod(mat_iron, 'setBertottiParams')
    mat_iron = MaterialLib.setBertottiParams(mat_iron, 0.35e-3, 2.0e6, 0, 2.0, 0); 
else
    mat_iron.Lam_d = 0.35e-3; mat_iron.Lam_sigma = 2.0e6;
end
matLib(tag_core) = mat_iron;

% 2.2 线性材料
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% 2.3 电导率设置
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');

% [关键] 铁芯设为工程等效值 (5e5 S/m)，产生物理阻尼
sigmaMap(tag_core) = 5e5; 
% 线圈设为 0，避免趋肤效应干扰
sigmaMap(NEW_TAG_PRIM) = 0.0; 
sigmaMap(tag_air) = 0.0;
sigmaMap(NEW_TAG_SEC) = 0.0;

% 2.4 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);
space_P = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_P);

%% --- 3. 线圈与软关断电路设置 ---
fprintf('[Step 3] Setup Field-Circuit Coupling (Soft Switching)...\n');
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

NP = 3000;
RP = 10;
winding = Winding('Primary', NEW_TAG_PRIM, NP, RP, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

% [关键] 定义软关断电压源
circuit = struct();
circuit.R = RP;
circuit.L = 0;
Vamp = 1500;

% 定义参数：0.01s 开始关断，耗时 1ms (0.001s)
T_trigger = 0.01;
T_ramp = 0.001; 

circuit.V_source_func = @(t) voltage_ramp(t, Vamp, T_trigger, T_ramp);

%% --- 4. 求解器配置 ---
fprintf('[Step 4] Configuring Solver...\n');
assembler = Assembler(mesh, dofHandler);
solver = BDF2BreakerSolver(assembler);

% [关键] 禁用强制断路逻辑
% 设置为一个不可能达到的时间，让求解器始终处于 "CLOSED" 模式
% 这样电流就会完全由电压源 V(t) 和阻抗决定
solver.BreakerTime = 100.0;   
solver.R_winding_Audit = RP; 

% 求解器参数
solver.Tolerance = 1e-4;
solver.RelTolerance = 1e-4;
solver.LinearSolver.MumpsSymmetry = 0;
solver.LinearSolver.MumpsICNTL.i14 = 300;

% 4.3 时间步进 (均匀步长即可，因为过程是平滑的)
dt = 2e-4;
timeSim = 0.015; % 总时长 15ms
timeSteps = repmat(dt, round(timeSim/dt), 1);

% 边界条件
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);

% 监视器
probePoint = [0, 0, 0]; 
plotFunc = @(t, I, t_vec, I_vec, B_curr, B_hist) monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist, T_trigger);

%% --- 5. 求解 ---
fprintf('[Step 5] Solving (Soft Switching Mode)...\n');
[sols, info] = solver.solve(space_A, space_P, matLib, sigmaMap, ...
                            circuit, winding, timeSteps, ...
                            fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);

%% --- 6. 结果分析 (Energy Audit) ---
fprintf('[Step 6] Analysis Complete.\n');

% 1. 确定分析区间 (从斜坡开始到结束)
idx_trigger = find(info.times >= T_trigger, 1);
if isempty(idx_trigger), idx_trigger = 1; end

% 分析从触发前一点点开始，直到仿真结束
idx_span = (idx_trigger - 5) : length(info.times);
if idx_span(1) < 1, idx_span = 1:length(info.times); end
t_span = info.times(idx_span);

% 2. 提取累计能量 (使用 solve 函数中计算好的累积值)
%    注意：solve 中是累加的，我们需要取区间内的增量
%    但为了简单，我们取整个过程的平衡，或者利用功率积分

P_source_vec = zeros(size(t_span));
% 重算源功率: P = V(t) * I(t)
I_vec = info.CurrentHistory(idx_span);
for k=1:length(t_span)
    V_val = circuit.V_source_func(t_span(k));
    P_source_vec(k) = V_val * I_vec(k);
end

P_eddy_vec = info.EnergyLog.P_eddy(idx_span); % 包含 Coil Ohm + Core Eddy
P_num_vec  = info.EnergyLog.P_num(idx_span);

% 磁能变化 dW
W_vec = info.EnergyLog.W_mag(idx_span);
Delta_W = W_vec(end) - W_vec(1);

% 积分能量 (J)
E_input_total = trapz(t_span, P_source_vec); % 输入的总能量 (可能是负的，如果电感在放电回馈电源)
% 注意：在关断阶段，电压下降，I依然为正，P=VI > 0，电源依然在做功(或者变小)。
% 这里的审计更关注平衡方程： Input = Delta_W + E_loss + E_error

E_loss_total = trapz(t_span, P_eddy_vec);
E_num_total  = trapz(t_span, P_num_vec);

% 3. 平衡检查
% Balance: E_input - Delta_W - E_loss = E_num
LHS = E_input_total - Delta_W - E_loss_total;

% 归一化基准：取各项绝对值的最大值，避免分母太小
denom = max([abs(E_input_total), abs(Delta_W), abs(E_loss_total)]);
ratio_error = (abs(E_num_total) / denom) * 100;

fprintf('\n==================================================\n');
fprintf('   SOFT SWITCHING ENERGY AUDIT (Whole Window)     \n');
fprintf('==================================================\n');
fprintf('Time Span: %.4fs -> %.4fs\n', t_span(1), t_span(end));
fprintf('--------------------------------------------------\n');
fprintf('1. Net Energy Input (Source):     %+10.6f J\n', E_input_total);
fprintf('2. Change in Mag Energy (dW):     %+10.6f J\n', Delta_W);
fprintf('3. Physical Loss (Ohm+Eddy):      %+10.6f J\n', E_loss_total);
fprintf('--------------------------------------------------\n');
fprintf('4. Numerical Error (Residual):    %+10.6f J\n', E_num_total);
fprintf('   [Relative Error]: %.4f %%\n', ratio_error);
fprintf('==================================================\n');

% 4. 绘图
figure('Name', 'Soft Switching Power Balance', 'Position', [150, 150, 800, 600]);

% 计算 dW/dt
P_mag_change = gradient(W_vec, t_span); % +dW/dt

subplot(2,1,1);
plot(t_span, P_source_vec, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Source Input (P_{in})');
hold on;
plot(t_span, P_mag_change + P_eddy_vec, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Total Consumption (dW/dt + P_{loss})');
xline(T_trigger, 'k--', 'Ramp Start');
xline(T_trigger + T_ramp, 'k-.', 'Ramp End');
grid on; legend('Location', 'Best');
ylabel('Power (W)'); title('Power Supply vs Consumption');

subplot(2,1,2);
plot(t_span, P_num_vec, 'r-', 'LineWidth', 1.5);
grid on; ylabel('Error Power (W)'); xlabel('Time (s)');
title('Numerical Residual (P_{num})');
subtitle(sprintf('Integrated Error: %.2e J (%.2f%%)', E_num_total, ratio_error));

toc;

%% --- 局部函数 ---

function V = voltage_ramp(t, Vamp, t_trigger, t_ramp)
    % 软关断电压函数
    if t <= t_trigger
        V = Vamp * sin(2 * pi * 50 * t);
    elseif t <= t_trigger + t_ramp
        % 线性衰减系数 (1 -> 0)
        scale = 1.0 - (t - t_trigger) / t_ramp;
        % 保持波形连续性：这里简单地让正弦波振幅衰减
        % 或者更简单：直接让当前时刻的电压值线性降到0
        % 为了物理平滑，我们采用 "V_trigger * scale" 策略
        
        % 策略A: 振幅衰减 (保持交流特性) -> V = (Vamp * scale) * sin(...)
        % 策略B: 直流式衰减 (模拟开关拉弧电压降) -> 从 V(t_trig) 线性降到 0
        
        % 这里使用策略B，更像分闸
        V_trigger_val = Vamp * sin(2 * pi * 50 * t_trigger);
        V = V_trigger_val * scale;
    else
        V = 0;
    end
end

function monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist, t_break)
    persistent hFig;
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Real-time Monitor', 'Position', [50, 50, 700, 500]);
    end
    set(0, 'CurrentFigure', hFig);
    
    subplot(2, 1, 1);
    plot(t_vec, I_vec, 'r.-', 'LineWidth', 1);
    xline(t_break, 'k--');
    grid on; ylabel('Current (A)');
    title(sprintf('Time: %.5fs | I: %.2fA', t, I));
    
    subplot(2, 1, 2);
    if isempty(B_hist), B_plot = 0; t_plot = 0; else, B_plot = B_hist; t_plot = t_vec; end
    plot(t_plot, B_plot, 'b.-', 'LineWidth', 1);
    xline(t_break, 'k--');
    grid on; ylabel('|B| (T)'); xlabel('Time (s)');
    title(sprintf('Probe |B|: %.4f T', B_curr));
    drawnow limitrate;
end