% run_ecore_breaker_energy_audit.m
% 
% 描述:
%   基于 BDF2BreakerSolver 的分闸与能量审计测试脚本。
%   用于验证:
%   1. 分闸逻辑: t > BreakerTime 后电流是否严格为 0。
%   2. 能量审计: 监测算法在强制电流截断时引入了多少数值耗散。
%   3. 剩磁锁定: 在无涡流情况下，分闸后的磁能是否保持守恒。

clear; clc; tic;
fprintf('=========================================================\n');
fprintf('   E-Core Breaker & Energy Audit Benchmark               \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
fprintf('[Step 1] Loading Mesh and Geometry...\n');
meshFile = 'data/Ecore.mphtxt'; % 请确保路径正确
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

% 区域 ID 重映射
tags_primary   = [5, 6, 8, 9];
tags_secondary = [3, 4, 7, 10];
tag_core       = 2;
tag_air        = 1;
NEW_TAG_PRIM = 50; 
NEW_TAG_SEC  = 60; 
mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料配置 ---
fprintf('[Step 2] Defining Materials and Physics...\n');
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 2.1 铁芯 (非线性)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 2.2 线性材料
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% % 2.3 电导率 (Sigma=0, 忽略涡流，聚焦于算法本身的耗散)
% sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
% all_tags = unique(mesh.RegionTags);
% for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

% 2.3 电导率设置
% [修改] 给铁芯赋予真实的电导率 (例如硅钢片等效电导率)
% 如果是实体铁芯，sigma 约为 1e6 ~ 1e7 S/m
% 如果是叠片铁芯，等效 sigma 较小，可设为 1e3 ~ 1e4 S/m 用于演示
% 空气和线圈区域依然为 0 (或线圈区也设为 0 以避免趋肤效应增加计算量)
sigma_core = 2e6;
sigmaMap(tag_core) = sigma_core;
sigmaMap(tag_air) = 0.0;
sigmaMap(NEW_TAG_PRIM) = 0.0;
sigmaMap(NEW_TAG_SEC) = 0.0;

% 2.4 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);
space_P = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_P);

%% --- 3. 线圈设置 ---
fprintf('[Step 3] Setup Field-Circuit Coupling...\n');
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

NP = 3000;
RP = 10; % [Ohm] 线圈电阻
winding = Winding('Primary', NEW_TAG_PRIM, NP, RP, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

% 电路参数
circuit = struct();
circuit.R = RP;
circuit.L = 0;
Vamp = 1500;
circuit.V_source_func = @(t) Vamp * sin(2 * pi * 50 * t);

%% --- 4. 求解器配置 (使用 BDF2BreakerSolver) ---
fprintf('[Step 4] Configuring Breaker Solver...\n');
assembler = Assembler(mesh, dofHandler);

% >>>>>>>>> 关键修改：实例化 Breaker Solver <<<<<<<<<
solver = BDF2BreakerSolver(assembler);

% 4.1 设置断路器参数
solver.BreakerTime = 0.01;   % [s] 在 xx s分闸
solver.R_winding_Audit = RP; % [Ohm] 必须与物理电阻一致，用于计算 ohm loss

% 4.2 求解器参数
solver.Tolerance = 1e-4;
solver.RelTolerance = 1e-4;
solver.LinearSolver.MumpsSymmetry = 0;
solver.LinearSolver.MumpsICNTL.i14 = 300;

% 4.3 时间步进
dt = 5e-5;
timeSim = 0.015;    % 总时长
% timeSteps = repmat(dt, round(timeSim/dt), 1);
timeSteps = [repmat(5e-4, round(solver.BreakerTime/(5e-4)), 1);
             repmat(dt, round((timeSim-solver.BreakerTime)/dt), 1)];


% 4.4 边界条件
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);

% 4.5 监视器 (使用铁芯中心作为探测点)
probePoint = [0, 0, 0]; 
plotFunc = @(t, I, t_vec, I_vec, B_curr, B_hist) monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist, solver.BreakerTime);

%% --- 5. 求解 ---
fprintf('[Step 5] Solving with Breaker Logic...\n');
% 注意：solver.solve 方法会自动调用 plotEnergyAudit 绘制审计图
[sols, info] = solver.solve(space_A, space_P, matLib, sigmaMap, ...
                            circuit, winding, timeSteps, ...
                            fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);

%% --- 6. 结果分析 (Finalized Energy Audit Logic) ---
fprintf('[Step 6] Analysis Complete.\n');

% 1. 寻找分闸后的第一个时刻 (First OPEN step)
%    原理：Solver 中使用了 currentTime > (BreakerTime + 1e-9)
%    因此，idx_first_open 对应的是电流强制归零的第一步 (t = 0.01 + dt)
idx_first_open = find(info.times > (solver.BreakerTime + 1e-9), 1);

if ~isempty(idx_first_open) && idx_first_open > 1
    % 2. 确定"分闸瞬间"的基准点 (The Exact Break Time, CLOSED state)
    %    这是 t = 0.01 时刻，此时电流尚未归零，磁场处于满载状态
    idx_break_moment = idx_first_open - 1; 
    
    % 打印状态检查，确保时间对齐逻辑正确
    fprintf('\n--------------------------------------------------\n');
    fprintf('   STATUS CHECK (Time Alignment Verification)\n');
    fprintf('--------------------------------------------------\n');
    fprintf('   T_break     (idx=%d): t=%.6fs | I = %.4f A (CLOSED, Max Energy)\n', ...
            idx_break_moment, info.times(idx_break_moment), info.CurrentHistory(idx_break_moment));
    fprintf('   T_break+dt  (idx=%d): t=%.6fs | I = %.4f A (OPEN,   Current Cut)\n', ...
            idx_first_open, info.times(idx_first_open), info.CurrentHistory(idx_first_open));
    
    % 3. 定义积分区间 (从分闸瞬间一直到仿真结束)
    idx_span = idx_break_moment : length(info.times);
    t_span = info.times(idx_span);
    
    % 4. 提取能量数据
    %    总释放磁能 (Source) = 分闸前一刻磁能 - 最终剩余磁能
    W_mag_before = info.EnergyLog.W_mag(idx_break_moment); 
    W_mag_end    = info.EnergyLog.W_mag(end);              
    Delta_W_mag  = W_mag_before - W_mag_end;
    
    % 5. 积分计算物理损耗与数值误差 (Sinks)
    %    使用 trapz (梯形积分) 计算累计能量 (J)
    P_eddy_vec = info.EnergyLog.P_eddy(idx_span);
    P_num_vec  = info.EnergyLog.P_num(idx_span);
    
    E_eddy_total = trapz(t_span, P_eddy_vec); % 物理耗散 (涡流热)
    E_num_total  = trapz(t_span, P_num_vec);  % 数值耗散 (算法误差)
    
    % =====================================================================
    % [核心修改] 评价指标归一化逻辑优化
    % =====================================================================
    % 分母基准：总释放的磁场能量 (Total Source Energy)
    total_source = abs(Delta_W_mag) + 1e-12; % 避免除零
    
    % 指标 1: 物理捕捉率 (有多少能量被识别为涡流热?)
    ratio_phys = (E_eddy_total / total_source) * 100;
    
    % 指标 2: 数值耗散率 (有多少能量被算法吃掉了?)
    ratio_num = (E_num_total / total_source) * 100;
    
    % 6. 打印详细审计报表
    fprintf('\n==================================================\n');
    fprintf('   DETAILED ENERGY AUDIT (Post-Breaker Phase)     \n');
    fprintf('==================================================\n');
    fprintf('1. Magnetic Energy Release (Source):  %10.6f J\n', Delta_W_mag);
    fprintf('   (The "Pie" to be divided)\n');
    fprintf('--------------------------------------------------\n');
    fprintf('2. Physical Eddy Heat (Captured):     %10.6f J\n', E_eddy_total);
    fprintf('   [Ratio]: %6.2f %% (Physics Fidelity)\n', ratio_phys);
    fprintf('--------------------------------------------------\n');
    fprintf('3. Numerical Dissipation (Lost):      %10.6f J\n', E_num_total);
    fprintf('   [Ratio]: %6.2f %% (Numerical Error)\n', ratio_num);
    fprintf('==================================================\n');
    fprintf('   Check Sum: %.2f %% (Should be ~100%%)\n', ratio_phys + ratio_num);
    
    % 7. 自动评级
    if ratio_num < 5.0
        fprintf('   [PASS] High Fidelity. Numerical loss is negligible.\n');
    elseif ratio_num < 20.0
        fprintf('   [ACCEPTABLE] Moderate Fidelity. Main physics captured.\n');
    else
        fprintf('   [WARNING] Low Fidelity (Undersampled). Algorithm dominates energy flow.\n');
        fprintf('             Action: Decrease time step (dt) to capture Eddy peaks.\n');
    end
    
    % ---------------------------------------------------------------------
    % 绘制分闸瞬间的功率平衡细节 (Zoom-in Plot)
    % ---------------------------------------------------------------------
    figure('Name', 'Power Balance Zoom (Breaker Moment)', 'Position', [150, 150, 700, 500]);
    
    % 定义显示窗口: 分闸前 5 步 ~ 分闸后 25 步
    idx_start_zoom = max(1, idx_break_moment - 5);
    idx_end_zoom   = min(length(info.times), idx_break_moment + 25);
    range_zoom = idx_start_zoom : idx_end_zoom;
    
    t_zoom = info.times(range_zoom);
    P_eddy_zoom = info.EnergyLog.P_eddy(range_zoom);
    
    % 计算瞬时磁场释放功率 P_mag = -dW/dt (使用差分估算)
    W_zoom = info.EnergyLog.W_mag(range_zoom);
    % 为了对齐，使用中心差分或简单的后向差分
    P_mag_release = -gradient(W_zoom, t_zoom); 
    
    % 绘制曲线
    plot(t_zoom, P_mag_release, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Magnetic Release (-dW/dt)');
    hold on;
    plot(t_zoom, P_eddy_zoom, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Physical Eddy Heat (P_{eddy})');
    
    % 绘制误差区域 (填充图)
    % 误差 = 释放 - 消耗。如果在某时刻黑线远高于绿线，中间的区域就是数值耗散。
    % 这里处理一下数据以防止绘图交叉
    fill_y_top = P_mag_release;
    fill_y_bot = P_eddy_zoom;
    
    fill_x = [t_zoom; flipud(t_zoom)];
    fill_y = [fill_y_top; flipud(fill_y_bot)];
    fill(fill_x, fill_y, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Numerical Dissipation Area');
    
    xline(solver.BreakerTime, 'r--', 'Breaker Trigger', 'LineWidth', 1.2);
    
    grid on; legend('Location', 'Best');
    title('Instantaneous Power Balance @ Break');
    xlabel('Time (s)'); ylabel('Power (W)');
    subtitle(sprintf('Gap between Black and Green lines indicates Numerical Error (dt=%.1e)', t_span(2)-t_span(1)));
end
toc;

%% --- 辅助绘图函数 ---
function monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist, t_break)
    persistent hFig;
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Real-time Monitor', 'Position', [50, 50, 700, 500]);
    end
    set(0, 'CurrentFigure', hFig);
    
    % 电流
    subplot(2, 1, 1);
    plot(t_vec, I_vec, 'r-', 'LineWidth', 1.5);
    xline(t_break, 'k--', 'Breaker Open');
    grid on; ylabel('Current (A)');
    title(sprintf('Time: %.4fs | I: %.2fA', t, I));
    xlim([0, max(t_vec(end), 0.01)]);
    
    % 磁密
    subplot(2, 1, 2);
    if isempty(B_hist), B_plot = 0; t_plot = 0; else, B_plot = B_hist; t_plot = t_vec; end
    plot(t_plot, B_plot, 'b-', 'LineWidth', 1.5);
    xline(t_break, 'k--', 'Breaker Open');
    grid on; ylabel('|B| (T)'); xlabel('Time (s)');
    title(sprintf('Probe |B|: %.4f T', B_curr));
    xlim([0, max(t_vec(end), 0.01)]);
    
    drawnow;
end