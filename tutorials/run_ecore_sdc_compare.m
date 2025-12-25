% run_ecore_sdc_compare.m
% =========================================================================
% E-Core 变压器瞬态分析对比测试脚本
% 对比对象：
%   1. 标准 BDF2 求解器 (TransientCoupledSolver) - 作为基准 (Reference)
%   2. 高阶 SDC 求解器 (SDCSolver)               - 待验证的高阶算法
%
% 测试目标：
%   验证 SDC 算法在粗时间步长下，能否达到与细步长 BDF2 相当甚至更高的精度，
%   特别是针对变压器涌流这种强非线性刚性问题。
% =========================================================================

clear; clc;
fprintf('=========================================================\n');
fprintf('   Comparison: Standard BDF2 vs High-Order SDC           \n');
fprintf('=========================================================\n');

%% --- 1. 模型共享设置 (Mesh & Physics) ---
% 目的：确保两个求解器使用完全相同的网格、材料和电路参数，保证对比公平性。

fprintf('[Setup] Loading Mesh and Materials...\n');

% 1.1 加载网格
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found: %s', meshFile); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

% 1.2 区域重映射 (Re-tagging)
% 原网格标签分散，需统一为：空气(1), 铁芯(2), 原边线圈(50), 副边线圈(60)
tags_primary = [5, 6, 8, 9]; 
tags_secondary = [3, 4, 7, 10];
tag_core = 2; 
tag_air = 1;
NEW_TAG_PRIM = 50; 
NEW_TAG_SEC = 60;

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

% 1.3 材料定义 (MaterialLib v4.0 PCHIP Smoothness)
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% B-H 曲线数据 (DW465-50 硅钢片)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164, 3.37480104];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500, 1112500];

% 创建非线性材料 (自动应用 PCHIP 平滑和真空修正)
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 线性材料 (空气和线圈)
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air) = matLinear; 
matLib(NEW_TAG_PRIM) = matLinear; 
matLib(NEW_TAG_SEC) = matLinear;

% 电导率设置 (本例忽略涡流，设为 0)
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

% 1.4 有限元空间与自由度分发
space_A = FunctionSpace('Nedelec', 1); % 磁矢位 A (边缘元)
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

space_P = FunctionSpace('Lagrange', 1); % 标量电位 P (节点元，用于规范)
dofHandler.distributeDofs(space_P);

% 1.5 线圈几何参数自动提取
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);

R = 8;
winding = Winding('Primary', NEW_TAG_PRIM, 3000, R, area_S, [0,0,0]);
winding.setDirectionField(dir_map);

% 1.6 外电路参数
voltage = 1500;
circuit = struct();
circuit.R = R;  % 电阻 (Ohm)
circuit.L = 0;  % 漏电感 (H)
% 电压源: voltage, 50Hz 正弦波
% circuit.V_source_func = @(t) voltage * sin(2 * pi * 50 * t);
% 定义分段激励函数
% 物理含义：50Hz 正弦波，幅值 voltage，在 t=0.01s 时刻切断（归零）
timeDisc = 100;
circuit.V_source_func = @(t) voltage * sin(2 * pi * 50 * t) .* (t < timeDisc);

% 1.7 边界条件与组装器
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);
assembler = Assembler(mesh, dofHandler);

% 探针位置 (铁芯中心)
probePoint = [0, 0, 0]; 

% %% --- 2. 运行 SDC (High-Order) ---
% fprintf('\n[Run 2] High-Order SDC Solver (Coarse Step)...\n');
% 
% % SDC 策略: 
% timeSim = 0.013;     % 仿真时长
% dt_slab = 3.5e-3;    % 4.5e-3; 3.5e-3;
% timeSteps_sdc = [dt_slab,timeSim-dt_slab];
% 
% solver_sdc = AcISDCSolver(assembler);
% solver_sdc.PolyOrder = 4;       % 多项式阶数 P=4 (5个节点)
% solver_sdc.MaxSDCIters = 10;    % SDC 最大修正次数
% solver_sdc.SDCTolerance = 1e-3; % SDC 收敛容差 (相对/绝对混合判据)
% 
% tic;
% [~, info_sdc] = solver_sdc.solve(space_A, space_P, matLib, sigmaMap, ...
%                                  circuit, winding, timeSteps_sdc, ...
%                                  fixedDofs_A, fixedDofs_P, [], [], probePoint);
% time_sdc = toc;
% fprintf('-> SDC Completed in %.2f seconds.\n', time_sdc);

%% --- 3. 运行 BDF2 (Reference) ---
fprintf('\n[Run 1] Standard BDF2 Solver (Fine Step)...\n');

dt_bdf = 5e-3;
timeSim = 0.02*5;
timeSteps_bdf = repmat(dt_bdf, round(timeSim/dt_bdf), 1);
% timeSteps_bdf = [repmat(5e-4, round(timeDisc/(5e-4)), 1);
%                  repmat(dt_bdf, round((timeSim-timeDisc)/dt_bdf), 1)];

tic;
solver_bdf = TransientBDF2Solver(assembler);
solver_bdf.Tolerance = 1e-4;
solver_bdf.RelTolerance = 1e-4;

% 定义绘图回调函数句柄
% 注意: 求解器会自动传入 (t, I, t_vec, I_vec, B_curr, B_hist)
plotFunc = @(t, I, t_vec, I_vec, B_curr, B_hist) monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist);

[~, info_bdf] = solver_bdf.solve(space_A, space_P, matLib, sigmaMap, ...
                                 circuit, winding, timeSteps_bdf, ...
                                 fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);
time_bdf = toc;
fprintf('-> BDF2 Completed in %.2f seconds.\n', time_bdf);

% %% --- 4. 结果对比绘图 ---
% fprintf('\n[Post] Plotting comparison...\n');
% 
% % BDF2 (Reference)
% % t_bdf = cumsum(timeSteps_bdf);
% % I_bdf = info_bdf.CurrentHistory;
% % B_bdf = info_bdf.ProbeB_History;
% 
% % SDC (High-Res Interpolated)
% t_sdc_full = info_sdc.Time_Full;
% I_sdc_full = info_sdc.Current_Full;
% B_sdc_full = info_sdc.ProbeB_Full;
% 
% figure('Name', 'BDF2 vs SDC Comparison', 'Position', [200, 200, 1000, 600]);
% 
% % 子图 1: 电流对比
% subplot(2, 1, 1);
% % 绘制 BDF2 参考线
% % plot(t_bdf*1e3, I_bdf, 'k-', 'LineWidth', 1.5, 'DisplayName', 'BDF2 (Ref)');
% hold on;
% % 绘制 SDC 全波形 (插值后的光滑曲线)
% plot(t_sdc_full*1e3, I_sdc_full, 'r-', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('SDC (P=%d, Polynomial Fit)', solver_sdc.PolyOrder));
% 
% % 标记真实的 GLL 计算节点，展示 SDC 是如何"以少胜多"的
% % t_sdc_ends = cumsum(timeSteps_sdc);
% % I_sdc_ends = info_sdc.CurrentHistory;
% % plot(t_sdc_ends*1e3, I_sdc_ends, 'ro', 'MarkerFaceColor', 'w', 'MarkerSize', 6, ...
% %     'DisplayName', 'SDC Slab Endpoints');
% 
% grid on;
% xlim([0,timeSim*1e3]);
% % legend('Location', 'best');
% xlabel('Time (ms)'); ylabel('Current (A)');
% title(sprintf('Inrush Current: Dense BDF2 vs Coarse SDC (dt=%.1fms)', dt_slab*1e3));
% 
% % 子图 2: 磁密对比
% subplot(2, 1, 2);
% % plot(t_bdf*1e3, B_bdf, 'k-', 'LineWidth', 1.5, 'DisplayName', 'BDF2');
% hold on;
% plot(t_sdc_full*1e3, B_sdc_full, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SDC (Smooth)');
% grid on;
% xlim([0,timeSim*1e3]);
% % legend('Location', 'best');
% xlabel('Time (ms)'); ylabel('|B| at Origin (T)');
% title('B-Field at Core Center');
% 
% % 误差分析
% % I_bdf_interp = interp1([0; t_bdf], [0; I_bdf], t_sdc_full, 'linear', 'extrap');
% % rel_err_norm = norm(I_sdc_full - I_bdf_interp) / norm(I_bdf_interp);
% % 
% % fprintf('\n--- Accuracy Report ---\n');
% % fprintf('Relative Error (Curve Match): %.2e%%\n', rel_err_norm * 100);
% 
% fprintf('\nDone.\n');


%% --- 辅助局部函数 (请将此函数放在脚本的最末尾) ---
function monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist)
    % 创建或激活图形窗口
    persistent hFig;
    if isempty(hFig) || ~isvalid(hFig)
        hFig = figure('Name', 'Real-time Monitor', 'Position', [100, 100, 800, 600]);
    end
    set(0, 'CurrentFigure', hFig);
    
    % --- 子图 1: 电流波形 ---
    subplot(2, 1, 1);
    plot(t_vec, I_vec, 'r', 'LineWidth', 1);
    grid on;
    ylabel('Current (A)', 'FontSize', 10);
    title(sprintf('Time: %.4f s | Current: %.4f A', t, I), 'FontSize', 11);
    xlim([0, max(t_vec(end), 1e-6)]); % 动态调整轴范围
    
    % --- 子图 2: 磁密波形 ---
    subplot(2, 1, 2);
    if isempty(B_hist)
        % 兼容性处理: 如果尚未生成 B 数据 (第一步)
        B_plot = 0; t_plot = 0;
    else
        B_plot = B_hist; t_plot = t_vec;
    end
    plot(t_plot, B_plot, 'b', 'LineWidth', 1);
    grid on;
    ylabel('|B| Magnitude (T)', 'FontSize', 10);
    xlabel('Time (s)', 'FontSize', 10);
    title(sprintf('Probe |B|: %.4f T', B_curr), 'FontSize', 11);
    xlim([0, max(t_vec(end), 1e-6)]);
    
    drawnow; % 限制刷新率以保证计算速度
end