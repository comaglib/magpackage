% run_D800KVA_compare.m
% =========================================================================
% D800KVA 变压器瞬态分析对比测试脚本 (开路涌流)
% 对比对象：
%   1. 标准 BDF2 求解器 (TransientBDF2Solver) - 基准 (Reference)
%   2. 高阶 SDC 求解器 (AcISDCSolver)         - 待验证的高阶算法
%
% 物理工况：
%   - 高压侧 (Primary) 并联励磁，电压 49497V, 50Hz, 初始相位 0。
%   - 低压侧 (Secondary) 开路。
%   - 铁芯非线性 B-H 曲线。
% =========================================================================

clear; clc;
fprintf('=========================================================\n');
fprintf('   D800KVA Comparison: BDF2 vs High-Order SDC            \n');
fprintf('=========================================================\n');

%% --- 1. 模型设置 (Mesh & Physics) ---

fprintf('[Setup] Loading Mesh and Materials...\n');

% 1.1 加载网格
meshFile = 'data/D800KVA.mphtxt'; 
if ~exist(meshFile, 'file')
    warning('Mesh file %s not found. Please ensure the file exists.', meshFile);
end
mesh = Mesh.load(meshFile, 'm');    % 假设单位是 m
mesh.generateEdges();

% 1.2 区域 ID 定义
tags_p1 = [5, 6, 10, 13];    % Primary 1
tags_p2 = [14, 15, 19, 22];  % Primary 2
tags_s1 = [7, 8, 11, 12];    % Secondary 1
tags_s2 = [16, 17, 20, 21];  % Secondary 2
tag_core = [2, 3, 4, 9, 18, 23];
tag_air = 1;

% 区域标签重映射
% 将线圈区域合并为逻辑组用于材料赋值和线圈分析
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

% 线性材料 (空气和线圈) - 假设相对磁导率为 1
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air) = matLinear;
for t = [NEW_TAG_PRIM1, NEW_TAG_PRIM2], matLib(t) = matLinear; end
for t = [NEW_TAG_SEC1, NEW_TAG_SEC2], matLib(t) = matLinear; end

% 电导率设置 (忽略涡流，设为 0)
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

% 1.4 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

space_P = FunctionSpace('Lagrange', 1); 
dofHandler.distributeDofs(space_P);

% 1.5 线圈绕组配置 (关键步骤)
% 物理描述：Primary1 与 Primary2 并联，电压相等，电流方向相反。
% 建模策略：采用"串联等效"模拟对称并联。
%   - 包含所有区域
%   - 施加 2 倍电压
%   - 设置 2 倍电阻 (R_series = R1 + R2)
%   - 计算得到的 Current 为单支路电流 I_sim。
%   - 最终的总电流 I_total = 2 * I_sim。

% 参数
N_turns_single = 2933;
R_single = 23.42;
V_amp_single = 49497;

% 自动计算方向场
fprintf('[Setup] Computing coil directions...\n');
[c1, ~, S1, ax1] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM1);
dir_map1 = CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM1, c1, ax1);

[c2, ~, S2, ax2] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM2);
dir_map2 = CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM2, c2, ax2);

% **关键**：Primary2 电流方向相反，将方向场翻转
dir_map2(1,:) = -dir_map2(1,:);
dir_map2(2,:) = -dir_map2(2,:);
dir_map2(3,:) = -dir_map2(3,:);

% 合并方向场
full_dir_map = dir_map1 + dir_map2; % 因为区域ID不重叠，直接相加即可合并

% 构建合并的 Winding 对象
% Turns = N1 + N2 (串联总匝数)
% Area = S1 + S2 (总截面积，Assembler 会用 Turns/Area 计算电流密度 J，保证 J 一致)
total_area = S1 + S2;
total_turns = N_turns_single * 2;
total_resistance = R_single * 2; 

% 将高压线圈区域合并带入求解器
tags_hv_all = 71; 
mesh.RegionTags(ismember(mesh.RegionTags, NEW_TAG_PRIM1)) = tags_hv_all;
mesh.RegionTags(ismember(mesh.RegionTags, NEW_TAG_PRIM2)) = tags_hv_all;

% 构造线圈
winding = Winding('Primary_Parallel_Eq', tags_hv_all, total_turns, total_resistance, total_area, [0,0,0]);
winding.setDirectionField(full_dir_map);

% 1.6 外电路参数
circuit = struct();
circuit.R = total_resistance; 
circuit.L = 0; 
% 电压源: 2 * V_single (串联等效)
circuit.V_source_func = @(t) (2 * V_amp_single) * sin(2 * pi * 50 * t);

% 1.7 边界条件与组装器
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);
assembler = Assembler(mesh, dofHandler);

% 1.8 对比磁密探测点
probePoint = [0, 0, -0.6]; 

% 1.9 绘图回调函数句柄
plotFunc = @(t, I, t_vec, I_vec, B_curr, B_hist) monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist);

% %% --- 2. 运行 SDC (High-Order) ---
% fprintf('\n[Run 2] High-Order SDC Solver (Coarse Step)...\n');
% 
% dt_slab = 5e-3;      % 粗步长
% timeSim = 0.02;      % 仿真时长
% timeSteps_sdc = repmat(dt_slab, 1, round(timeSim/dt_slab));
% 
% solver_sdc = AcISDCSolver(assembler);
% solver_sdc.PolyOrder = 4;       
% solver_sdc.MaxSDCIters = 10;    
% solver_sdc.SDCTolerance = 1e-3; 
% 
% tic;
% [~, info_sdc] = solver_sdc.solve(space_A, space_P, matLib, sigmaMap, ...
%                                  circuit, winding, timeSteps_sdc, ...
%                                  fixedDofs_A, fixedDofs_P, [], [], probePoint);
% time_sdc = toc;
% fprintf('-> SDC Completed in %.2f seconds.\n', time_sdc);

%% --- 3. 运行 BDF2 (Reference) ---
fprintf('\n[Run 1] Standard BDF2 Solver (Fine Step)...\n');

dt_bdf = 5e-4;       % 细步长
timeSim = 0.02;      % 仿真时长
timeSteps_bdf = repmat(dt_bdf, round(timeSim/dt_bdf), 1);

solver_bdf = TransientBDF2Solver(assembler);
solver_bdf.Tolerance = 1e-3; 
solver_bdf.RelTolerance = 1e-3;

tic;
[~, info_bdf] = solver_bdf.solve(space_A, space_P, matLib, sigmaMap, ...
                                 circuit, winding, timeSteps_bdf, ...
                                 fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);
time_bdf = toc;
fprintf('-> BDF2 Completed in %.2f seconds.\n', time_bdf);

%% --- 4. 结果对比绘图 ---
fprintf('\n[Post] Plotting comparison...\n');

% 处理 BDF2 数据
t_bdf = cumsum(timeSteps_bdf);
I_bdf_branch = info_bdf.CurrentHistory; % 这是单支路电流
I_bdf_total = 2 * I_bdf_branch;         % 总励磁电流 = 2 * 支路
B_bdf = info_bdf.ProbeB_History;

% % 处理 SDC 数据
% t_sdc_full = info_sdc.Time_Full;
% I_sdc_branch = info_sdc.Current_Full;
% I_sdc_total = 2 * I_sdc_branch;         % 总励磁电流
% B_sdc_full = info_sdc.ProbeB_Full;

figure('Name', 'D800KVA Analysis: BDF2 vs SDC', 'Position', [200, 200, 1000, 600]);

% 子图 1: 总励磁电流对比
subplot(2, 1, 1);
plot(t_bdf*1e3, I_bdf_total, 'k-', 'LineWidth', 1.5, 'DisplayName', 'BDF2 (Fine)');
hold on;
% plot(t_sdc_full*1e3, I_sdc_total, 'r--', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('SDC (P=%d, dt=%.1fms)', solver_sdc.PolyOrder, dt_slab*1e3));
grid on;
xlim([0, timeSim*1e3]);
% legend('Location', 'best');
xlabel('Time (ms)'); ylabel('Total Excitation Current (A)');
title('D800KVA Inrush Current (Parallel Primary Coils)');

% 子图 2: 磁密对比
subplot(2, 1, 2);
plot(t_bdf*1e3, B_bdf, 'k-', 'LineWidth', 1.5, 'DisplayName', 'BDF2');
hold on;
% plot(t_sdc_full*1e3, B_sdc_full, 'b--', 'LineWidth', 1.5, 'DisplayName', 'SDC');
grid on;
xlim([0, timeSim*1e3]);
% legend('Location', 'best');
xlabel('Time (ms)'); ylabel('|B| at [0,0,-0.6] (T)');
title('Magnetic Flux Density at Reference Point');

fprintf('\nDone.\n');


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
    plot(t_vec, 2*I_vec, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
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
    plot(t_plot, B_plot, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 4);
    grid on;
    ylabel('|B| Magnitude (T)', 'FontSize', 10);
    xlabel('Time (s)', 'FontSize', 10);
    title(sprintf('Probe |B|: %.4f T', B_curr), 'FontSize', 11);
    xlim([0, max(t_vec(end), 1e-6)]);
    
    drawnow limitrate; % 限制刷新率以保证计算速度
end