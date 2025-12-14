% run_ecore_transient.m - E-Core Transformer Transient Benchmark
% 
% 描述:
%   本脚本模拟 E 型变压器在正弦电压源激励下的瞬态响应。
%   这是一个典型的场路耦合(Field-Circuit Coupled)问题。
% 
% 物理模型特征:
%   1. 几何: E 型铁芯 + 跑道型/圆形线圈。
%   2. 材料: 铁芯采用非线性 B-H 曲线 (模拟饱和效应)。
%   3. 物理场: 磁准静态场 (Magnetoquasistatic)。
%      - 设置全局电导率 sigma=0，即忽略涡流效应 (No Eddy)。
%      - 引入 Lagrange 乘子 P 场以施加库伦规范 (div A = 0)，保证唯一解。
%   4. 激励: 电压源驱动 (Voltage Driven)，需联立求解线圈电流 I。

clear; clc; tic;
fprintf('=========================================================\n');
fprintf('   E-Core Transient Analysis (No Eddy, No V-Dofs)        \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
fprintf('[Step 1] Loading Mesh and Geometry...\n');

% 1.1 加载网格
% 网格文件包含四面体单元信息。'm' 表示单位为米。
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 

% 1.2 生成拓扑
% 生成棱边(Edges)信息，这是使用 Nedelec 棱边元的前提。
mesh.generateEdges();

% 1.3 区域 ID 重映射 (Retagging)
% 原始网格中线圈可能被切分为多个部分 (如 ID 5,6,8,9)。
% 为了方便物理定义，将其合并为统一的 ID (NEW_TAG_PRIM)。
tags_primary   = [5, 6, 8, 9];   % 原边线圈各部分 ID
tags_secondary = [3, 4, 7, 10];  % 副边线圈各部分 ID
tag_core       = 2;              % 铁芯 ID
tag_air        = 1;              % 空气/背景 ID

NEW_TAG_PRIM = 50; % 新的原边 ID
NEW_TAG_SEC  = 60; % 新的副边 ID

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料与物理属性配置 ---
fprintf('[Step 2] Defining Materials and Physics...\n');

% 初始化材料库容器
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 2.1 铁芯材料 (非线性 B-H)
% 定义 B(T) 与 H(A/m) 的对应关系，模拟磁饱和特性。
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 2.2 线性材料 (空气/线圈)
% 这里假设空气和铜线的相对磁导率均为 1.0
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% 2.3 电导率设置 (Sigma)
% 本案例不分析涡流损耗，因此将全场电导率设为 0。
% 注意：即使 Sigma=0，瞬态求解器仍需计算 d/dt 项(来自外电路耦合)。
sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags)
    sigmaMap(all_tags(i)) = 0.0;
end

% 2.4 有限元空间定义
% space_A: 使用一阶 Nedelec 棱边元求解磁矢位 A (Curl-Conforming)
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

% space_P: 使用一阶 Lagrange 节点元求解标量位 P (用于规范固定)
% 在非导电区域，需要 P 场来保证库伦规范 div A = 0，防止刚度矩阵奇异。
space_P = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_P);

%% --- 3. 线圈与电路耦合设置 ---
fprintf('[Step 3] Setup Field-Circuit Coupling...\n');

% 3.1 自动计算线圈几何信息
% 辅助工具自动识别线圈的中心、半径、截面积和轴向。
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);

% 3.2 计算线圈方向场
% 生成符合右手螺旋定则的切向单位向量场，用于定义电流流向。
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);
fprintf('   -> Coil Config: R=%.4fm, Clockwise Direction Set.\n', radius);

% 3.3 创建绕组对象
% 定义物理属性：匝数=3000, 直流电阻=10欧姆
NP = 3000;
RP = 10;
winding = Winding('Primary', NEW_TAG_PRIM, NP, RP, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map); % 赋予空间分布的方向场

% 3.4 外电路定义 (Netlist)
% 电路方程: V_source(t) = i(t)*R + L*di/dt + dPsi/dt (感应电动势)
circuit = struct();
circuit.R = RP;  % 外电路串联电阻
circuit.L = 0;   % 外电路串联电感 (不包括线圈自感)
% 定义电压源函数 (50Hz 正弦波, 幅值 500V)
Vamp = 500;
circuit.V_source_func = @(t) Vamp * sin(2 * pi * 50 * t);

%% --- 4. 求解器配置 ---
fprintf('[Step 4] Configuring Solver...\n');

assembler = Assembler(mesh, dofHandler);
solver = TransientCoupledSolver(assembler);

% 4.1 迭代控制参数
solver.Tolerance = 1e-2;    % 绝对收敛容差
solver.RelTolerance = 1e-2; % 相对收敛容差 (瞬态分析通常 1% 足够)

% 4.2 线性求解器配置 (MUMPS)
% 场路耦合矩阵是非对称且不定的 (Saddle Point)，必须设置 Symmetry=0。
solver.LinearSolver.MumpsSymmetry = 0;
solver.LinearSolver.MumpsICNTL.i14 = 300; % 增加内存预分配以防溢出

% 4.3 时间步进设置
dt = 5e-4;         % 时间步长
timeSim = 0.1-dt;  % 模拟时长，观察启动瞬态
timeSteps = repmat(dt, round(timeSim/dt), 1);

% 4.4 边界条件 (Dirichlet)
% 在外边界处强制 A x n = 0 (磁壁边界) 和 P = 0
fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);
fixedDofs_P = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_P);

% 4.5 实时绘图回调
% 定义磁密探测点 (例如选在铁芯柱中心 [0,0,0] 或其他感兴趣的位置)
probePoint = [0, 0, 0];

% 定义绘图回调函数句柄
% 注意: 求解器会自动传入 (t, I, t_vec, I_vec, B_curr, B_hist)
plotFunc = @(t, I, t_vec, I_vec, B_curr, B_hist) monitor_callback(t, I, t_vec, I_vec, B_curr, B_hist);

%% --- 5. 求解 ---
fprintf('[Step 5] Solving Transient System...\n');

% 执行求解
% solve 方法会自动处理时间步进、非线性牛顿迭代及场路耦合更新
[sols, info] = solver.solve(space_A, space_P, matLib, sigmaMap, ...
                            circuit, winding, timeSteps, ...
                            fixedDofs_A, fixedDofs_P, [], plotFunc, probePoint);

%% --- 6. 后处理 ---
fprintf('[Step 6] Post-Processing...\n');

post = PostProcessor(assembler);
viz  = Visualizer(post);

% 6.1 数据提取
num_steps = length(sols);
time_vec = cumsum(timeSteps);
current_vec = info.CurrentHistory; % 获取电流历史数据

% 6.2 绘制电流响应曲线
figure('Name', 'E-Core Transient Analysis', 'Position', [100, 100, 1000, 500]);

subplot(1, 2, 1);
plot([0;time_vec], [0;current_vec], 'b-o', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('Current (A)');
title('Primary Coil Current (Nonlinear Response)');
legend('Transient Current');

% 6.3 绘制磁通密度 B 分布 (3D 云图)
% 选择 t = 0.006s 时刻进行观察
target_t = 0.006;
[~, step_idx] = min(abs(time_vec - target_t));

% 注意: sols 是最后一步的解向量，若需特定步长的全场解，需在 solve 中保存快照。
% 此处假设 sols 存储的是最终时刻解，或者为了演示，我们直接使用由 solve 返回的最后状态。
% (注: 标准 solver.solve 返回的是最后一步的 x_curr。若需中间步结果，需修改 solver 保存逻辑)
% 这里我们使用 solver 结束时的状态 A_sol 进行绘制。
A_sol = sols(1:dofHandler.NumGlobalDofs);

subplot(1, 2, 2);
% 计算单元中心的 B 场
B_elems = post.computeElementB(A_sol, 'Nedelec_P1');
% 计算模值 |B|
B_mag = sqrt(sum(abs(B_elems).^2, 1));
% 平滑映射到节点以便绘图
B_nodes = post.mapElementsToNodes(B_mag);
% 仅提取铁芯区域的数据
B_core = post.combineRegionData(B_nodes, tag_core);

% 在铁芯表面绘制 B 模值云图
viz.plotFieldOnSurface(tag_core, B_core(:), 'FaceAlpha', 1.0, 'EdgeColor', 'none');
title(sprintf('|B| Magnitude at t=%.3fs', time_vec(end)));
view(3); axis equal;

toc;

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
    plot(t_vec, I_vec, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
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