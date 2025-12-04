% run_team7.m - TEAM Problem 7 Benchmark
%
% 描述:
%   该脚本用于运行 TEAM Problem 7 基准测试案例。
%   案例包含一个非对称导体板（铝板）和一个励磁线圈。
%   求解器使用 A-V 方法计算磁场和涡流分布，并将结果与基准实验数据进行对比。
%
% 流程:
%   1. 加载网格 -> 2. 自动线圈建模 -> 3. 配置求解器 (材料/边界) -> 4. 求解 -> 5. 后处理与可视化

clear; clc;

fprintf('=========================================================\n');
fprintf('   TEAM Problem 7 Benchmark (v19.0 - Final Viz)          \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
% 加载四面体网格文件，单位转换为毫米
meshFile = 'data/meshData/Team7.mphtxt';
fprintf('[Step 1] Loading Mesh from: %s (Unit: mm)\n', meshFile);
if ~exist(meshFile, 'file'), error('Mesh file not found: %s', meshFile); end

mesh = Mesh.load(meshFile, 'mm'); % 加载并自动缩放单位
mesh.generateEdges();             % 生成棱边数据 (Nedelec 单元必需)

%% --- 2. 自动线圈建模 ---
% 为了方便计算源电流密度 J，暂时将所有线圈区域合并为一个临时 Tag
fprintf('[Step 2] Auto-Detecting Coil Geometry...\n');
coil_tag = 9999;
coil_regions = [3, 4, 5, 6]; % 线圈各部分的原始 Region ID

% 临时修改 Tag 以便统一识别几何特征
original_tags = mesh.RegionTags; 
mesh.RegionTags(ismember(mesh.RegionTags, coil_regions)) = coil_tag;

% 调用线圈分析函数，自动计算电流方向和大小 [3 x Ne]
[J_field] = analyzeCoil(mesh, coil_tag);

% 恢复原始 Tag，以免影响后续材料属性赋值
mesh.RegionTags = original_tags; 

%% --- 3. 求解器配置 ---
fprintf('[Step 3] Setup Solver...\n');
f = 50;              % 频率 (Hz)
omega = 2*pi*f;      % 角频率
sigma_Al = 3.526e7;  % 铝板电导率 (S/m)

% 3.1 配置磁导率 (本例假设全场线性材料 mu_r = 1)
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');
matLinear = MaterialLib.createLinear(1.0); 
for i = 1:max(mesh.RegionTags), matLib(i) = matLinear; end

% 3.2 配置电导率 (仅铝板 Region 2 导电，其余为空气)
SigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:max(mesh.RegionTags)
    if i == 2, SigmaMap(i) = sigma_Al; else, SigmaMap(i) = 0.0; end
end

% 3.3 自由度管理 (A-V 公式)
dofHandler = DofHandler(mesh);
space_A = FunctionSpace('Nedelec', 1); % 磁矢位 A (棱边元)
space_V = FunctionSpace('Lagrange', 1);% 电标量位 V (节点元，用于涡流区)
dofHandler.distributeDofs(space_A); 

% 3.4 组装器与求解器初始化
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);
solver.LinearSolver.Method = 'Auto';     % 自动选择线性求解器
solver.LinearSolver.MumpsICNTL.i14 = 60; % 增加内存预估 (MUMPS)

% 3.5 边界条件 (自然边界条件: 外边界切向 A = 0)
is_bnd_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

%% --- 4. 求解 ---
% 执行频域求解 A-V 方程
fprintf('[Step 4] Solving...\n');
% J_field 作为源项传入，is_bnd_A 为狄利克雷边界
[A_sol, V_sol] = solver.solve(space_A, space_V, matLib, SigmaMap, J_field, is_bnd_A, []);

%% --- 5. 后处理与可视化 ---
fprintf('[Step 5] Post-Processing & Visualization...\n');
post = PostProcessor(assembler); % 后处理核心
viz = Visualizer(post);          % 可视化核心

lines = define_measurement_lines(); % 定义测线 (A1-B1 等)

% 创建综合绘图窗口 (2x2 布局)
figure('Name', 'TEAM 7 Comprehensive Analysis', 'Position', [50, 50, 1200, 900]);

% Subplot 1 & 2: 绘制 Bz 沿测线的分布 (仿真 vs 实验数据)
plotBzLine(post, viz, A_sol, lines);

% Subplot 3: 绘制表面磁通密度 |B| 云图
plotBCloudMap(post, viz, A_sol);

% Subplot 4: 绘制铝板表面涡流密度 |J| 云图
plotJCloudMap(post, viz, A_sol, V_sol, omega, SigmaMap, J_field);

fprintf('[Done] Calculation Completed.\n');


%% ========================================================================
%  本地辅助函数
% =========================================================================

function [J_field] = analyzeCoil(mesh, coil_tag)
    % ANALYZECOIL 自动分析线圈几何并生成源电流密度场
    
    % 1. 识别跑道型线圈的几何参数 (中心、半径、长宽、轴向)
    [center, R_in, R_out, Lx, Ly, axis_idx] = ...
        CoilGeometryUtils.autoDetectRoundedRect(mesh, coil_tag);
    
    % 2. 计算电流密度幅值 (安匝数 / 平均截面积)
    AmpereTurns = 2742; 
    tempAssembler = Assembler(mesh, DofHandler(mesh));
    lossCalc = LossCalculator(tempAssembler);
    all_vols = lossCalc.computeElementVolumes();
    
    vol_total = sum(all_vols(mesh.RegionTags == coil_tag)); % 线圈总体积
    R_avg = (R_in + R_out) / 2;
    len_avg = 2*Lx + 2*Ly + 2*pi*R_avg; % 平均周长
    area_avg = vol_total / len_avg;     % 平均截面积
    J_mag = AmpereTurns / area_avg;
    
    % 3. 生成单位方向向量场
    dir_map = CoilGeometryUtils.computeRoundedRectDirection(...
        mesh, coil_tag, center, Lx, Ly, axis_idx);
    
    % 4. 最终电流密度矢量
    J_field = dir_map * J_mag;
end

function lines = define_measurement_lines()
    % 定义 TEAM 7 标准测线坐标 (Bz 和 Jy)
    lines = struct();
    lines(1).name = 'A1-B1(y=72mm,z=34mm)'; lines(1).type = 'Bz';
    lines(1).start = [0, 0.072, 0.034]; lines(1).end = [0.288, 0.072, 0.034]; 
    lines(2).name = 'A2-B2(y=144mm,z=34mm)'; lines(2).type = 'Bz';
    lines(2).start = [0, 0.144, 0.034]; lines(2).end = [0.288, 0.144, 0.034]; 
    lines(3).name = 'A3-B3(y=72mm,z=19mm)'; lines(3).type = 'Jy'; 
    lines(3).start = [0, 0.072, 0.019]; lines(3).end = [0.288, 0.072, 0.019]; 
    lines(4).name = 'A4-B4(y=72mm,z=0mm)'; lines(4).type = 'Jy';
    lines(4).start = [0, 0.072, 0.000]; lines(4).end = [0.288, 0.072, 0.000]; 
end

function [x, re, im] = get_benchmark_data(idx)
    % 获取 TEAM 7 官方基准实验数据 (硬编码)
    X_full = [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288];
    X_sparse = [0, 18, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288];
    if idx == 1 
        x = X_full;
        re = [-4.9, -17.88, -22.13, -20.19, -15.67, 0.36, 43.64, 78.11, 71.55, 60.44, 53.91, 52.62, 53.81, 56.91, 59.24, 52.78, 27.61];
        im = [-1.16, 2.84, 4.15, 4, 3.07, 2.31, 1.89, 4.97, 12.61, 14.15, 13.04, 12.4, 12.05, 12.27, 12.66, 9.96, 2.36];
    elseif idx == 2
        x = X_full;
        re = [-1.83, -8.5, -13.6, -15.21, -14.48, -5.62, 28.77, 60.34, 61.84, 56.64, 53.4, 52.36, 53.93, 56.82, 59.48, 52.08, 26.56];
        im = [-1.63, -0.6, -0.43, 0.11, 1.26, 3.4, 6.53, 10.25, 11.83, 11.83, 11.01, 10.58, 10.8, 10.54, 10.62, 9.03, 1.79];
    elseif idx == 3
        x = X_sparse; % Jy 数据点较少
        re = [0.249, 0.685, -0.015, -0.103, -0.061, -0.004, 0.051, 0.095, 0.135, 0.104, -0.321, -0.687];
        im = [-0.629, -0.873, -0.593, -0.249, -0.101, -0.001, 0.087, 0.182, 0.322, 0.555, 0.822, 0.855];
    elseif idx == 4
        x = X_sparse;
        re = [0.461, 0.621, 1.573, 0.556, 0.237, 0.097, -0.034, -0.157, -0.305, -0.478, -0.66, -1.217];
        im = [-0.662, -0.644, -1.027, -0.757, -0.364, -0.149, 0.015, 0.154, 0.311, 0.508, 0.747, 1.034];
    else
        x=[]; re=[]; im=[];
    end
end

function [x_out, val_re, val_im] = extract_line_data_robust(post, line_info, A_sol, N_pts)
    % 沿线采样并计算 B 场数值
    pts = [linspace(line_info.start(1), line_info.end(1), N_pts)', ...
           linspace(line_info.start(2), line_info.end(2), N_pts)', ...
           linspace(line_info.start(3), line_info.end(3), N_pts)'];
    vals_re = zeros(N_pts, 1); vals_im = zeros(N_pts, 1);
    
    for k = 1:N_pts
        [B_val, ~] = post.probeB(A_sol, pts(k,:)); 
        val = B_val(3) * 1e4; % Tesla -> Gauss (1T = 10000G)
        vals_re(k) = real(val); 
        vals_im(k) = -imag(val); % 注意相位的符号约定
    end
    x_out = pts(:, 1) * 1000; % m -> mm
    val_re = vals_re; val_im = vals_im;
end

function plotBzLine(post, viz, A_sol, lines)
    % 绘制 Bz 对比曲线 (Subplot 1 & 2)
    for i = 1:2
        ln = lines(i);
        % 获取基准数据
        [bench_x, bench_re, bench_im] = get_benchmark_data(i);
        % 获取仿真数据
        [plot_x, plot_re, plot_im] = extract_line_data_robust(post, ln, A_sol, 100);
        
        simData.x = plot_x; simData.re = plot_re; simData.im = plot_im;
        refData.x = bench_x; refData.re = bench_re; refData.im = bench_im;
        
        meta.title = sprintf('%s (%s)', ln.name, ln.type);
        meta.xlabel = 'x (mm)';
        meta.ylabel = 'Bz (Gauss)';
        
        subplot(2, 2, i);
        viz.plotLineComparison(simData, refData, meta);
    end
end

function plotBCloudMap(post, viz, A_sol)
    % 绘制 B 场云图 (Subplot 3)
    subplot(2, 2, 3);
    fprintf('[Visual] Plotting B-Field Cloud Map...\n');
    
    % 1. 计算单元 B 场
    B_elems = post.computeElementB(A_sol, 'Nedelec_P1'); 
    B_mag_elems = sqrt(sum(abs(B_elems).^2, 1)); 
    
    % 2. 映射到节点 (解决材料界面不连续问题)
    B_node_map = post.mapElementsToNodes(B_mag_elems);
    
    % 3. 合并需要绘制的区域 (导体 + 空气包层)
    target_regions = [2, 3, 4, 5, 6];
    B_mag_nodes_combined = post.combineRegionData(B_node_map, target_regions);
    
    % 4. 调用 Visualizer 绘图
    viz.plotFieldOnSurface(target_regions, B_mag_nodes_combined(:), ...
                           'FaceAlpha', 1.0, 'EdgeColor', 'none'); 
    title('Surface Magnetic Flux Density |B| (Tesla)');
end

function plotJCloudMap(post, viz, A_sol, V_sol, omega, SigmaMap, J_field)
    % 绘制涡流 J 场云图 (Subplot 4)
    subplot(2, 2, 4);
    fprintf('[Visual] Plotting Eddy Current J-Field...\n');
    
    % 1. 计算全场电流密度 J = Js + J_eddy
    J_elems_all = post.computeElementJ(A_sol, V_sol, omega, SigmaMap, J_field);
    
    % 2. 计算模值并平滑到节点
    J_mag_elems_all = sqrt(sum(abs(J_elems_all).^2, 1));
    J_node_map = post.mapElementsToNodes(J_mag_elems_all);
    
    % 3. 仅提取铝板区域 (Tag=2)
    plate_tag = 2;
    J_mag_nodes_plate = post.combineRegionData(J_node_map, plate_tag);
    
    % 4. 绘制铝板表面涡流
    viz.plotFieldOnSurface(plate_tag, J_mag_nodes_plate(:), ...
                           'FaceAlpha', 1.0, 'EdgeColor', 'none');
    title('Plate Surface Eddy Current J (A/m^2)');
end