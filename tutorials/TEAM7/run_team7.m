% run_team7.m - TEAM Problem 7 Benchmark (v10.0 - Field Viz Updated)
%
% 更新日志:
%   - v10.0: [Mesh] 适配单位为 'm' 的新网格文件。
%   - v10.0: [Viz] 新增 B 场云图，与网格图合并在同一窗口。

clear; clc;

fprintf('=========================================================\n');
fprintf('   TEAM Problem 7 Benchmark (v10.0 - Field Viz)          \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
meshFile = 'tutorials/TEAM7/Team7.mphtxt';
fprintf('[Step 1] Loading Mesh from: %s (Unit: m)\n', meshFile);
mesh = Mesh.load(meshFile, 'm');
mesh.generateEdges();

%% --- 2. 自动线圈建模 ---
fprintf('[Step 2] Auto-Detecting Coil Geometry...\n');

coil_regions = [3, 4, 5, 6];
temp_coil_tag = 9999;
original_tags = mesh.RegionTags; 
mesh.RegionTags(ismember(mesh.RegionTags, coil_regions)) = temp_coil_tag;

try
    [center, R_in, R_out, Lx, Ly, axis_idx] = ...
        CoilGeometryUtils.autoDetectRoundedRect(mesh, temp_coil_tag);
    
    AmpereTurns = 2742; 
    tempAssembler = Assembler(mesh, DofHandler(mesh));
    lossCalc = LossCalculator(tempAssembler);
    all_vols = lossCalc.computeElementVolumes();
    vol_total = sum(all_vols(mesh.RegionTags == temp_coil_tag));
    
    R_avg = (R_in + R_out) / 2;
    len_avg = 2*Lx + 2*Ly + 2*pi*R_avg;
    area_avg = vol_total / len_avg;
    J_mag = AmpereTurns / area_avg;
    
    dir_map = CoilGeometryUtils.computeRoundedRectDirection(...
        mesh, temp_coil_tag, center, Lx, Ly, axis_idx);
    J_field = dir_map * J_mag;
    
catch ME
    mesh.RegionTags = original_tags; rethrow(ME);
end
mesh.RegionTags = original_tags; 

%% --- 3. 求解器配置 ---
fprintf('[Step 3] Setup Solver...\n');
f = 50; 
omega = 2*pi*f;
sigma_Al = 3.526e7; 

matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');
matLinear = MaterialLib.createLinear(1.0); 
for i = 1:max(mesh.RegionTags), matLib(i) = matLinear; end

SigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:max(mesh.RegionTags)
    if i == 2, SigmaMap(i) = sigma_Al; else, SigmaMap(i) = 0.0; end
end

dofHandler = DofHandler(mesh);
space_A = FunctionSpace('Nedelec', 1);
space_V = FunctionSpace('Lagrange', 1);
dofHandler.distributeDofs(space_A); 

assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);
solver.LinearSolver.Method = 'Auto';
solver.LinearSolver.MumpsICNTL.i14 = 80;

is_bnd_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

%% --- 4. 求解 ---
fprintf('[Step 4] Solving...\n');
[A_sol, V_sol] = solver.solve(space_A, space_V, matLib, SigmaMap, J_field, is_bnd_A, []);

%% --- 5. 后处理与可视化 ---
fprintf('[Step 5] Post-Processing & Visualization...\n');
post = PostProcessor(assembler);
viz = Visualizer(post); 

lines = define_measurement_lines();

% --- 图1: 数据对比 ---
figure('Name', 'TEAM 7 Results Comparison', 'Position', [50, 100, 1200, 800]);
for i = 1:4
    ln = lines(i);
    [bench_x, bench_re, bench_im] = get_benchmark_data(i);
    [plot_x, plot_re, plot_im] = extract_line_data_robust(post, ln, A_sol, V_sol, J_field, omega, SigmaMap, 50);
    
    simData.x = plot_x; simData.re = plot_re; simData.im = plot_im;
    refData.x = bench_x; refData.re = bench_re; refData.im = bench_im;
    
    meta.title = sprintf('%s (%s)', ln.name, ln.type);
    meta.xlabel = 'x (mm)';
    if strcmp(ln.type, 'Bz'), meta.ylabel = 'Bz (Gauss)'; else, meta.ylabel = 'Jy (MA/m^2)'; end
    
    subplot(2, 2, i);
    viz.plotLineComparison(simData, refData, meta);
    
    if i == 3 || i == 4
        print_error_table_robust(post, ln, bench_x, bench_re, bench_im, A_sol, V_sol, J_field, omega, SigmaMap);
    end
end

% --- 图2: 3D 几何与场量可视化 (合并显示) ---
fprintf('[Visual] Generating 3D Field Plots...\n');
figure('Name', '3D Model & B-Field', 'Position', [100, 100, 1400, 600]);

% 子图 1: 铝板网格 + 探测点 (半透明，检查几何)
subplot(1, 2, 1);
[bench_x, ~, ~] = get_benchmark_data(3);
line_info = lines(3); % A3-B3
y_line = line_info.start(2);
z_line = line_info.start(3);
probe_pts = zeros(length(bench_x), 3);
for k = 1:length(bench_x), probe_pts(k, :) = [bench_x(k)/1000, y_line, z_line]; end

% 调用 Visualizer 绘制表面和点
viz.plotSurfaceWithProbes(0.019, probe_pts, 2, 1e-3);
title('Mesh Inspection: Plate Top Surface & Probes');

% 子图 2: 线圈与铝板的表面磁密 B (不透明，物理场)
subplot(1, 2, 2);

% 2.1 计算单元 B 场
B_elems = post.computeElementB(A_sol);
B_mag_elems = sqrt(sum(abs(B_elems).^2, 1)); % 取模值

% 2.2 映射到节点 (用于平滑云图)
B_mag_nodes = post.mapElementsToNodes(B_mag_elems);

% 2.3 使用 Visualizer 绘制云图
% 绘制铝板 (2) 和 线圈 (3,4,5,6)
viz.plotFieldOnSurface([2, 3, 4, 5, 6], B_mag_nodes(:), ...
                       'FaceAlpha', 1.0, ... % 不透明
                       'EdgeColor', 'none'); 
title('Surface Magnetic Flux Density |B| (Tesla)');

fprintf('[Done] Calculation Completed.\n');


%% ========================================================================
%  本地数据提取与物理计算函数
% =========================================================================

function lines = define_measurement_lines()
    lines = struct();
    lines(1).name = 'A1-B1'; lines(1).type = 'Bz';
    lines(1).start = [0, 0.072, 0.034]; lines(1).end = [0.288, 0.072, 0.034]; 
    lines(2).name = 'A2-B2'; lines(2).type = 'Bz';
    lines(2).start = [0, 0.144, 0.034]; lines(2).end = [0.288, 0.144, 0.034]; 
    lines(3).name = 'A3-B3'; lines(3).type = 'Jy';
    lines(3).start = [0, 0.072, 0.019]; lines(3).end = [0.288, 0.072, 0.019]; 
    lines(4).name = 'A4-B4'; lines(4).type = 'Jy';
    lines(4).start = [0, 0.072, 0.000]; lines(4).end = [0.288, 0.072, 0.000]; 
end

function [re, im] = probe_wrapper(post, pt, type, A_sol, V_sol, omega, SigmaMap, J_field)
    [val, elem_idx] = try_probe(post, pt, type, A_sol, V_sol, omega, SigmaMap, J_field);
    need_retry = isnan(elem_idx);
    if ~need_retry && strcmp(type, 'Jy')
        tag = post.Mesh.RegionTags(elem_idx);
        if tag ~= 2, need_retry = true; end
    end
    if need_retry
        d1 = 1e-4; d2 = 5e-4; 
        offsets = [0,0,-d1; 0,0,d1; -d1,0,0; d1,0,0; 0,-d1,0; 0,d1,0;
                   0,0,-d2; 0,0,d2; -d2,0,0; d2,0,0; 0,-d2,0; 0,d2,0;
                   -d1,-d1,-d1; d1,d1,d1];
        for k = 1:size(offsets, 1)
            pt_try = pt + offsets(k, :);
            [val_try, idx_try] = try_probe(post, pt_try, type, A_sol, V_sol, omega, SigmaMap, J_field);
            success = ~isnan(idx_try);
            if success && strcmp(type, 'Jy')
                if post.Mesh.RegionTags(idx_try) ~= 2, success = false; end
            end
            if success, val = val_try; break; end
        end
    end
    re = real(val); im = -imag(val); 
end

function [val, elem_idx] = try_probe(post, pt, type, A_sol, V_sol, omega, SigmaMap, J_field)
    if strcmp(type, 'Bz')
        [B_val, elem_idx] = post.probeB(A_sol, pt); val = B_val(3) * 1e4; 
    else
        [J_tot, elem_idx] = post.probeTotalCurrent(pt, A_sol, V_sol, omega, SigmaMap, J_field); val = J_tot(2) / 1e6; 
    end
end

function [x_out, val_re, val_im] = extract_line_data_robust(post, line_info, A_sol, V_sol, J_field, omega, SigmaMap, N_pts)
    pts = [linspace(line_info.start(1), line_info.end(1), N_pts)', ...
           linspace(line_info.start(2), line_info.end(2), N_pts)', ...
           linspace(line_info.start(3), line_info.end(3), N_pts)'];
    vals_re = zeros(N_pts, 1); vals_im = zeros(N_pts, 1);
    for k = 1:N_pts
        [re, im] = probe_wrapper(post, pts(k,:), line_info.type, A_sol, V_sol, omega, SigmaMap, J_field);
        vals_re(k) = re; vals_im(k) = im;
    end
    x_out = pts(:, 1) * 1000; val_re = vals_re; val_im = vals_im;
end

function print_error_table_robust(post, line_info, bench_x, bench_re, bench_im, A_sol, V_sol, J_field, omega, SigmaMap)
    fprintf('\n----------------------------------------------------------------------\n');
    fprintf(' Error Report: %s (%s)\n', line_info.name, line_info.type);
    fprintf('----------------------------------------------------------------------\n');
    fprintf('%-8s | %-20s | %-20s | %-20s\n', 'X(mm)', 'Real (Sim/Ref)', 'Imag (Sim/Ref)', 'AbsDiff (Re/Im)');
    fprintf('----------------------------------------------------------------------\n');
    y = line_info.start(2); z = line_info.start(3);
    for k = 1:length(bench_x)
        pt = [bench_x(k)/1000, y, z];
        [sim_re, sim_im] = probe_wrapper(post, pt, line_info.type, A_sol, V_sol, omega, SigmaMap, J_field);
        fprintf('%-8.1f | %6.3f / %6.3f        | %6.3f / %6.3f        | %.1e / %.1e\n', ...
            bench_x(k), sim_re, bench_re(k), sim_im, bench_im(k), abs(sim_re-bench_re(k)), abs(sim_im-bench_im(k)));
    end
    fprintf('----------------------------------------------------------------------\n');
end

function [x, re, im] = get_benchmark_data(idx)
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
        x = X_sparse;
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