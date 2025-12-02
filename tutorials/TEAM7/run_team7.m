% run_team7_v3.m - TEAM Problem 7 Benchmark
clear; clc;
% 请根据实际情况调整 src 路径，这里假设在 tutorials/TEAM7 目录下运行
addpath(genpath('../../src')); 

fprintf('=========================================================\n');
fprintf('   TEAM Problem 7 Benchmark (v3.0 - User Spec)           \n');
fprintf('=========================================================\n');

% --- 1. 网格加载 ---
meshFile = 'C:\Users\haoha\code\magpackage\tutorials\TEAM7\Team7.mphtxt';

fprintf('[Step 1] Loading Mesh from: %s\n', meshFile);
if exist(meshFile, 'file')
    mesh = Mesh.load(meshFile);
else
    error('Mesh file not found: %s', meshFile);
end

% 显式更新统计信息 (防止 mphtxt 读取后未更新)
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

fprintf('   Nodes: %d, Elements: %d, Edges: %d\n', ...
    mesh.NumNodes, mesh.NumElements, size(mesh.Edges, 2));

% --- 2. 物理参数设置 ---
f = 50; 
omega = 2*pi*f;
mu0 = 4*pi*1e-7;
sigma_Al = 3.526e7; % Aluminum S/m

% 线圈电流密度 (TEAM 7 标准值: 2742 AT,截面 100x100mm?)
% 假设总安匝数 2742, 截面需根据几何确认。这里使用近似值 1e6 A/m^2 用于定性验证。
% 或者如果网格是 Quarter 模型，需注意边界条件。这里假设是 Full 或 Half。
J_coil_mag = 2742 * 1; % 假设值，请根据具体几何调整

% 区域定义
% Region 1: Air
% Region 2: Plate (Al)
% Region 3-6: Coil
maxID = max(mesh.RegionTags);
NuMap = ones(maxID, 1) * (1/mu0); % 全域非磁性

SigmaMap = zeros(maxID, 1);
SigmaMap(2) = sigma_Al; % Region 2 是铝板

SourceMap = zeros(maxID, 3);
% 设置线圈电流方向 (需要根据几何确认，这里假设绕 Z 轴旋转或简单的 Y 方向)
% TEAM 7 线圈通常是跑道型。简单起见，这里设为 Y 方向，实际应根据几何计算切向。
SourceMap(3:6, :) = repmat([0, J_coil_mag, 0], 4, 1); 

% --- 3. 自由度与组装 ---
fprintf('[Step 2] Setup Solver...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);
solver.LinearSolver.Method = 'Auto';
solver.LinearSolver.MumpsICNTL.i14 = 60;

% 边界条件 (假设外边界为 PEC)
is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% --- 4. 求解 ---
fprintf('[Step 3] Solving...\n');
A_sol = solver.solve(space, NuMap, SigmaMap, SourceMap, is_bnd);

% --- 5. 后处理 (数据提取) ---
fprintf('[Step 4] Extracting Results (Lines A1B1 - A4B4)...\n');
post = PostProcessor(assembler);

% 定义测量线坐标 (TEAM 7 标准定义，单位 m)
% A1(0, 72, 34) -> B1(288, 72, 34)  (Bz)
% A2(0, 144, 34) -> B2(288, 144, 34) (Bz)
% A3(0, 72, 19) -> B3(288, 72, 19)  (Jy in Plate)
% A4(0, 144, 19) -> B4(288, 144, 19) (Jy in Plate)

lines = struct();
lines(1).name = 'A1-B1'; lines(1).start = [0, 0.072, 0.034]; lines(1).end = [0.288, 0.072, 0.034]; lines(1).type = 'Bz';
lines(2).name = 'A2-B2'; lines(2).start = [0, 0.144, 0.034]; lines(2).end = [0.288, 0.144, 0.034]; lines(2).type = 'Bz';
lines(3).name = 'A3-B3'; lines(3).start = [0, 0.072, 0.019]; lines(3).end = [0.288, 0.072, 0.019]; lines(3).type = 'Jy';
lines(4).name = 'A4-B4'; lines(4).start = [0, 0.144, 0.019]; lines(4).end = [0.288, 0.144, 0.019]; lines(4).type = 'Jy';

N_pts = 50;

figure('Name', 'TEAM 7 Results Comparison', 'Position', [100, 100, 1000, 800]);

for i = 1:4
    ln = lines(i);
    pts = [linspace(ln.start(1), ln.end(1), N_pts)', ...
           linspace(ln.start(2), ln.end(2), N_pts)', ...
           linspace(ln.start(3), ln.end(3), N_pts)'];
    
    vals_real = zeros(N_pts, 1);
    vals_imag = zeros(N_pts, 1);
    
    for k = 1:N_pts
        pt = pts(k, :);
        
        if strcmp(ln.type, 'Bz')
            % 提取 Bz
            [B_val, ~] = post.probeB(A_sol, pt);
            vals_real(k) = real(B_val(3));
            vals_imag(k) = imag(B_val(3));
        else
            % 提取 Jy (使用新封装的 probeEddyCurrent)
            % 注意: SigmaMap 需要传入
            [J_val, ~] = post.probeEddyCurrent(A_sol, pt, omega, SigmaMap);
            vals_real(k) = real(J_val(2));
            vals_imag(k) = imag(J_val(2));
        end
    end
    
    % 绘图
    subplot(2, 2, i);
    x_axis = pts(:, 1) * 1000; % mm
    
    if strcmp(ln.type, 'Bz')
        % Bz 通常显示为 Gauss (1 T = 10000 G)
        plot(x_axis, vals_real*1e4, 'b-', 'LineWidth', 1.5); hold on;
        plot(x_axis, vals_imag*1e4, 'r--', 'LineWidth', 1.5);
        ylabel('Bz (Gauss)');
    else
        % Jy 显示为 A/mm^2 或 MA/m^2
        plot(x_axis, vals_real/1e6, 'b-', 'LineWidth', 1.5); hold on;
        plot(x_axis, vals_imag/1e6, 'r--', 'LineWidth', 1.5);
        ylabel('Jy (MA/m^2)');
    end
    
    title(sprintf('%s (%s)', ln.name, ln.type));
    xlabel('x (mm)');
    legend('Real', 'Imag'); grid on;
end

fprintf('[Done] TEAM 7 Calculation and Plotting Completed.\n');