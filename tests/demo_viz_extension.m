% demo_viz_extension.m
% 演示 magpackage 的可视化与后处理新功能 (Precision Fix v3)
% ---------------------------------------------------------
clear; clc; close all;

% 1. 环境初始化
if exist('../src', 'dir')
    addpath(genpath('../src'));
    fprintf('[Info] Added ../src to path.\n');
elseif exist('magpackage/src', 'dir')
    addpath(genpath('magpackage/src')); 
    fprintf('[Info] Added magpackage/src to path.\n');
else
    warning('Please ensure "src" folder is in the path.');
end

%% --- Step 1: 准备数据 (生成简单网格与解析场) ---
fprintf('\n[Step 1] Preparing Mock Problem...\n');

% 1.1 生成网格 (加密步长以提高可视化质量)
% 使用 h=0.1 保证网格质量
[P, T] = generate_simple_box_mesh(1.0, 1.0, 0.5, 0.1); 
mesh = Mesh();
mesh.P = P; mesh.T = T;
mesh.RegionTags = ones(1, size(T, 2)); 
mesh.NumNodes = size(P, 2);
mesh.NumElements = size(T, 2);
mesh.generateEdges();

fprintf('   Generated Mesh: Nodes=%d, Elems=%d, Edges=%d\n', ...
    mesh.NumNodes, mesh.NumElements, size(mesh.Edges, 2));

% 1.2 初始化 DofHandler
fprintf('   Initializing DofHandler...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

% 1.3 构造解析解 (使用精确代数公式)
% 目标: A = (-0.5*y, 0.5*x, 0) => B = (0, 0, 1)
fprintf('   Projecting analytical field A to edges (Exact Formula)...\n');

numDofs = dofHandler.NumGlobalDofs;
A_sol = zeros(numDofs, 1);
is_set = false(numDofs, 1); 

% 获取 DoF 映射表
key = space.toString();
if ~dofHandler.DofMaps.isKey(key)
    error('Key "%s" not found in DofMaps.', key);
end
dofMap = dofHandler.DofMaps(key); 

% 局部边定义 [StartNode, EndNode]
local_edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

for e = 1:mesh.NumElements
    elem_dofs = dofMap(:, e);
    node_ids = mesh.T(:, e);
    P_nodes = mesh.P(:, node_ids); % [3x4]
    signs = double(mesh.T2E_Sign(:, e));
    
    for i = 1:6
        gid = elem_dofs(i);
        if is_set(gid), continue; end
        
        % 获取边端点
        n1 = local_edge_defs(i, 1);
        n2 = local_edge_defs(i, 2);
        P1 = P_nodes(:, n1);
        P2 = P_nodes(:, n2);
        
        % 坐标差分
        dx = P2(1) - P1(1);
        dy = P2(2) - P1(2);
        % dz = P2(3) - P1(3); % Az=0, 所以 dz 对点积无贡献
        
        % 精确积分公式: Integral A*dl
        % A = (-0.5y, 0.5x, 0)
        % Int = 0.5 * (x1*dy - y1*dx)
        % (这是该线性场沿直线段积分的解析解，无任何数值近似)
        val_local = 0.5 * (P1(1) * dy - P1(2) * dx);
        
        % 转换到全局方向
        A_sol(gid) = val_local * signs(i);
        is_set(gid) = true;
    end
end

fprintf('   Data prepared. A_sol norm = %.4f\n', norm(A_sol));

%% --- Step 2: 初始化可视化模块 ---
fprintf('\n[Step 2] Initializing Visualizer...\n');
assembler = Assembler(mesh, dofHandler); 
post = PostProcessor(assembler);
viz = Visualizer(post);

% --- 自检: 验证 B 场一致性 ---
fprintf('   [Diagnostic] Checking Element B values...\n');
B_elems = post.computeElementB(A_sol);
B_mag = post.computeMagnitude(B_elems);
fprintf('   -> Mean |B|: %.6f (Expected 1.0)\n', mean(B_mag));
fprintf('   -> Min  |B|: %.6f\n', min(B_mag));
fprintf('   -> Max  |B|: %.6f\n', max(B_mag));
if std(B_mag) > 1e-4
    warning('B field has significant variance! Check mesh or projection.');
end

%% --- Step 3: 新功能演示 ---

% === 功能 1: 网格可视化 ===
fprintf('[Demo 1] Visualizing Mesh (Outer Surface)...\n');
figure; subplot(2,2,1);
viz.plotMesh('Alpha', 0.6, 'FaceColor', [0.8, 0.8, 0.8]);
title('Mesh Outer Surface');
drawnow;

% === 功能 2: 矢量场绘制 ===
fprintf('[Demo 2] Visualizing Vector Field...\n');
subplot(2,2,2);
pts_random = rand(50, 3) .* [0.8, 0.8, 0.4] + [0.1, 0.1, 0.05]; 
viz.plotVectorField(A_sol, pts_random);
title('Magnetic Flux Density B (Corrected)');
drawnow;

% === 功能 3: 平滑切片云图 ===
fprintf('[Demo 3] Visualizing Smoothed Slice...\n');
subplot(2,2,3);
viz.plotSliceSmoothed(A_sol, 'z', 0.25, 100); 
title('Smoothed |B| Slice at Z=0.25m');
clim([0.95, 1.05]); % 聚焦颜色范围
drawnow;

% === 功能 4: 沿线数据提取 ===
fprintf('[Demo 4] Probing Data along Line...\n');
subplot(2,2,4);
% 避免正好落在网格线上的采样点 (加 eps)
startPt = [0, 0.5 + 1e-6, 0.25 + 1e-6];
endPt   = [1, 0.5 + 1e-6, 0.25 + 1e-6];
N_samples = 100;

data_B = post.probeLine(A_sol, startPt, endPt, N_samples, 'B');

dist_axis = linspace(0, 1, N_samples);
% 绘制 Bz
plot(dist_axis, data_B(:,3), 'b.-', 'LineWidth', 1.0, 'DisplayName', 'FEM Bz');
hold on;
yline(1.0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
% 设置 Y 轴范围以清晰显示精度
ylim([0.9, 1.1]); 
xlabel('X (m)'); ylabel('Bz (T)');
title('Bz Distribution (Should be flat 1.0)');
legend; grid on;

fprintf('\n[Done] All demonstrations completed successfully.\n');

%% --- 辅助函数 ---
function [P, T] = generate_simple_box_mesh(Lx, Ly, Lz, h)
    % 生成稍微扰动的点阵，避免完美对称导致的数值简并
    x = 0:h:Lx; y = 0:h:Ly; z = 0:h:Lz;
    [X, Y, Z] = meshgrid(x, y, z);
    
    % 添加微小随机扰动 (1% h) 以提高 Delaunay 鲁棒性
    noise = (rand(size(X)) - 0.5) * h * 0.01;
    X = X + noise; Y = Y + noise; Z = Z + noise;
    
    % 边界点强制归位 (保持盒子形状)
    X(X < 0) = 0; X(X > Lx) = Lx;
    Y(Y < 0) = 0; Y(Y > Ly) = Ly;
    Z(Z < 0) = 0; Z(Z > Lz) = Lz;
    
    P = [X(:)'; Y(:)'; Z(:)'];
    
    try
        dt = delaunayTriangulation(P');
        T = dt.ConnectivityList';
        P = dt.Points';
    catch
        error('Delaunay failed.');
    end
end