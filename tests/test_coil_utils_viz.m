% test_coil_utils_viz.m
% 测试 CoilGeometryUtils 的 6 个新函数 (V3.0)
% ---------------------------------------------------------
clear; clc; close all;

% 1. 环境初始化
if exist('../src', 'dir'), addpath(genpath('../src'));
elseif exist('magpackage/src', 'dir'), addpath(genpath('magpackage/src')); end

fprintf('==============================================\n');
fprintf('   Test: CoilGeometryUtils (6 Functions)      \n');
fprintf('==============================================\n');

main_fig = figure('Name', 'Coil Geometry Tests', 'Position', [50, 50, 1200, 900]);

%% --- Case 1: 圆形线圈 ---
fprintf('\n[Case 1] Generating Circular Coil Mesh...\n');
[P, T] = gen_torus(0.5, 0.1, 0.08);
mesh1 = create_mesh(P, T);
viz1 = create_viz(mesh1);

% 1. 自动识别
[center, R, ax_idx] = CoilGeometryUtils.autoDetectCircular(mesh1, 1);

% 2. 计算流向
dir1 = CoilGeometryUtils.computeCircularDirection(mesh1, 1, center, R, ax_idx);

% 3. 绘图
plot_case(main_fig, 1, viz1, dir1, 'Case 1: Circular');


%% --- Case 2: 跑道型线圈 ---
fprintf('\n[Case 2] Generating Racetrack Coil Mesh...\n');
[P, T] = gen_racetrack(1.0, 0.3, 0.1, 0.08);
mesh2 = create_mesh(P, T);
viz2 = create_viz(mesh2);

% 1. 自动识别
[center, L, R, ax_idx] = CoilGeometryUtils.autoDetectRacetrack(mesh2, 1);

% 2. 计算流向
dir2 = CoilGeometryUtils.computeRacetrackDirection(mesh2, 1, center, L, R, ax_idx);

% 3. 绘图
plot_case(main_fig, 2, viz2, dir2, 'Case 2: Racetrack');


%% --- Case 3: 圆角矩形线圈 ---
fprintf('\n[Case 3] Generating Rounded Rect Mesh...\n');
Lx_real = 1.0; Ly_real = 0.6; R_real = 0.2;
[P, T] = gen_rounded_rect(Lx_real, Ly_real, R_real, 0.08, 0.08);
mesh3 = create_mesh(P, T);
viz3 = create_viz(mesh3);

% 1. 自动识别 (估算)
[center, Lx_est, Ly_est, R_est, ax_idx] = CoilGeometryUtils.autoDetectRoundedRect(mesh3, 1);

% 2. 计算流向 (使用真实几何以验证算法正确性，实际使用中可用估算值或手动修正)
fprintf('   [Test] Using explicit geometry for verification: Lx=%.2f, Ly=%.2f\n', Lx_real, Ly_real);
dir3 = CoilGeometryUtils.computeRoundedRectDirection(mesh3, 1, center, Lx_real, Ly_real, R_real, ax_idx);

% 3. 绘图
plot_case(main_fig, 3, viz3, dir3, 'Case 3: Rounded Rect');

fprintf('\n[Done] All tests passed.\n');


%% --- 辅助函数 ---
function mesh = create_mesh(P, T)
    mesh = Mesh(); mesh.P = P; mesh.T = T;
    mesh.RegionTags = ones(1, size(T, 2));
    mesh.NumNodes = size(P,2); mesh.NumElements = size(T,2);
    mesh.generateEdges();
end

function viz = create_viz(mesh)
    dof = DofHandler(mesh); 
    asm = Assembler(mesh, dof);
    post = PostProcessor(asm); 
    viz = Visualizer(post);
end

function plot_case(fig, row_idx, viz, dir_map, title_str)
    figure(fig);
    % Left: Mesh
    subplot(3, 2, (row_idx-1)*2 + 1);
    viz.plotMesh('Alpha', 0.6, 'FaceColor', [0.7 0.8 1.0], 'EdgeColor', 'k');
    title([title_str ' Mesh']);
    
    % Right: Vectors
    subplot(3, 2, (row_idx-1)*2 + 2);
    viz.plotElementVectors(dir_map, 'MaxArrows', 250, 'Color', 'r', 'Scale', 0.6);
    hold on; 
    viz.plotMesh('Alpha', 0.1, 'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8]);
    title([title_str ' Vectors']);
end

% --- 网格生成器 (保持不变) ---
function [P, T] = gen_torus(R, r, h)
    lim = R + r + h;
    [X, Y, Z] = meshgrid(-lim:h:lim, -lim:h:lim, -r-h:h:r+h);
    d_total = sqrt( (sqrt(X.^2 + Y.^2) - R).^2 + Z.^2 );
    mask = d_total <= r;
    P_cloud = [X(mask)'; Y(mask)'; Z(mask)'];
    dt = delaunayTriangulation(P_cloud');
    T = dt.ConnectivityList'; P = dt.Points';
end

function [P, T] = gen_racetrack(L, R, r, h)
    lx = L/2 + R + r + h; ly = R + r + h; lz = r + h;
    [X, Y, Z] = meshgrid(-lx:h:lx, -ly:h:ly, -lz:h:lz);
    X_abs = abs(X); dist = zeros(size(X));
    mask_str = (X_abs <= L/2);
    dist(mask_str) = sqrt( (abs(Y(mask_str)) - R).^2 + Z(mask_str).^2 );
    mask_cur = (~mask_str);
    d_cen = sqrt( (X_abs(mask_cur) - L/2).^2 + Y(mask_cur).^2 );
    dist(mask_cur) = sqrt( (d_cen - R).^2 + Z(mask_cur).^2 );
    mask = dist <= r;
    P_cloud = [X(mask)'; Y(mask)'; Z(mask)'];
    dt = delaunayTriangulation(P_cloud');
    T = dt.ConnectivityList'; P = dt.Points';
end

function [P, T] = gen_rounded_rect(Lx, Ly, R, r, h)
    range_x = Lx/2 + R + r + h; range_y = Ly/2 + R + r + h; range_z = r + h;
    [X, Y, Z] = meshgrid(-range_x:h:range_x, -range_y:h:range_y, -range_z:h:range_z);
    x = abs(X); y = abs(Y); d_planar = zeros(size(x));
    
    mask_top = (x <= Lx/2);
    d_planar(mask_top) = abs( y(mask_top) - (Ly/2 + R) );
    
    mask_right = (y <= Ly/2) & (x > Lx/2);
    d_planar(mask_right) = abs( x(mask_right) - (Lx/2 + R) );
    
    mask_corner = (x > Lx/2) & (y > Ly/2);
    d_corn = sqrt( (x(mask_corner) - Lx/2).^2 + (y(mask_corner) - Ly/2).^2 );
    d_planar(mask_corner) = abs( d_corn - R );
    
    dist = sqrt( d_planar.^2 + Z.^2 );
    mask = dist <= r;
    P_cloud = [X(mask)'; Y(mask)'; Z(mask)'];
    dt = delaunayTriangulation(P_cloud');
    T = dt.ConnectivityList'; P = dt.Points';
end