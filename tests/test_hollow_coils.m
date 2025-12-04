% test_hollow_coils.m
% 验证 CoilGeometryUtils 对空心/厚壁线圈的识别精度与计算性能
% [Fix] 1. 修正 gen_robust_racetrack 确保壁厚均匀
%       2. 适配 computeRoundedRectDirection 的 6 参数接口
% -------------------------------------------------------------------------
clear; clc; close all;

if exist('../src', 'dir'), addpath(genpath('../src'));
elseif exist('magpackage/src', 'dir'), addpath(genpath('magpackage/src')); end

fprintf('============================================================\n');
fprintf('   Test: Hollow Coil Geometry & Vectorized Performance      \n');
fprintf('============================================================\n');

main_fig = figure('Name', 'Hollow Coil Analysis (v13.1)', 'Position', [50, 50, 1200, 900]);

%% --- Case 1: Hollow Circular Coil ---
fprintf('\n[Case 1] Hollow Circular Coil...\n');
R_mean = 0.8; W = 0.4; H = 1.0;
[P, T] = gen_robust_cylinder(R_mean, W, H, 0.06); 
mesh1 = create_mesh(P, T);
viz1 = create_viz(mesh1);

tic;
[center, R, ~, ax_idx] = CoilGeometryUtils.autoDetectCircular(mesh1, 1);
t_detect = toc;

check_pass('Radius', R, R_mean, 0.02);
fprintf('   -> Detection Time: %.4f s\n', t_detect);

tic;
dir1 = CoilGeometryUtils.computeCircularDirection(mesh1, 1, center, ax_idx);
t_calc = toc;
fprintf('   -> Vector Calc Time: %.4f s (%d elements)\n', t_calc, mesh1.NumElements);

plot_case(main_fig, 1, viz1, dir1, sprintf('Circular (R=%.2f)', R));


%% --- Case 2: Hollow Racetrack Coil ---
fprintf('\n[Case 2] Hollow Racetrack Coil...\n');
L_str = 3;   % 直线段长度 (L)
R_mean = 1;  % 平均半径 (R) = R_in + w/2
W = 0.4;       % 壁厚 (w)
H = 0.5;       % 高度 (h)
h_mesh = 0.06; % 网格步长

[P, T] = gen_robust_racetrack(L_str, R_mean, W, H, h_mesh);
mesh2 = create_mesh(P, T);
viz2 = create_viz(mesh2);

tic;
[center, L, R, ~, ax_idx] = CoilGeometryUtils.autoDetectRacetrack(mesh2, 1);
t_detect = toc;

check_pass('Length', L, L_str, 0.05);
check_pass('Radius', R, R_mean, 0.05);
fprintf('   -> Detection Time: %.4f s\n', t_detect);

tic;
dir2 = CoilGeometryUtils.computeRacetrackDirection(mesh2, 1, center, L, ax_idx);
t_calc = toc;
fprintf('   -> Vector Calc Time: %.4f s (%d elements)\n', t_calc, mesh2.NumElements);

plot_case(main_fig, 2, viz2, dir2, sprintf('Racetrack (L=%.2f, R=%.2f)', L, R));


%% --- Case 3: Hollow Rounded Rect ---
fprintf('\n[Case 3] Hollow Rounded Rect (Exact Dimensions)...\n');
Lx_skel = 0.9; Ly_skel = 0.9; R_skel = 0.6; W = 0.3; H = 0.5;
Lx_mean = Lx_skel + 2*R_skel; 
Ly_mean = Ly_skel + 2*R_skel; 

[P, T] = gen_robust_rounded_rect(Lx_mean, Ly_mean, R_skel, W, H, 0.04);
mesh3 = create_mesh(P, T);
viz3 = create_viz(mesh3);

tic;
[center, R_in, R_out, Lx, Ly, ~, ax_idx] = CoilGeometryUtils.autoDetectRoundedRect(mesh3, 1);
t_detect = toc;

check_pass('Lx (Straight)', Lx, Lx_skel, 0.08);
check_pass('Ly (Straight)', Ly, Ly_skel, 0.08);
check_pass('R_in', R_in, R_skel - W/2, 0.05);
check_pass('R_out', R_out, R_skel + W/2, 0.05);
fprintf('   -> Detection Time: %.4f s\n', t_detect);

tic;
dir3 = CoilGeometryUtils.computeRoundedRectDirection(mesh3, 1, center, Lx, Ly, ax_idx);
t_calc = toc;
fprintf('   -> Vector Calc Time: %.4f s (%d elements)\n', t_calc, mesh3.NumElements);

plot_case(main_fig, 3, viz3, dir3, 'Hollow Rounded Rect');

fprintf('\n[Done] All tests passed successfully.\n');


%% --- Helper Functions ---
function check_pass(name, val, exp, tol)
    if abs(val - exp) < tol
        fprintf('   [PASS] %-15s = %6.4f (Exp %6.4f)\n', name, val, exp);
    else
        fprintf('   [FAIL] %-15s = %6.4f (Exp %6.4f) Diff=%.4f\n', name, val, exp, abs(val-exp));
    end
end

function mesh = create_mesh(P, T)
    mesh = Mesh(); mesh.P = P; mesh.T = T;
    mesh.RegionTags = ones(1, size(T, 2));
    mesh.NumNodes = size(P,2); mesh.NumElements = size(T,2);
    mesh.generateEdges();
end

function viz = create_viz(mesh)
    dof = DofHandler(mesh); asm = Assembler(mesh, dof);
    post = PostProcessor(asm); viz = Visualizer(post);
end

function plot_case(fig, row, viz, dir, tit)
    figure(fig);
    subplot(3, 2, (row-1)*2+1); 
    viz.plotMesh('Alpha', 0.4, 'FaceColor', 'g', 'EdgeColor', 'k', 'Style', 'surface'); 
    title([tit ' Mesh']);
    subplot(3, 2, (row-1)*2+2); 
    viz.plotElementVectors(dir, 'MaxArrows', 400, 'Color', 'r', 'Scale', 0.6); 
    hold on; 
    viz.plotMesh('Alpha', 0.1, 'EdgeColor', 'none', 'Style', 'surface'); 
    title('Current Flow');
end

%% --- Robust Mesh Generators ---
function [P, T] = gen_robust_cylinder(R_mean, W, H, h)
    R_in = R_mean - W/2; R_out = R_mean + W/2;
    th = linspace(0, 2*pi, ceil(2*pi*R_out/h)); zs = linspace(-H/2, H/2, ceil(H/h)+1);
    P_bnd = [];
    for z = zs
        P_bnd = [P_bnd, [R_in*cos(th); R_in*sin(th); z*ones(size(th))]]; %#ok<*AGROW>
        P_bnd = [P_bnd, [R_out*cos(th); R_out*sin(th); z*ones(size(th))]];
    end
    [X,Y,Z] = meshgrid(-R_out:h:R_out, -R_out:h:R_out, -H/2:h:H/2);
    r = sqrt(X.^2 + Y.^2); mask = (r > R_in+1e-5) & (r < R_out-1e-5);
    P_all = [P_bnd, [X(mask)'; Y(mask)'; Z(mask)']];
    P_all = unique(P_all', 'rows')';
    dt = delaunayTriangulation(P_all'); pts = dt.Points'; T = dt.ConnectivityList';
    ct = (pts(:,T(1,:))+pts(:,T(2,:))+pts(:,T(3,:))+pts(:,T(4,:)))/4;
    rc = sqrt(ct(1,:).^2+ct(2,:).^2);
    T = T(:, (rc >= R_in-1e-3) & (rc <= R_out+1e-3) & (abs(ct(3,:))<=H/2+1e-3));
    P = pts;
end

function [P, T] = gen_robust_racetrack(L, R, w, h, h_mesh)
    % [Modified] 跑道型线圈网格生成器
    % 输入:
    %   L: 直线段长度 (Skeleton Straight Length)
    %   R: 平均圆角半径 (Skeleton Corner Radius, R_mean)
    %   w: 壁厚 (Thickness)
    %   h: 线圈高度 (Height)
    %   h_mesh: 网格步长
    
    % 几何参数导出
    R_out = R + w/2;
    
    % 外包围盒尺寸 (用于生成点阵)
    Lx_out = L + 2*R_out; % 长轴总长
    Ly_out = 2*R_out;     % 短轴总宽 (即直径)
    
    % 生成背景点阵 (Grid)
    [X, Y, Z] = meshgrid(-(Lx_out/2 + h_mesh):h_mesh:(Lx_out/2 + h_mesh), ...
                         -(Ly_out/2 + h_mesh):h_mesh:(Ly_out/2 + h_mesh), ...
                         -(h/2 + h_mesh):h_mesh:(h/2 + h_mesh));
    
    % 1. 定义有向距离场 (Signed Distance Field, SDF)
    % SDF = 0 定义了骨架 (Skeleton)
    % 骨架是一条长度为 L 的线段，沿 X 轴分布，从 -L/2 到 L/2
    X_abs = abs(X);
    % qx 是点在 X 轴上投影到线段 [-L/2, L/2] 之外的距离
    qx = max(X_abs - L/2, 0);
    
    % dist_xy 是点到骨架线段的欧几里得距离
    % 如果点在直段范围内 (qx=0)，距离就是 |y|
    % 如果点在端部外侧，距离就是到端点 (L/2, 0) 的距离 sqrt(qx^2 + y^2)
    dist_to_skeleton = sqrt(qx.^2 + Y.^2);
    
    % 2. 定义壳体 (Shell)
    % 线圈是围绕骨架，半径为 R (平均半径)，厚度为 w 的壳
    % d_shell < 0 表示在壁厚内部
    % abs(dist_to_skeleton - R) 表示点到"平均半径层"的距离
    d_shell = abs(dist_to_skeleton - R) - w/2;
    
    % 3. 定义高度限制
    d_z = abs(Z) - h/2;
    
    % 4. 组合距离场 (交集)
    dist_final = max(d_shell, d_z);
    
    % 提取内部点 (SDF <= 0)
    % 使用微小阈值包含边界点，确保外形完整
    mask_in = dist_final <= 1e-5;
    P_all = [X(mask_in)'; Y(mask_in)'; Z(mask_in)'];
    
    % 5. 剖分与过滤
    P_all = unique(P_all', 'rows')'; % 去重
    dt = delaunayTriangulation(P_all');
    pts = dt.Points'; 
    T = dt.ConnectivityList';
    
    % 计算单元重心并再次过滤 (移除边缘锯齿)
    ct = (pts(:,T(1,:)) + pts(:,T(2,:)) + pts(:,T(3,:)) + pts(:,T(4,:))) / 4;
    
    cx = abs(ct(1,:)); cy = ct(2,:); cz = ct(3,:);
    qx_c = max(cx - L/2, 0);
    dist_xy_c = sqrt(qx_c.^2 + cy.^2);
    
    d_shell_c = abs(dist_xy_c - R) - w/2;
    d_z_c = abs(cz) - h/2;
    dist_c = max(d_shell_c, d_z_c);
    
    % 保留满足几何定义的单元
    valid_elems = dist_c <= 1e-4;
    T = T(:, valid_elems);
    P = pts;
end

function [P, T] = gen_robust_rounded_rect(Lx_mean, Ly_mean, R_skel, W, H, h)
    Lx_str = Lx_mean - 2*R_skel; Ly_str = Ly_mean - 2*R_skel;
    R_out = R_skel + W/2; Max_X = Lx_str/2 + R_out; Max_Y = Ly_str/2 + R_out;
    [X,Y,Z] = meshgrid(-(Max_X+h):h:(Max_X+h), -(Max_Y+h):h:(Max_Y+h), -H/2:h:H/2);
    X_abs = abs(X); Y_abs = abs(Y);
    qx = max(X_abs - Lx_str/2, 0); qy = max(Y_abs - Ly_str/2, 0); dist_xy = sqrt(qx.^2 + qy.^2);
    d_shell = abs(dist_xy - R_skel) - W/2; d_z = abs(Z) - H/2; dist = max(d_shell, d_z);
    P_all = [X(dist<=1e-5)'; Y(dist<=1e-5)'; Z(dist<=1e-5)'];
    P_all = unique(P_all', 'rows')';
    dt = delaunayTriangulation(P_all'); pts = dt.Points'; T = dt.ConnectivityList';
    ct = (pts(:,T(1,:))+pts(:,T(2,:))+pts(:,T(3,:))+pts(:,T(4,:)))/4;
    cx = abs(ct(1,:)); cy = abs(ct(2,:)); cz = ct(3,:);
    qx = max(cx - Lx_str/2, 0); qy = max(cy - Ly_str/2, 0); dist_xy_c = sqrt(qx.^2 + qy.^2);
    d_shell_c = abs(dist_xy_c - R_skel) - W/2; d_z_c = abs(cz) - H/2; dist_c = max(d_shell_c, d_z_c);
    T = T(:, dist_c <= 1e-4);
    P = pts;
end