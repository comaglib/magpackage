% tests/test_coil_utils_viz.m
clear; clc; close all;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Coil Geometry Utils & Visualization (Fixed v3)  \n');
fprintf('=========================================================\n');

%% === Case 1: 圆角矩形线圈 (Racetrack Coil) ===
fprintf('\n[Case 1] Generating Racetrack Coil Mesh...\n');

% 定义参数
Params1.L_straight = 0.2; 
Params1.R_mean = 0.05;    
Params1.Width = 0.02;     
Params1.Height = 0.02;    
Params1.Center = [0, 0, 0];

% 生成网格 (修复了空洞过滤逻辑)
mesh1 = generate_racetrack_mesh(Params1);

% 1. 自动识别
fprintf('[Analysis] Auto-detecting geometry...\n');
[center, L, R, axis_idx] = CoilGeometryUtils.autoDetectRacetrack(mesh1, 1);
fprintf('   Detected: Center=[%.3f %.3f %.3f], L=%.4f, R=%.4f, Axis=%d\n', ...
    center, L, R, axis_idx);

% 2. 计算流向
fprintf('[Analysis] Computing current direction field...\n');
dir_map1 = CoilGeometryUtils.computeRacetrackDirection(mesh1, 1, center, L, R, axis_idx);

% 3. 可视化
figure('Name', 'Case 1: Racetrack Coil', 'Color', 'w', 'Position', [100, 100, 800, 600]);
visualize_coil(mesh1, dir_map1, 'Racetrack Coil Current Flow');


%% === Case 2: 圆形线圈 (Circular Coil) ===
fprintf('\n[Case 2] Generating Circular Coil Mesh...\n');

Params2.L_straight = 0.0; 
Params2.R_mean = 0.08;    
Params2.Width = 0.02;
Params2.Height = 0.02;
Params2.Center = [0.2, 0, 0.1]; 

mesh2 = generate_racetrack_mesh(Params2); 

fprintf('[Analysis] Auto-detecting geometry...\n');
[center2, L2, R2, axis2] = CoilGeometryUtils.autoDetectRacetrack(mesh2, 1);
fprintf('   Detected: Center=[%.3f %.3f %.3f], L=%.4f, R=%.4f, Axis=%d\n', ...
    center2, L2, R2, axis2);

fprintf('[Analysis] Computing current direction field...\n');
dir_map2 = CoilGeometryUtils.computeRacetrackDirection(mesh2, 1, center2, L2, R2, axis2);

figure('Name', 'Case 2: Circular Coil', 'Color', 'w', 'Position', [150, 150, 800, 600]);
visualize_coil(mesh2, dir_map2, 'Circular Coil Current Flow');


%% === 辅助函数: 网格生成器 ===
function mesh = generate_racetrack_mesh(p)
    L = p.L_straight;
    R = p.R_mean;
    w = p.Width;
    h = p.Height;
    c = p.Center;
    
    n_radial = 3;
    n_height = 3;
    n_theta = 40; 
    n_straight = 15; 
    
    points = [];
    
    r_range = linspace(-w/2, w/2, n_radial);
    z_range = linspace(-h/2, h/2, n_height);
    [RR, ZZ] = meshgrid(r_range, z_range);
    RR = RR(:); ZZ = ZZ(:);
    
    % 右半圆
    theta_right = linspace(-pi/2, pi/2, n_theta);
    for i = 1:length(theta_right)
        th = theta_right(i);
        for k = 1:length(RR)
            x = L/2 + (R + RR(k)) * cos(th);
            y = 0   + (R + RR(k)) * sin(th);
            z = ZZ(k);
            points = [points; x, y, z]; %#ok<AGROW>
        end
    end
    
    % 左半圆
    theta_left = linspace(pi/2, 3*pi/2, n_theta);
    for i = 1:length(theta_left)
        th = theta_left(i);
        for k = 1:length(RR)
            x = -L/2 + (R + RR(k)) * cos(th);
            y = 0    + (R + RR(k)) * sin(th);
            z = ZZ(k);
            points = [points; x, y, z]; %#ok<AGROW>
        end
    end
    
    % 上直段
    x_straight = linspace(-L/2, L/2, n_straight);
    for i = 1:length(x_straight)
        for k = 1:length(RR)
            x = x_straight(i);
            y = R + RR(k);
            z = ZZ(k);
            points = [points; x, y, z]; %#ok<AGROW>
        end
    end
    
    % 下直段
    for i = 1:length(x_straight)
        for k = 1:length(RR)
            x = x_straight(i);
            y = -R - RR(k); 
            z = ZZ(k);
            points = [points; x, y, z]; %#ok<AGROW>
        end
    end
    
    points = points + c;
    
    % Delaunay
    DT = delaunayTriangulation(points);
    T_all = DT.ConnectivityList;
    P_all = DT.Points;
    
    % 过滤内孔
    centers = (P_all(T_all(:,1),:) + P_all(T_all(:,2),:) + P_all(T_all(:,3),:) + P_all(T_all(:,4),:)) / 4.0;
    rel_centers = centers - c;
    
    R_inner = R - w/2;
    Threshold = R_inner * 0.5; % 保守阈值
    
    valid_mask = false(size(T_all, 1), 1);
    
    for i = 1:size(T_all, 1)
        px = rel_centers(i, 1);
        py = rel_centers(i, 2);
        if abs(px) <= L/2
            dist = abs(py);
        else
            dist = sqrt( (abs(px) - L/2)^2 + py^2 );
        end
        if dist >= Threshold
            valid_mask(i) = true;
        end
    end
    
    T_final = T_all(valid_mask, :);
    
    % 构建 Mesh 对象 (标准: P 为 3xN, T 为 4xM)
    mesh = Mesh();
    mesh.P = P_all'; 
    mesh.T = T_final'; 
    mesh.NumNodes = size(mesh.P, 2);
    mesh.NumElements = size(mesh.T, 2);
    mesh.RegionTags = ones(1, mesh.NumElements); 
    
    mesh.generateEdges();
end

%% === 辅助函数: 可视化 (修复版) ===
function visualize_coil(mesh, dir_map, title_str)
    % 使用 MagPackage 标准数据结构: P [3xN], T [4xM]
    P = mesh.P;
    T = mesh.T;
    
    % 1. 绘制外表面云图
    % triangulation 需要 [M x 4] 和 [N x 3]
    TR = triangulation(T', P');
    [F_bnd, ~] = freeBoundary(TR);
    
    % 着色数据: Z 坐标
    C_data = P(3, :); 
    
    % trisurf 需要 x, y, z 向量
    trisurf(F_bnd, P(1,:)', P(2,:)', P(3,:)', C_data, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    
    % 2. 绘制流向矢量
    % 计算重心 [3 x M]
    % P(:, T(1,:)) 是 [3 x M]
    centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
    
    % 降采样
    num_arrows = 500;
    total_elems = mesh.NumElements;
    step = max(1, floor(total_elems / num_arrows));
    idx = 1:step:total_elems;
    
    % 提取采样点和向量 [3 x K]
    pts = centers(:, idx);
    vecs = dir_map(:, idx);
    
    % quiver3 需要 u, v, w (向量)
    quiver3(pts(1,:), pts(2,:), pts(3,:), ...
            vecs(1,:), vecs(2,:), vecs(3,:), ...
            0.5, 'k', 'LineWidth', 1.2);
        
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(title_str);
    view(3);
    colorbar;
    drawnow;
end