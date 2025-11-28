function Params = analyze_coil_geometry(Mesh, coil_tags)
% ANALYZE_COIL_GEOMETRY [全象限鲁棒版 + 截面分析]
%
% 功能:
% 1. 提取线圈的几何中心 (Center) 和 平均路径尺寸 (Lx_mean, Ly_mean, R_mean)。
% 2. 自动提取线圈的截面尺寸 (Width, Height)，用于多丝建模。
%
% 输入: Mesh, coil_tags
% 输出: Params (.Center, .Lx_mean, .Ly_mean, .R_mean, .Width, .Height)

    % --- 1. 数据准备 ---
    elem_idx = ismember(Mesh.RegionTags, coil_tags);
    if ~any(elem_idx), error('未找到线圈 Region: %s', num2str(coil_tags)); end
    node_idx = unique(Mesh.T(:, elem_idx));
    P = Mesh.P(:, node_idx);
    
    % --- 2. 确定中心与包围盒 ---
    min_p = min(P, [], 2);
    max_p = max(P, [], 2);
    dims = max_p - min_p;
    
    Params.Center = (min_p + max_p) / 2;
    P_centered = P - Params.Center;
    
    % --- 3. 截面尺寸提取 ---
    % 假设线圈平躺在 XY 平面 (TEAM 7 标准)，Z 方向高度即为厚度
    Params.Height = dims(3);
    
    % --- 4. 确定平均跨度 (Lx, Ly) ---
    tol_slice = max(dims) * 0.05;
    
    % [计算 Lx_mean]: 取 y ≈ 0 的切片
    mask_x_arm = abs(P_centered(2,:)) < tol_slice;
    if ~any(mask_x_arm)
        warning('未检测到 X 方向切片，使用 Box 估算。');
        Params.Lx_mean = dims(1) - Params.Height; 
    else
        x_positions = abs(P_centered(1, mask_x_arm));
        Params.Lx_mean = mean(x_positions) * 2;
    end
    
    % [计算 Ly_mean]: 取 x ≈ 0 的切片
    mask_y_arm = abs(P_centered(1,:)) < tol_slice;
    if ~any(mask_y_arm)
        warning('未检测到 Y 方向切片，使用 Box 估算。');
        Params.Ly_mean = dims(2) - Params.Height;
    else
        y_positions = abs(P_centered(2, mask_y_arm));
        Params.Ly_mean = mean(y_positions) * 2;
    end
    
    % [计算径向宽度 Width]
    % Width = Outer_Dimension - Mean_Dimension
    width_x = dims(1) - Params.Lx_mean;
    width_y = dims(2) - Params.Ly_mean;
    Params.Width = (width_x + width_y) / 2;
    
    % --- 5. 拟合圆角半径 R ---
    X_lim = Params.Lx_mean / 2;
    Y_lim = Params.Ly_mean / 2;
    Max_R = min(X_lim, Y_lim); 
    
    mask_corner = (abs(P_centered(1,:)) > X_lim * 0.5) & ...
                  (abs(P_centered(2,:)) > Y_lim * 0.5);
              
    if sum(mask_corner) < 10
        Params.R_mean = 0;
    else
        P_corner = P_centered(:, mask_corner);
        err_fun = @(r) compute_fit_error_full(r, P_corner, X_lim, Y_lim);
        options = optimset('TolX', 1e-4, 'Display', 'off');
        [best_R, ~] = fminbnd(err_fun, 0, Max_R, options);
        if best_R < 1e-3, best_R = 0; end
        Params.R_mean = best_R;
    end
    
    % 调试输出
    fprintf('  [GeoAnalysis] Box: %.1f x %.1f x %.1f mm\n', dims(1)*1000, dims(2)*1000, dims(3)*1000);
    fprintf('  [GeoAnalysis] Mean Path: %.1f x %.1f mm, R=%.1f mm\n', ...
        Params.Lx_mean*1000, Params.Ly_mean*1000, Params.R_mean*1000);
    fprintf('  [GeoAnalysis] Section: Width=%.1f mm, Height=%.1f mm\n', ...
        Params.Width*1000, Params.Height*1000);
end

function err = compute_fit_error_full(R, P, X_lim, Y_lim)
    x = abs(P(1,:)); y = abs(P(2,:));
    cx = X_lim - R; cy = Y_lim - R;
    dist = zeros(1, length(x));
    
    mask_top = (x <= cx);
    dist(mask_top) = abs(y(mask_top) - Y_lim);
    
    mask_right = (y <= cy) & (~mask_top);
    dist(mask_right) = abs(x(mask_right) - X_lim);
    
    mask_arc = (~mask_top) & (~mask_right);
    if any(mask_arc)
        dx = x(mask_arc) - cx; dy = y(mask_arc) - cy;
        dist(mask_arc) = abs(sqrt(dx.^2 + dy.^2) - R);
    end
    err = sqrt(mean(dist.^2));
end