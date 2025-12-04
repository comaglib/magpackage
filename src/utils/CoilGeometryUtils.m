classdef CoilGeometryUtils
    % COILGEOMETRYUTILS 线圈几何辅助工具 (v13.1 - Vectorized & Commented)
    % 
    % 功能:
    %   自动识别有限元网格中线圈区域的几何特征（中心、半径、长宽），
    %   并基于几何特征高效计算电流密度方向矢量。
    %
    % 核心算法:
    %   1. PCA + MABR: 主成分分析结合最小面积包围盒，自动确定任意空间姿态线圈的局部坐标系 (u,v,n)。
    %   2. 平面切片法 (Plane Intersection): 通过构造切平面与网格求交，精确计算截面面积和重心，
    %      基于物理守恒原理推导骨架尺寸和壁厚，对空心/厚壁线圈极度鲁棒。
    %   3. 向量化计算: 电流方向计算完全向量化，支持百万级单元的毫秒级处理。
    
    methods (Static)
        
        %  1. 圆形线圈 (Circular Coil)
        function [center, R, area, normal_axis_idx] = autoDetectCircular(mesh, regionID)
            % AUTODETECTCIRCULAR 自动识别圆形线圈
            % 输入: mesh对象, 区域ID
            % 输出: 圆心, 半径, 法向轴索引(1=X,2=Y,3=Z)
            
            fprintf('   [CoilUtils] 正在识别圆形线圈 (Region %d)...\n', regionID);
            
            % 借用通用矩形检测算法获取骨架尺寸
            % 对于圆形，Lx_mean 和 Ly_mean 理论上相等且等于直径
            [center, H, Lx_mean, Ly_mean, ~, normal_axis_idx, WallWidths] = ...
                CoilGeometryUtils.detectGeneralRect(mesh, regionID);
                
            Diameter = (Lx_mean + Ly_mean) / 2;
            R = Diameter / 2;
            area = mean(WallWidths) * H;
            
            fprintf('   -> 识别结果: Center=[%.3f %.3f %.3f], R=%.4f\n', ...
                center(1), center(2), center(3), R);
        end
        
        function dir_map = computeCircularDirection(mesh, regionID, center, axis_idx)
            % COMPUTECIRCULARDIRECTION 计算圆形线圈电流方向 (向量化版)
            % 算法: 切向矢量 = (-v, u) / norm
            
            % 1. 获取所有单元重心及所属区域掩码
            centers = CoilGeometryUtils.getElementCenters(mesh); % [3 x Ne]
            mask = (mesh.RegionTags == regionID);
            
            % 2. 预分配结果矩阵
            dir_map = zeros(3, mesh.NumElements);
            if ~any(mask), return; end
            
            % 3. 提取线圈单元的坐标并去中心化
            % pts_local: [3 x N_coil]
            pts_local = centers(:, mask) - center(:);
            
            % 4. 获取局部平面坐标轴索引 (u, v)
            [u_idx, v_idx] = CoilGeometryUtils.getUVIndices(axis_idx);
            u = pts_local(u_idx, :);
            v = pts_local(v_idx, :);
            
            % 5. 向量化计算切向
            mag = sqrt(u.^2 + v.^2);
            valid = mag > 1e-12; % 避免除以零
            
            % 构造切向矢量 (在局部平面内旋转 90 度)
            % Tangent_U = -V / R
            % Tangent_V =  U / R
            
            % 仅处理有效点
            u_valid = u(valid); 
            v_valid = v(valid); 
            inv_mag = 1 ./ mag(valid);
            
            dir_vecs = zeros(3, sum(valid));
            dir_vecs(u_idx, :) = -v_valid .* inv_mag;
            dir_vecs(v_idx, :) =  u_valid .* inv_mag;
            
            % 填回结果矩阵 (利用部分索引)
            % 注意: mask 是全域掩码，valid 是 mask 子集内的掩码
            full_indices = find(mask);
            dir_map(:, full_indices(valid)) = dir_vecs;
        end
        
        %  2. 跑道型线圈 (Racetrack Coil)
        function [center, L_straight, R_curve, area, normal_axis_idx] = autoDetectRacetrack(mesh, regionID)
            % AUTODETECTRACETRACK 自动识别跑道型线圈
            % 跑道型由矩形直段和两端的半圆组成。
            
            fprintf('   [CoilUtils] 正在识别跑道型线圈 (Region %d)...\n', regionID);
            
            % 1. 通用检测获取骨架尺寸
            [center, H, Lx_mean, Ly_mean, R_mean, normal_axis_idx, WallWidths] = ...
                CoilGeometryUtils.detectGeneralRect(mesh, regionID);
            
            % 2. 几何推导
            % 假设 U 轴为长轴 (Lx >= Ly)。如果检测结果反了，逻辑会自动适应。
            if Lx_mean > Ly_mean
                L_straight = Lx_mean - 2 * R_mean;
                R_curve = R_mean;
            else
                L_straight = Ly_mean - 2 * R_mean;
                R_curve = R_mean;
            end
            
            if L_straight < 0, L_straight = 0; end
            
            area = mean(WallWidths) * H;
            
            fprintf('   -> 识别结果: Center=[%.3f %.3f %.3f], L=%.4f, R=%.4f\n', ...
                center(1), center(2), center(3), L_straight, R_curve);
        end
        
        function dir_map = computeRacetrackDirection(mesh, regionID, center, straight_len, axis_idx)
            % COMPUTERACETRACKDIRECTION 计算跑道型电流方向 (向量化版)
            % 算法: 分区处理。直段区沿直线，弯曲区沿圆弧。
            
            centers = CoilGeometryUtils.getElementCenters(mesh); 
            mask = (mesh.RegionTags == regionID);
            dir_map = zeros(3, mesh.NumElements);
            if ~any(mask), return; end
            
            pts_local = centers(:, mask) - center(:);
            [u_idx, v_idx] = CoilGeometryUtils.getUVIndices(axis_idx);
            u = pts_local(u_idx, :);
            v = pts_local(v_idx, :);
            
            % --- 分区掩码 ---
            % 直段区域: |u| <= L/2
            mask_straight = abs(u) <= (straight_len / 2);
            mask_curve = ~mask_straight;
            
            % --- 1. 直段区域计算 ---
            if any(mask_straight)
                % 切向平行于 V 轴，方向由 V 的符号决定 (逆时针假设: -sign(v))
                % 注意: 这里假设长轴是 U。
                v_str = v(mask_straight);
                vecs_str = zeros(3, sum(mask_straight));
                vecs_str(u_idx, :) = -sign(v_str); 
                vecs_str(v_idx, :) = 0;
                
                full_idx = find(mask);
                dir_map(:, full_idx(mask_straight)) = vecs_str;
            end
            
            % --- 2. 弯曲区域计算 ---
            if any(mask_curve)
                u_cur = u(mask_curve);
                v_cur = v(mask_curve);
                
                % 将坐标平移到对应的端部圆心 (±L/2, 0)
                u_shift = u_cur - sign(u_cur) * (straight_len / 2);
                
                mag = sqrt(u_shift.^2 + v_cur.^2);
                valid = mag > 1e-12;
                
                % 圆弧切向
                u_valid = u_shift(valid); 
                v_valid = v_cur(valid);
                inv_mag = 1 ./ mag(valid);
                
                vecs_cur = zeros(3, sum(valid));
                vecs_cur(u_idx, :) = -v_valid .* inv_mag;
                vecs_cur(v_idx, :) =  u_valid .* inv_mag;
                
                full_idx = find(mask);
                sub_idx = find(mask_curve);
                dir_map(:, full_idx(sub_idx(valid))) = vecs_cur;
            end
        end
        
        %  3. 圆角矩形线圈 (Rounded Rectangle Coil)
        function [center, R_in, R_out, Lx, Ly, area, normal_axis_idx] = autoDetectRoundedRect(mesh, regionID)
            % AUTODETECTROUNDEDRECT 自动识别圆角矩形
            % 输出:
            %   Lx, Ly: 骨架直段长度
            %   R_in, R_out: 内外圆角半径
            %   area: 平均截面积
            
            fprintf('   [CoilUtils] 正在识别圆角矩形 (Region %d)...\n', regionID);
            
            % 1. 核心检测 (返回骨架尺寸和双向壁厚)
            [center, H, Lx_mean, Ly_mean, R_mean, normal_axis_idx, WallWidths] = ...
                CoilGeometryUtils.detectGeneralRect(mesh, regionID);
            
            % 2. 导出参数
            Lx = Lx_mean - 2 * R_mean;
            Ly = Ly_mean - 2 * R_mean;
            
            avg_wall = mean(WallWidths);
            R_out = R_mean + avg_wall / 2;
            R_in  = R_mean - avg_wall / 2;
            area = avg_wall * H;
            
            if Lx < 0, Lx = 0; end
            if Ly < 0, Ly = 0; end
            if R_in < 0, R_in = 0; end
            
            fprintf('   -> 识别结果: Center=[%.3f %.3f %.3f]\n', center);
            fprintf('      R_in=%.4f, R_out=%.4f\n', R_in, R_out);
            fprintf('      Lx=%.4f, Ly=%.4f (Wall_U=%.4f, Wall_V=%.4f)\n', ...
                Lx, Ly, WallWidths(1), WallWidths(2));
        end
        
        function dir_map = computeRoundedRectDirection(mesh, regionID, center, Lx, Ly, axis_idx)
            % COMPUTEROUNDEDRECTDIRECTION 计算圆角矩形电流方向 (向量化版)
            
            centers = CoilGeometryUtils.getElementCenters(mesh); 
            mask = (mesh.RegionTags == regionID);
            dir_map = zeros(3, mesh.NumElements);
            if ~any(mask), return; end
            
            pts_local = centers(:, mask) - center(:);
            [u_idx, v_idx] = CoilGeometryUtils.getUVIndices(axis_idx);
            u = pts_local(u_idx, :);
            v = pts_local(v_idx, :);
            
            au = abs(u); av = abs(v);
            
            % --- 分区掩码 ---
            % 1. 上下直段: |u| <= Lx/2
            mask_tb = au <= (Lx / 2);
            
            % 2. 左右直段: |v| <= Ly/2 (且不是上下区)
            mask_lr = (av <= (Ly / 2)) & (~mask_tb);
            
            % 3. 四个角区: 其余部分
            mask_corner = (~mask_tb) & (~mask_lr);
            
            full_idx = find(mask);
            
            % --- 计算: 上下直段 ---
            if any(mask_tb)
                v_sub = v(mask_tb);
                vecs = zeros(3, length(v_sub));
                vecs(u_idx, :) = -sign(v_sub);
                vecs(v_idx, :) = 0;
                dir_map(:, full_idx(mask_tb)) = vecs;
            end
            
            % --- 计算: 左右直段 ---
            if any(mask_lr)
                u_sub = u(mask_lr);
                vecs = zeros(3, length(u_sub));
                vecs(u_idx, :) = 0;
                vecs(v_idx, :) = sign(u_sub);
                dir_map(:, full_idx(mask_lr)) = vecs;
            end
            
            % --- 计算: 角区 ---
            if any(mask_corner)
                u_c = u(mask_corner);
                v_c = v(mask_corner);
                
                % 平移到对应象限的圆心
                C_u = sign(u_c) * (Lx / 2);
                C_v = sign(v_c) * (Ly / 2);
                
                du = u_c - C_u;
                dv = v_c - C_v;
                
                mag = sqrt(du.^2 + dv.^2);
                valid = mag > 1e-12;
                
                du = du(valid); dv = dv(valid); inv_mag = 1./mag(valid);
                
                vecs = zeros(3, sum(valid));
                vecs(u_idx, :) = -dv .* inv_mag;
                vecs(v_idx, :) =  du .* inv_mag;
                
                sub_idx = find(mask_corner);
                dir_map(:, full_idx(sub_idx(valid))) = vecs;
            end
        end
        
        %  核心算法: 通用矩形检测 (General Rect Detection)
        function [center, H, Lx_mean, Ly_mean, R_mean, normal_axis_idx, WallWidths] = detectGeneralRect(mesh, regionID)
            % DETECTGENERALRECT 通过平面求交法精确测量线圈几何参数
            %
            % 步骤:
            % 1. 坐标系对齐: PCA 找法向 -> MABR (最小包围盒) 找平面主轴 -> 强制 U 为长轴。
            % 2. 高度测量: 法向投影范围。
            % 3. 切片分析: 
            %    - 构造垂直于 U 轴和 V 轴的切平面。
            %    - 计算切面与四面体网格的交集面积和重心。
            %    - 利用面积推算壁厚 (WallWidth = Area / 2H)。
            %    - 利用重心推算骨架跨度 (Lx_mean = Dist(Center_L, Center_R))。
            % 4. 圆角拟合: 提取角部点云，拟合最佳骨架半径 R。
            
            % --- 1. 获取点云并对齐 ---
            [pts, ~] = CoilGeometryUtils.getPts(mesh, regionID);
            center = mean(pts, 2);
            pts_c = pts - center;
            
            % PCA 确定法向
            C = cov(pts_c'); [V, D] = eig(C); [~, idx] = sort(diag(D));
            axis_n = V(:, idx(1)); % 法向
            [~, normal_axis_idx] = max(abs(axis_n));
            
            % MABR 确定平面主轴 (消除 PCA 在正方形附近的旋转不确定性)
            base_v = cross(axis_n, [1;0;0]);
            if norm(base_v) < 1e-3, base_v = cross(axis_n, [0;1;0]); end
            base_v = base_v / norm(base_v);
            base_u = cross(axis_n, base_v);
            
            u0 = base_u' * pts_c; v0 = base_v' * pts_c;
            
            % 寻找最小包围盒角度
            cost_func = @(theta) CoilGeometryUtils.bboxArea(theta, u0, v0);
            options = optimset('TolX', 1e-4, 'Display', 'off');
            best_theta = fminbnd(cost_func, 0, pi/2, options);
            
            c = cos(best_theta); s = sin(best_theta);
            axis_u_temp = c * base_u + s * base_v;
            axis_v_temp = -s * base_u + c * base_v;
            
            % 强制 U 为长轴
            u_temp = abs(axis_u_temp' * pts_c);
            v_temp = abs(axis_v_temp' * pts_c);
            
            if (max(v_temp)-min(v_temp)) > (max(u_temp)-min(u_temp))
                axis_u = axis_v_temp; axis_v = axis_u_temp;
            else
                axis_u = axis_u_temp; axis_v = axis_v_temp;
            end
            
            % --- 2. 测量高度 ---
            w = axis_n' * pts_c;
            H = max(w) - min(w);
            if H < 1e-9, H = 1.0; end % 2D 兼容
            
            % --- 3. 切片分析 (Lx_mean, WallWidth_u) ---
            % 切平面法向 V (截取 U 方向的臂)
            % 测量轴 U
            [Area_U, Centroid_U_Pos, Centroid_U_Neg] = ...
                CoilGeometryUtils.computeSliceProperties(mesh, regionID, center, axis_v, axis_u);
            
            WallWidth_u = Area_U / (2 * H);
            
            if isempty(Centroid_U_Pos) || isempty(Centroid_U_Neg)
                Lx_mean = max(abs(axis_u' * pts_c)) * 2 - WallWidth_u; % 回退
            else
                Lx_mean = abs(Centroid_U_Pos - Centroid_U_Neg);
            end
            
            % --- 4. 切片分析 (Ly_mean, WallWidth_v) ---
            % 切平面法向 U (截取 V 方向的臂)
            [Area_V, Centroid_V_Pos, Centroid_V_Neg] = ...
                CoilGeometryUtils.computeSliceProperties(mesh, regionID, center, axis_u, axis_v);
            
            WallWidth_v = Area_V / (2 * H);
            
            if isempty(Centroid_V_Pos) || isempty(Centroid_V_Neg)
                Ly_mean = max(abs(axis_v' * pts_c)) * 2 - WallWidth_v;
            else
                Ly_mean = abs(Centroid_V_Pos - Centroid_V_Neg);
            end
            
            WallWidths = [WallWidth_u, WallWidth_v];
            
            % --- 5. 拟合圆角 R ---
            u = abs(axis_u' * pts_c); v = abs(axis_v' * pts_c);
            X_lim = Lx_mean / 2;
            Y_lim = Ly_mean / 2;
            Max_R = min(X_lim, Y_lim);
            
            mask_corner = (u > X_lim * 0.5) & (v > Y_lim * 0.5);
            if sum(mask_corner) < 10
                R_mean = 0;
            else
                u_c = u(mask_corner); v_c = v(mask_corner);
                err_fun = @(r) CoilGeometryUtils.skeletonFitError(r, u_c, v_c, X_lim, Y_lim);
                options = optimset('TolX', 1e-4, 'Display', 'off');
                R_mean = fminbnd(err_fun, 0, Max_R, options);
                if R_mean < 1e-3, R_mean = 0; end
            end
        end
        
        function [TotalArea, Centroid_Pos, Centroid_Neg] = computeSliceProperties(mesh, regionID, center, plane_n, measure_axis)
            % COMPUTECSLICEPROPERTIES 计算截面属性 (面积和重心)
            % 算法: 找到所有与平面相交的四面体，计算截面多边形(三角形或四边形)的面积和重心。
            % 优化: 预计算距离场，利用向量化筛选。
            
            TotalArea = 0;
            w_pos = 0; a_pos = 0;
            w_neg = 0; a_neg = 0;
            
            mask = (mesh.RegionTags == regionID);
            T_sub = mesh.T(:, mask);
            P = mesh.P;
            
            % 计算距离场 D = (P - Center) . n
            D_all = plane_n' * (P - center);
            D_elem = D_all(T_sub); % [4 x N_sub]
            
            % 筛选跨越平面的单元 (min < 0 < max)
            is_cross = (min(D_elem) < -1e-9) & (max(D_elem) > 1e-9);
            T_cross = T_sub(:, is_cross);
            
            if isempty(T_cross)
                Centroid_Pos=[]; Centroid_Neg=[]; return;
            end
            
            num_cross = size(T_cross, 2);
            edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            
            % 循环处理截面 (虽然有循环，但仅处理截面附近的少量单元，速度可接受)
            % 若需进一步加速，可完全向量化所有交点计算，但这会让代码复杂度剧增。
            for i = 1:num_cross
                node_ids = T_cross(:, i);
                d_loc = D_all(node_ids);
                p_loc = P(:, node_ids);
                
                % 收集交点
                pts = zeros(3, 4); % 最多4个交点
                cnt = 0;
                for k = 1:6
                    idx1 = edges(k,1); idx2 = edges(k,2);
                    d1 = d_loc(idx1); d2 = d_loc(idx2);
                    if sign(d1) ~= sign(d2)
                        cnt = cnt + 1;
                        t = -d1 / (d2 - d1);
                        pts(:, cnt) = p_loc(:,idx1) + t*(p_loc(:,idx2)-p_loc(:,idx1));
                    end
                end
                
                if cnt < 3, continue; end
                active_pts = pts(:, 1:cnt);
                
                % 简单去重 (防止非常接近的点)
                % 由于性能敏感，仅当 cnt=4 时检查
                if cnt == 4
                    if norm(active_pts(:,1)-active_pts(:,4)) < 1e-6
                        cnt=3; active_pts=active_pts(:,1:3);
                    end
                end
                
                % 计算面积 (基于重心三角剖分)
                c_poly = sum(active_pts, 2) / cnt;
                
                % 顶点排序 (投影到局部平面)
                vecs = active_pts - c_poly;
                v_base1 = measure_axis;
                v_base2 = cross(plane_n, measure_axis);
                ang = atan2(v_base2' * vecs, v_base1' * vecs);
                [~, sort_idx] = sort(ang);
                sorted_pts = active_pts(:, sort_idx);
                
                area_poly = 0;
                for k = 1:cnt
                    p1 = sorted_pts(:, k);
                    p2 = sorted_pts(:, mod(k, cnt)+1);
                    area_poly = area_poly + 0.5 * norm(cross(p1-c_poly, p2-c_poly));
                end
                
                TotalArea = TotalArea + area_poly;
                
                % 投影重心到测量轴
                proj = measure_axis' * (c_poly - center);
                if proj > 0
                    w_pos = w_pos + proj * area_poly;
                    a_pos = a_pos + area_poly;
                else
                    w_neg = w_neg + proj * area_poly;
                    a_neg = a_neg + area_poly;
                end
            end
            
            Centroid_Pos = []; Centroid_Neg = [];
            if a_pos > 0, Centroid_Pos = w_pos / a_pos; end
            if a_neg > 0, Centroid_Neg = w_neg / a_neg; end
        end
        
        function area = bboxArea(theta, u0, v0)
            c = cos(theta); s = sin(theta);
            u_rot = c*u0 + s*v0;
            v_rot = -s*u0 + c*v0;
            area = (max(u_rot) - min(u_rot)) * (max(v_rot) - min(v_rot));
        end
        
        function err = skeletonFitError(R, u, v, X_lim, Y_lim)
            cx = X_lim - R; cy = Y_lim - R;
            dist = zeros(size(u));
            mask_top = (u <= cx);
            dist(mask_top) = abs(v(mask_top) - Y_lim);
            mask_right = (v <= cy) & (~mask_top);
            dist(mask_right) = abs(u(mask_right) - X_lim);
            mask_arc = (~mask_top) & (~mask_right);
            if any(mask_arc)
                dx = u(mask_arc) - cx; dy = v(mask_arc) - cy;
                dist(mask_arc) = abs(sqrt(dx.^2 + dy.^2) - R);
            end
            err = sqrt(mean(dist.^2));
        end
        
        %% --- 内部辅助函数 ---
        function [pts, unique_node_indices] = getPts(mesh, regionID)
            mask = (mesh.RegionTags == regionID);
            if ~any(mask), error('Region %d is empty!', regionID); end
            relevant_T = mesh.T(:, mask);
            unique_node_indices = unique(relevant_T(:));
            pts = mesh.P(:, unique_node_indices);
        end
        
        function [pts_centered, pts] = getCenteredPts(mesh, regionID, center)
            [pts, ~] = CoilGeometryUtils.getPts(mesh, regionID);
            pts_centered = pts - center;
        end
        
        function centers = getElementCenters(mesh)
            P = mesh.P; T = mesh.T;
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
        end
        
        function [u_idx, v_idx] = getUVIndices(axis_idx)
            if axis_idx == 1, u_idx = 3; v_idx = 2;
            elseif axis_idx == 2, u_idx = 1; v_idx = 3;
            else, u_idx = 1; v_idx = 2; 
            end
        end
    end
end