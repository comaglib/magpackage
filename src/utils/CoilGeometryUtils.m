classdef CoilGeometryUtils
    % COILGEOMETRYUTILS 线圈几何辅助工具 (v2.0 - Auto Detection)
    % 
    % 更新:
    %   新增 autoDetectRacetrack 方法，通过 PCA 自动分析网格点云，
    %   提取线圈中心、法向、长短轴尺寸，实现参数自适应。
    
    methods (Static)
        function [center, L_straight, R_curve, normal_axis_idx] = autoDetectRacetrack(mesh, regionID)
            % AUTODETECTRACETRACK 自适应识别跑道型线圈参数
            %
            % 输入:
            %   mesh     - 网格对象
            %   regionID - 线圈区域 ID
            %
            % 输出:
            %   center   - 线圈中心 [x, y, z]
            %   L_straight - 直段长度
            %   R_curve    - 圆角半径 (平均值)
            %   normal_axis_idx - 法线轴索引 (1=X, 2=Y, 3=Z)
            
            fprintf('   [CoilUtils] Auto-detecting geometry for Region %d...\n', regionID);
            
            % 1. 提取线圈区域的所有单元重心
            % 为了提高精度，直接使用节点坐标可能更好，但重心足以
            P = mesh.P;
            T = mesh.T;
            
            % 筛选属于 regionID 的单元
            mask = (mesh.RegionTags == regionID);
            if ~any(mask)
                error('Region %d is empty!', regionID);
            end
            
            % 提取相关单元的节点索引
            relevant_T = T(:, mask);
            unique_node_indices = unique(relevant_T(:));
            
            % 获取点云 (Point Cloud)
            pts = P(:, unique_node_indices); % [3 x N]
            
            % 2. 计算几何中心
            center = mean(pts, 2); % [3 x 1]
            
            % 3. 主成分分析 (PCA) 确定方向
            % 去中心化
            pts_centered = pts - center;
            
            % 计算协方差矩阵 [3 x 3]
            C = cov(pts_centered');
            
            % 特征分解 (V 的列是特征向量，D 是特征值)
            [V, D] = eig(C);
            eig_vals = diag(D);
            
            % 排序特征值 (从小到大: Min, Mid, Max)
            % 假设线圈是扁平的跑道型：
            % - 最小特征值对应 "法线" (厚度方向)
            % - 最大特征值对应 "长轴" (直段方向)
            % - 中间特征值对应 "短轴" (宽度方向)
            
            [~, sort_idx] = sort(eig_vals);
            axis_min = V(:, sort_idx(1)); % 法线
            axis_mid = V(:, sort_idx(2)); % 短轴
            axis_max = V(:, sort_idx(3)); % 长轴
            
            % 确定法线最接近哪个全局坐标轴 (X, Y, Z)
            [~, normal_axis_idx] = max(abs(axis_min));
            
            % 4. 投影点云以测量尺寸
            % 投影到 (长轴, 短轴) 平面
            u_coords = axis_max' * pts_centered;
            v_coords = axis_mid' * pts_centered;
            
            % 计算包围盒范围
            % 由于我们要计算的是"平均"电流路径的尺寸，而不是外轮廓
            % 我们假设点云均匀分布在导体截面上，那么平均路径对应分布的中心区域？
            % 不，通常我们需要电流方向场的"骨架"。
            % 简单做法：取点云的 Max - Min，然后减去截面宽度？
            % 
            % 更稳健的做法：假设线圈是对称的，取均方根或平均绝对偏差来估算特征尺寸
            % 对于均匀分布的矩形截面 [-W/2, W/2]，mean(abs(x)) = W/4
            % 因此 W_effective ≈ 4 * mean(abs(coordinate))
            
            % 估算长轴方向的总长度 (Length) 和短轴方向的总宽度 (Width)
            % 使用 4 * mean(abs) 对于环形分布可能不准。
            % 使用 quantile (例如 5% - 95%) 忽略离群点可能更准，或者直接取 max-min
            
            u_span = max(u_coords) - min(u_coords);
            v_span = max(v_coords) - min(v_coords);
            
            fprintf('   [CoilUtils] Point Cloud Span: Long=%.4f, Short=%.4f\n', u_span, v_span);
            
            % 假设是跑道形 (Racetrack):
            % 总长 L_total = L_straight + 2*R_outer
            % 总宽 W_total = 2*R_outer
            % 这计算的是外轮廓。
            
            % 但我们需要的是"平均路径"的尺寸用于计算切向向量
            % 平均半径 R_avg ≈ R_outer - Thickness/2 ??
            % 这是一个近似。通常取几何中心分布的范围。
            
            % 让我们尝试用简单的 bounding box 估算平均路径：
            % 假设截面较小，可以用 bounding box 代表路径
            L_path_approx = u_span;
            W_path_approx = v_span;
            
            % 如果截面很厚，上述值偏大。
            % 修正策略：我们不仅需要 bounding box，还需要截面厚度。
            % 法线方向的跨度即为厚度 T_thick
            w_coords = axis_min' * pts_centered;
            thickness = max(w_coords) - min(w_coords);
            
            % 假设线圈截面是正方形或圆形，宽度 ~= 厚度
            % 修正后的路径尺寸 (减去一半截面宽度的影响)
            % R_avg = (W_total - thickness) / 2
            R_curve = (W_path_approx - thickness) / 2;
            
            % L_straight = L_total - 2*R_outer = L_total - W_total
            L_straight = L_path_approx - W_path_approx;
            
            % 保护措施
            if L_straight < 0
                L_straight = 0; % 可能是圆形线圈
                R_curve = L_path_approx / 2;
            end
            
            fprintf('   [CoilUtils] Detected: Center=[%.3f %.3f %.3f], L_str=%.4f, R=%.4f, Axis=%d\n', ...
                center(1), center(2), center(3), L_straight, R_curve, normal_axis_idx);
        end
        
        function dir_map = computeRacetrackDirection(mesh, regionID, center, straight_len, radius, axis_idx)
            % (保持原有的计算逻辑不变，直接复制之前的实现)
            % ... [Previous Code for computeRacetrackDirection] ...
            
            numElems = mesh.NumElements;
            dir_map = zeros(3, numElems);
            P = mesh.P; T = mesh.T;
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
            tags = mesh.RegionTags;
            
            if axis_idx == 1, u_idx = 3; v_idx = 2; 
            elseif axis_idx == 2, u_idx = 1; v_idx = 3; 
            elseif axis_idx == 3, u_idx = 1; v_idx = 2; end
            
            for e = 1:numElems
                if tags(e) ~= regionID, continue; end
                p = centers(:, e) - center(:);
                u = p(u_idx); v = p(v_idx);
                
                if abs(u) <= straight_len/2
                    tangent_u = -sign(v); tangent_v = 0;
                else
                    u_shift = u - sign(u) * (straight_len/2); v_shift = v;
                    mag = sqrt(u_shift^2 + v_shift^2);
                    if mag > 0
                        tu = -v_shift / mag; tv = u_shift / mag;
                        tangent_u = tu; tangent_v = tv;
                    else
                        tangent_u = 0; tangent_v = 0;
                    end
                end
                
                vec = zeros(3, 1);
                vec(u_idx) = tangent_u; vec(v_idx) = tangent_v;
                dir_map(:, e) = vec;
            end
        end
    end
end