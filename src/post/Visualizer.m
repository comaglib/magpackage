classdef Visualizer < handle
    % VISUALIZER 有限元可视化模块 (v3.2 - Wireframe Fix)
    % 
    % 描述:
    %   该类封装了基于 MATLAB Graphics 的有限元绘图功能。
    %   它利用 PostProcessor 计算的数据，提供网格显示、场量云图、
    %   矢量图、以及切片图等多种可视化手段。
    %
    % 更新:
    %   1. 修复 plotMesh 的 wireframe 模式: 改为仅绘制外表面线框 (Outline)，
    %      通过提取 freeBoundary 消除了内部海量线条的视觉干扰。
    
    properties
        PostProcessor % 关联的后处理对象，用于获取 Mesh 和计算场量
    end
    
    methods
        function obj = Visualizer(postProcessor)
            obj.PostProcessor = postProcessor;
        end
        
        function plotMesh(obj, varargin)
            % PLOTMESH 绘制网格 (支持 Surface 实体模式和 Corrected Wireframe 线框模式)
            %
            % 参数 (Name-Value Pairs):
            %   'RegionTags': 数组，指定要绘制的区域 ID (默认为空，绘制所有)
            %   'FaceColor' : 颜色，[r g b]
            %   'EdgeColor' : 边线颜色
            %   'Alpha'     : 透明度
            %   'Style'     : 'surface' (默认) 或 'wireframe'
            
            p = inputParser;
            addParameter(p, 'RegionTags', [], @isnumeric);
            addParameter(p, 'FaceColor', [0.8 0.9 1.0]);
            addParameter(p, 'EdgeColor', 'k');
            addParameter(p, 'Alpha', 0.5, @isnumeric);
            addParameter(p, 'Style', 'surface', @ischar); 
            parse(p, varargin{:});
            
            targetTags = p.Results.RegionTags;
            fColor = p.Results.FaceColor;
            eColor = p.Results.EdgeColor;
            fAlpha = p.Results.Alpha;
            style  = p.Results.Style;
            
            mesh = obj.PostProcessor.Mesh;
            
            % 1. 筛选单元
            if isempty(targetTags)
                T_subset = mesh.T; 
            else
                mask = ismember(mesh.RegionTags, targetTags);
                if ~any(mask)
                    warning('Visualizer:NoElements', 'No elements found for the specified tags.');
                    return;
                end
                T_subset = mesh.T(:, mask);
            end
            
            % 2. 节点重索引 (Re-indexing)
            % 如果只绘制网格的一部分，直接引用 mesh.P 会导致产生大量未引用的悬空节点，
            % 这会导致 triangulation 计算量变大或绘图异常。
            % 此处建立局部节点索引映射。
            unique_node_idx = unique(T_subset(:));
            max_id = max(unique_node_idx);
            
            % 建立 Local ID 映射
            map_nodes = zeros(max_id, 1);
            map_nodes(unique_node_idx) = 1:length(unique_node_idx);
            
            % 提取并重组数据
            nodes_subset = mesh.P(:, unique_node_idx)'; % [N_sub x 3]
            elements_local = map_nodes(T_subset);       % [4 x Ne]
            
            % 3. 构建 Triangulation 并提取外表面
            % 3D网格内部面是看不到的，为了提高性能和视觉效果，只提取外表面 (Boundary)
            tr = triangulation(elements_local', nodes_subset);
            
            try
                boundaryFaces = freeBoundary(tr);
            catch ME
                warning('Visualizer:BoundaryFailed', ...
                    '无法提取外表面 (可能是非流形网格): %s', ME.message);
                return;
            end
            
            % 4. 绘图分支
            if strcmpi(style, 'wireframe')
                % [Fix] 仅绘制外表面的线框 (trimesh of boundaryFaces)
                % 将 FaceColor 设为 'none' 使其透明，只显示 EdgeColor
                trimesh(boundaryFaces, nodes_subset(:,1), nodes_subset(:,2), nodes_subset(:,3), ...
                    'FaceColor', 'none', ...
                    'EdgeColor', eColor, ...
                    'EdgeAlpha', fAlpha, ...
                    'LineWidth', 0.8); % 线宽可微调
            else
                % 默认 Surface 模式 (实体表面)
                trisurf(boundaryFaces, nodes_subset(:,1), nodes_subset(:,2), nodes_subset(:,3), ...
                    'FaceColor', fColor, ...
                    'FaceAlpha', fAlpha, ...
                    'EdgeColor', eColor, ...
                    'EdgeAlpha', 0.3); % 表面模式下边线淡一点
            end
            
            % 5. 视图优化
            axis equal; grid on; box on; axis vis3d;
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
            
            % 光照设置 (可选)
            % if strcmpi(style, 'surface')
            %     light('Position', [1 1 1], 'Style', 'infinite');
            %     lighting gouraud; material dull;
            % end
            
            view(3);
        end
        
        function plotFieldOnSurface(obj, regionTags, nodeValues, varargin)
            % PLOTFIELDONSURFACE 在指定区域表面绘制标量场云图
            % 
            % 输入:
            %   regionTags: 目标区域 Tag
            %   nodeValues: 全局节点值数组 (通常来自 PostProcessor.mapElementsToNodes)
            
            p = inputParser;
            addParameter(p, 'FaceAlpha', 1.0, @isnumeric);
            addParameter(p, 'EdgeColor', 'none');
            parse(p, varargin{:});
            fAlpha = p.Results.FaceAlpha;
            eColor = p.Results.EdgeColor;
            
            mesh = obj.PostProcessor.Mesh;
            mask_reg = ismember(mesh.RegionTags, regionTags);
            if ~any(mask_reg), warning('No elements found for tags'); return; end
            
            % 提取子网格并重索引
            T_sub = mesh.T(:, mask_reg);
            unique_nodes = unique(T_sub(:));
            max_id = max(unique_nodes);
            map_nodes = zeros(max_id, 1);
            map_nodes(unique_nodes) = 1:length(unique_nodes);
            
            T_local = map_nodes(T_sub);
            P_sub = mesh.P(:, unique_nodes);
            Val_sub = nodeValues(unique_nodes);
            
            tr = triangulation(T_local', P_sub');
            try
                % 提取外表面进行绘制 (比绘制所有四面体面快得多)
                fb = freeBoundary(tr);
            catch
                warning('Non-manifold mesh detected, plotting all faces.');
                trisurf(T_local', P_sub(1,:), P_sub(2,:), P_sub(3,:), Val_sub, ...
                    'EdgeColor', eColor, 'FaceAlpha', fAlpha);
                return;
            end
            
            % 绘制带颜色的表面
            trisurf(fb, P_sub(1,:), P_sub(2,:), P_sub(3,:), Val_sub, ...
                'EdgeColor', eColor, 'FaceAlpha', fAlpha);
            
            shading interp; % 插值平滑颜色
            colormap('jet'); 
            
            % 设置 Colorbar 并包含最大最小值
            % 默认 colorbar 有时范围不准，这里强制设置
            cb = colorbar;
            v_min = min(Val_sub);
            v_max = max(Val_sub);
            
            % 设置范围
            clim([v_min, v_max]);
            
            % 设置刻度: 保证 Min 和 Max 显示出来，方便读数
            ticks = linspace(v_min, v_max, 7);
            cb.Ticks = sort(unique([v_min, ticks, v_max]));
            
            axis equal; grid on; box on; axis vis3d; shading interp; 
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        end
        
        function plotVectorsOnSurface(obj, regionTags, vectorData, varargin)
            % PLOTVECTORSONSURFACE 在指定区域的外表面绘制矢量场 (Quiver图)
            %
            % 输入:
            %   regionTags: 区域 ID 数组 (例如 [2])
            %   vectorData: 对应区域内所有单元的矢量值 [3 x N_sub]
            %               注意：必须与 find(ismember(mesh.RegionTags, regionTags)) 的顺序一致
            
            p = inputParser;
            addParameter(p, 'MaxArrows', 500, @isnumeric); % 限制箭头最大数量，防卡顿
            addParameter(p, 'Color', 'k');
            addParameter(p, 'Scale', 0.5, @isnumeric);
            parse(p, varargin{:});
            
            n_arrows = p.Results.MaxArrows;
            c_arrow = p.Results.Color;
            s_arrow = p.Results.Scale;
            
            mesh = obj.PostProcessor.Mesh;
            
            % 1. 提取区域单元
            mask_reg = ismember(mesh.RegionTags, regionTags);
            if ~any(mask_reg), warning('No elements found for tags'); return; end
            
            T_sub = mesh.T(:, mask_reg);
            P = mesh.P;
            num_sub = size(T_sub, 2);
            
            if size(vectorData, 2) ~= num_sub
                error('Vector data columns (%d) must match number of elements in regions (%d)', ...
                    size(vectorData, 2), num_sub);
            end
            
            % 2. 提取所有面并构建索引映射
            % 我们需要找到哪些面是外表面，并找到该面所属的单元
            F_all = [T_sub(1,:)' T_sub(2,:)' T_sub(3,:)';
                     T_sub(1,:)' T_sub(2,:)' T_sub(4,:)';
                     T_sub(1,:)' T_sub(3,:)' T_sub(4,:)';
                     T_sub(2,:)' T_sub(3,:)' T_sub(4,:)'];
                 
            % 记录每个面来自哪个局部单元索引 (1...num_sub)
            elem_indices_map = repmat(1:num_sub, 4, 1);
            elem_indices_map = elem_indices_map(:);
            
            % 3. 查找边界面 (只出现一次的面)
            F_sorted = sort(F_all, 2);
            [~, ~, ic] = unique(F_sorted, 'rows');
            counts = accumarray(ic, 1);
            
            % 只保留外表面 (出现次数为1)
            mask_boundary = (counts(ic) == 1);
            
            boundary_elem_indices = elem_indices_map(mask_boundary);
            boundary_faces = F_all(mask_boundary, :);
            
            if isempty(boundary_faces), return; end
            
            % 4. 计算面重心
            C1 = P(:, boundary_faces(:,1));
            C2 = P(:, boundary_faces(:,2));
            C3 = P(:, boundary_faces(:,3));
            centers = (C1 + C2 + C3) / 3.0; 
            
            % 5. 提取对应单元的矢量
            vecs = vectorData(:, boundary_elem_indices);
            
            % 6. 随机降采样 (避免箭头太密看不清)
            num_faces = size(centers, 2);
            if num_faces > n_arrows
                perm = randperm(num_faces, n_arrows);
                draw_centers = centers(:, perm);
                draw_vecs = vecs(:, perm);
            else
                draw_centers = centers;
                draw_vecs = vecs;
            end
            
            % 7. 绘图 (稍微沿法向抬高一点通常更好，这里简单抬高Z)
            % 偏移量是为了防止箭头被表面遮挡
            offset = 1e-4; 
            
            quiver3(draw_centers(1,:), draw_centers(2,:), draw_centers(3,:) + offset, ...
                    draw_vecs(1,:), draw_vecs(2,:), draw_vecs(3,:), ...
                    s_arrow, 'Color', c_arrow, 'LineWidth', 1.0);
        end
        
        function plotLineComparison(~, simStruct, refStruct, metaStruct)
            % PLOTLINECOMPARISON 绘制 1D 仿真与基准数据的对比图 (实部/虚部)
            %
            % 输入:
            %   simStruct: 包含 .x, .re, .im (仿真曲线数据)
            %   refStruct: 包含 .x, .re, .im (基准散点数据，可为空)
            %   metaStruct: 包含 .title, .ylabel, .xlabel
            
            % 1. 绘制仿真曲线 (实线/虚线)
            h1 = plot(simStruct.x, simStruct.re, 'b-', 'LineWidth', 1.5); hold on;
            h2 = plot(simStruct.x, simStruct.im, 'r--', 'LineWidth', 1.5);
            
            % 2. 绘制基准点 (圆圈/叉号)
            if ~isempty(refStruct)
                h3 = plot(refStruct.x, refStruct.re, 'bo', 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
                h4 = plot(refStruct.x, refStruct.im, 'rx', 'LineWidth', 1.0);
                legends = {'Sim (Re)', 'Sim (Im)', 'Exp (Re)', 'Exp (Im)'};
                handles = [h1, h2, h3, h4];
            else
                legends = {'Sim (Re)', 'Sim (Im)'};
                handles = [h1, h2];
            end
            
            % 3. 装饰
            grid on;
            if isfield(metaStruct, 'xlabel'), xlabel(metaStruct.xlabel); else, xlabel('x'); end
            if isfield(metaStruct, 'ylabel'), ylabel(metaStruct.ylabel); end
            if isfield(metaStruct, 'title'), title(metaStruct.title); end
            
            legend(handles, legends, 'Location', 'best', 'FontSize', 8);
            hold off;
        end
        
        % OTHERS
        function plotSurfaceWithProbes(obj, targetZ, probePts, regionTag, tol)
            % PLOTSURFACEWITHPROBES 绘制特定高度的平面网格并叠加探测点
            % 
            % 输入:
            %   targetZ: 目标 Z 平面坐标 (如 0.019)
            %   probePts: 探测点坐标矩阵 [N x 3]
            %   regionTag: (可选) 限制区域 ID，默认 2 (通常是导体)
            %   tol: (可选) Z轴容差，默认 1mm (用于截取平面)
            
            if nargin < 5, tol = 1e-3; end
            if nargin < 4, regionTag = 2; end
            
            mesh = obj.PostProcessor.Mesh;
            
            % 1. 提取区域网格
            mask_reg = (mesh.RegionTags == regionTag);
            T_sub = mesh.T(:, mask_reg);
            
            % 2. 提取所有面并筛选 Z 面
            % (简单的面提取，不做严格的拓扑去重以提高速度，绘图会自动处理重叠)
            faces = [T_sub(1,:)' T_sub(2,:)' T_sub(3,:)';
                     T_sub(1,:)' T_sub(2,:)' T_sub(4,:)';
                     T_sub(1,:)' T_sub(3,:)' T_sub(4,:)';
                     T_sub(2,:)' T_sub(3,:)' T_sub(4,:)'];
            
            P = mesh.P;
            z_nodes = P(3, :);
            % 计算面中心Z坐标
            z_faces = (z_nodes(faces(:,1)) + z_nodes(faces(:,2)) + z_nodes(faces(:,3))) / 3.0;
            
            % 筛选接近 targetZ 的面
            mask_surf = abs(z_faces - targetZ) < tol;
            surf_faces = faces(mask_surf, :);
            
            if isempty(surf_faces)
                warning('Visualizer:NoSurface', 'No surface found at Z=%.3f +/- %.3f', targetZ, tol);
                return;
            end
            
            % 3. 绘图
            patch('Faces', surf_faces, 'Vertices', P', ...
                  'FaceColor', [0.8 0.9 1.0], 'EdgeColor', [0.5 0.5 0.5], 'FaceAlpha', 0.8);
            hold on;
            
            % 4. 叠加探测点
            if ~isempty(probePts)
                % 将点稍微抬高一点以免被面遮挡
                plot3(probePts(:,1), probePts(:,2), probePts(:,3) + tol/10, ...
                      'r.', 'MarkerSize', 15, 'LineWidth', 2);
            end
            
            axis equal; grid on; view(2);
            xlabel('X (m)'); ylabel('Y (m)');
            title(sprintf('Surface Mesh (Z=%.3f) & Probes', targetZ));
            hold off;
        end
        
        function plotElementVectors(obj, vectorData, varargin)
            % PLOTELEMENTVECTORS 绘制单元矢量 (基于单元中心)
            p = inputParser;
            addParameter(p, 'MaxArrows', 500, @isnumeric);
            addParameter(p, 'Color', 'r');
            addParameter(p, 'Scale', 0.5, @isnumeric);
            parse(p, varargin{:});
            
            n_arrows = p.Results.MaxArrows;
            c_arrow = p.Results.Color;
            s_arrow = p.Results.Scale;
            
            mesh = obj.PostProcessor.Mesh;
            
            % 筛选非零矢量 (避免绘制零向量)
            mags = sum(vectorData.^2, 1);
            valid_indices = find(mags > 1e-12);
            
            if isempty(valid_indices), return; end
            
            % 降采样
            if length(valid_indices) > n_arrows
                perm = randperm(length(valid_indices), n_arrows);
                draw_indices = valid_indices(perm);
            else
                draw_indices = valid_indices;
            end
            
            % 计算重心
            draw_T = mesh.T(:, draw_indices);
            centers = (mesh.P(:, draw_T(1,:)) + ...
                       mesh.P(:, draw_T(2,:)) + ...
                       mesh.P(:, draw_T(3,:)) + ...
                       mesh.P(:, draw_T(4,:))) / 4.0;
            
            vecs = vectorData(:, draw_indices);
            
            quiver3(centers(1,:), centers(2,:), centers(3,:), ...
                    vecs(1,:), vecs(2,:), vecs(3,:), ...
                    s_arrow, 'LineWidth', 1.5, 'Color', c_arrow, 'MaxHeadSize', 0.5);
            
            axis equal vis3d; grid on; box on; view(3);
        end
        
        function plotVectorField(obj, A_sol, points)
            % PLOTVECTORFIELD 绘制给定点集的矢量场 (通过 Probe 计算)
            N = size(points, 1);
            B_vals = zeros(N, 3);
            
            for i = 1:N
                val = obj.PostProcessor.probeB(A_sol, points(i,:));
                B_vals(i,:) = real(val'); 
            end
            
            quiver3(points(:,1), points(:,2), points(:,3), ...
                    B_vals(:,1), B_vals(:,2), B_vals(:,3), 1.5);
            axis equal; grid on;
            axis vis3d;
            
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Magnetic Flux Density B (Vectors)');
            view(3);
        end
        
        function plotSliceSmoothed(obj, A_sol, normal, offset, gridDensity)
            % PLOTSLICESMOOTHED 绘制平滑的 2D 切片云图
            %
            % 输入:
            %   normal: 切片法向 ('x', 'y', 'z')
            %   offset: 截距距离
            %   gridDensity: 插值网格密度
            
            if nargin < 5, gridDensity = 100; end
            
            % 1. 先计算全场单元值并平滑到节点
            B_elems = obj.PostProcessor.computeElementB(A_sol);
            B_nodes = obj.PostProcessor.mapElementsToNodes(B_elems);
            B_mag_nodes = sqrt(sum(abs(B_nodes).^2, 1)); 
            
            % 2. 准备插值器 scatteredInterpolant
            X = obj.PostProcessor.Mesh.P(1,:)';
            Y = obj.PostProcessor.Mesh.P(2,:)';
            Z = obj.PostProcessor.Mesh.P(3,:)';
            V = B_mag_nodes(:);
            
            F = scatteredInterpolant(X, Y, Z, V, 'natural', 'none');
            mesh_min = min(obj.PostProcessor.Mesh.P, [], 2);
            mesh_max = max(obj.PostProcessor.Mesh.P, [], 2);
            
            % 3. 生成切片网格
            switch lower(normal)
                case 'z', x_r=linspace(mesh_min(1),mesh_max(1),gridDensity); y_r=linspace(mesh_min(2),mesh_max(2),gridDensity); [G1,G2]=meshgrid(x_r,y_r); G3=ones(size(G1))*offset; xl='X'; yl='Y'; p1=G1; p2=G2;
                case 'y', x_r=linspace(mesh_min(1),mesh_max(1),gridDensity); z_r=linspace(mesh_min(3),mesh_max(3),gridDensity); [G1,G3]=meshgrid(x_r,z_r); G2=ones(size(G1))*offset; xl='X'; yl='Z'; p1=G1; p2=G3;
                case 'x', y_r=linspace(mesh_min(2),mesh_max(2),gridDensity); z_r=linspace(mesh_min(3),mesh_max(3),gridDensity); [G2,G3]=meshgrid(y_r,z_r); G1=ones(size(G2))*offset; xl='Y'; yl='Z'; p1=G2; p2=G3;
            end
            
            % 4. 插值并绘制
            V_grid = F(G1, G2, G3);
            surf(p1, p2, V_grid, 'EdgeColor', 'none');
            view(2); axis equal tight; shading interp; colormap('jet'); colorbar; xlabel(xl); ylabel(yl);
        end
    end
end