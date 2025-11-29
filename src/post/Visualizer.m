classdef Visualizer < handle
    % VISUALIZER 有限元可视化模块 (v3.2 - Wireframe Fix)
    % 
    % 更新:
    %   1. 修复 plotMesh 的 wireframe 模式: 改为仅绘制外表面线框 (Outline)，
    %      消除了内部线条的视觉干扰。
    
    properties
        PostProcessor
    end
    
    methods
        function obj = Visualizer(postProcessor)
            obj.PostProcessor = postProcessor;
        end
        
        function plotMesh(obj, varargin)
            % PLOTMESH 绘制网格 (支持 Surface 和 Corrected Wireframe)
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
            
            % 2. 节点重索引 (Re-indexing) - 避免未引用节点警告
            unique_node_idx = unique(T_subset(:));
            max_id = max(unique_node_idx);
            
            % 建立 Local ID 映射
            map_nodes = zeros(max_id, 1);
            map_nodes(unique_node_idx) = 1:length(unique_node_idx);
            
            % 提取并重组数据
            nodes_subset = mesh.P(:, unique_node_idx)'; % [N_sub x 3]
            elements_local = map_nodes(T_subset);       % [4 x Ne]
            
            % 3. 构建 Triangulation 并提取外表面
            % 注意: triangulation 输入需为 [Ne x 4]
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
                % FaceColor='none' 使其透明，只显示 EdgeColor
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
            
            if strcmpi(style, 'surface')
                % 仅在表面模式启用光照
                light('Position', [1 1 1], 'Style', 'infinite');
                lighting gouraud; material dull;
            end
            
            view(3);
        end
        
        function plotElementVectors(obj, vectorData, varargin)
            % PLOTELEMENTVECTORS 绘制单元矢量
            p = inputParser;
            addParameter(p, 'MaxArrows', 500, @isnumeric);
            addParameter(p, 'Color', 'r');
            addParameter(p, 'Scale', 0.5, @isnumeric);
            parse(p, varargin{:});
            
            n_arrows = p.Results.MaxArrows;
            c_arrow = p.Results.Color;
            s_arrow = p.Results.Scale;
            
            mesh = obj.PostProcessor.Mesh;
            numElems = mesh.NumElements;
            
            % 筛选非零矢量
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
            % PLOTVECTORFIELD 绘制矢量场
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
            % (保持原有代码不变)
            if nargin < 5, gridDensity = 100; end
            B_elems = obj.PostProcessor.computeElementB(A_sol);
            B_nodes = obj.PostProcessor.mapElementsToNodes(B_elems);
            B_mag_nodes = sqrt(sum(abs(B_nodes).^2, 1)); 
            
            X = obj.PostProcessor.Mesh.P(1,:)';
            Y = obj.PostProcessor.Mesh.P(2,:)';
            Z = obj.PostProcessor.Mesh.P(3,:)';
            V = B_mag_nodes(:);
            
            F = scatteredInterpolant(X, Y, Z, V, 'natural', 'none');
            mesh_min = min(obj.PostProcessor.Mesh.P, [], 2);
            mesh_max = max(obj.PostProcessor.Mesh.P, [], 2);
            
            switch lower(normal)
                case 'z', x_r=linspace(mesh_min(1),mesh_max(1),gridDensity); y_r=linspace(mesh_min(2),mesh_max(2),gridDensity); [G1,G2]=meshgrid(x_r,y_r); G3=ones(size(G1))*offset; xl='X'; yl='Y'; p1=G1; p2=G2;
                case 'y', x_r=linspace(mesh_min(1),mesh_max(1),gridDensity); z_r=linspace(mesh_min(3),mesh_max(3),gridDensity); [G1,G3]=meshgrid(x_r,z_r); G2=ones(size(G1))*offset; xl='X'; yl='Z'; p1=G1; p2=G3;
                case 'x', y_r=linspace(mesh_min(2),mesh_max(2),gridDensity); z_r=linspace(mesh_min(3),mesh_max(3),gridDensity); [G2,G3]=meshgrid(y_r,z_r); G1=ones(size(G2))*offset; xl='Y'; yl='Z'; p1=G2; p2=G3;
            end
            
            V_grid = F(G1, G2, G3);
            surf(p1, p2, V_grid, 'EdgeColor', 'none');
            view(2); axis equal tight; shading interp; colormap('jet'); colorbar; xlabel(xl); ylabel(yl);
        end
    end
end