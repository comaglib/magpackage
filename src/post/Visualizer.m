classdef Visualizer < handle
    % VISUALIZER 有限元可视化模块 (v1.0)
    % 
    % 功能:
    %   1. 提供基于插值的云图切片 (Slice Contour)。
    %   2. 将绘图逻辑与计算逻辑分离。
    
    properties
        PostProcessor
    end
    
    methods
        function obj = Visualizer(postProcessor)
            obj.PostProcessor = postProcessor;
        end
        
        function plotSliceMagnitude(obj, A_sol, normal, offset, gridDensity)
            % PLOTSLICEMAGNITUDE 绘制磁通密度模值 (|B|) 的切片云图
            % 
            % 输入:
            %   A_sol       - 解向量 (单列)
            %   normal      - 切片法向 ('x', 'y', 'z')
            %   offset      - 切片位置 (例如 z=0.25)
            %   gridDensity - (可选) 插值网格密度，默认 100
            
            if nargin < 5, gridDensity = 100; end
            
            % 1. 计算全场 B 及其模值
            % 注意: 这里假设 A_sol 是单列向量 (如果是 HBFEM，请传入某一列)
            B_elems = obj.PostProcessor.computeElementB(A_sol);
            B_mag = obj.PostProcessor.computeMagnitude(B_elems); % [Ne x 1]
            
            % 2. 获取几何中心
            centers = obj.PostProcessor.computeCentroids(); % [3 x Ne]
            X = centers(1,:)'; 
            Y = centers(2,:)'; 
            Z = centers(3,:)';
            V = B_mag(:);
            
            % 3. 建立插值器 (Scattered Interpolant)
            % 使用 'linear' 插值和 'none' 外插 (背景设为 NaN/白色)
            F = scatteredInterpolant(X, Y, Z, V, 'linear', 'none');
            
            % 4. 生成切片网格
            mesh_min = min(obj.PostProcessor.Mesh.P, [], 2);
            mesh_max = max(obj.PostProcessor.Mesh.P, [], 2);
            
            switch lower(normal)
                case 'z'
                    x_range = linspace(mesh_min(1), mesh_max(1), gridDensity);
                    y_range = linspace(mesh_min(2), mesh_max(2), gridDensity);
                    [G1, G2] = meshgrid(x_range, y_range);
                    G3 = ones(size(G1)) * offset;
                    xlabel_str = 'X (m)'; ylabel_str = 'Y (m)'; title_str = sprintf('|B| Slice at Z=%.3f', offset);
                    
                case 'y'
                    x_range = linspace(mesh_min(1), mesh_max(1), gridDensity);
                    z_range = linspace(mesh_min(3), mesh_max(3), gridDensity);
                    [G1, G3] = meshgrid(x_range, z_range);
                    G2 = ones(size(G1)) * offset;
                    xlabel_str = 'X (m)'; ylabel_str = 'Z (m)'; title_str = sprintf('|B| Slice at Y=%.3f', offset);
                    
                case 'x'
                    y_range = linspace(mesh_min(2), mesh_max(2), gridDensity);
                    z_range = linspace(mesh_min(3), mesh_max(3), gridDensity);
                    [G2, G3] = meshgrid(y_range, z_range);
                    G1 = ones(size(G2)) * offset;
                    xlabel_str = 'Y (m)'; ylabel_str = 'Z (m)'; title_str = sprintf('|B| Slice at X=%.3f', offset);
                    
                otherwise
                    error('Unknown normal direction: %s', normal);
            end
            
            % 5. 插值计算
            V_grid = F(G1, G2, G3);
            
            % 6. 绘图 (Surface Plot)
            % 根据法向选择绘图轴
            figure;
            if strcmpi(normal, 'z')
                surf(G1, G2, V_grid, 'EdgeColor', 'none');
            elseif strcmpi(normal, 'y')
                surf(G1, G3, V_grid, 'EdgeColor', 'none');
            else
                surf(G2, G3, V_grid, 'EdgeColor', 'none');
            end
            
            view(2); % 2D 视图
            axis equal tight;
            shading interp;
            colormap('jet'); 
            colorbar;
            xlabel(xlabel_str); ylabel(ylabel_str);
            title(title_str);
            
            fprintf('   [Visualizer] Slice plot generated.\n');
        end
    end
end