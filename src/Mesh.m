classdef Mesh < handle
    % MESH 有限元网格基础类 (v2.1 - Volume Calculation)
    % 
    % 更新:
    %   1. load 方法现在支持 unit_tag 参数 (例如 'mm')，用于坐标缩放。
    %   2. [New] 增加 CellVolumes 属性及自动计算逻辑。
    
    properties
        P           % 节点坐标 [3 x N_nodes]
        T           % 单元连接关系 [N_nodes_per_elem x N_elems] (Local Node IDs)
        RegionTags  % 单元所属区域 ID [1 x N_elems]
        
        % 拓扑实体
        Edges       % [2 x N_edges] 存储每条棱的两个端点 (全局ID)
        Faces       % [3 x N_faces]
        
        % 拓扑映射
        T2E         % [6 x N_elems] 单元到棱的映射
        T2E_Sign    % [6 x N_elems] 局部棱与全局棱的方向一致性
        
        % [New] 几何属性
        CellVolumes % [1 x N_elems] 单元体积 (m^3)
        
        NumNodes
        NumElements
    end
    
    methods (Static)
        function obj = load(filename, unit_tag)
            % LOAD 从文件加载网格
            % 输入:
            %   filename - 网格文件路径 (.mphtxt 或 .msh)
            %   unit_tag - (可选) 长度单位，如 'mm', 'cm', 'm' (默认)
            
            if nargin < 2, unit_tag = 'm'; end
            
            obj = Mesh();
            [~, ~, ext] = fileparts(filename);
            
            fprintf('[Mesh] Loading %s (Unit: %s) ...\n', filename, unit_tag);
            
            if strcmpi(ext, '.msh')
                raw = read_msh(filename); 
                if ~strcmpi(unit_tag, 'm')
                    warning('read_msh 暂不支持自动缩放，请手动确认单位。');
                end
            elseif strcmpi(ext, '.mphtxt')
                % 传递 unit_tag 给 read_mphtxt
                raw = read_mphtxt(filename, unit_tag);
            else
                error('Unsupported mesh format: %s', ext);
            end
            
            obj.P = raw.P;
            obj.T = raw.T;
            obj.RegionTags = raw.RegionTags;
            
            obj.NumNodes = size(obj.P, 2);
            obj.NumElements = size(obj.T, 2);
            
            % [New] 加载后立即计算单元体积
            obj.computeCellVolumes();
        end
    end
    
    methods
        function computeCellVolumes(obj)
            % COMPUTECELLVOLUMES 计算所有单元的体积
            if isempty(obj.P) || isempty(obj.T), return; end
            
            fprintf('[Mesh] Computing element volumes...\n');
            
            dim = size(obj.P, 1);
            nodesPerElem = size(obj.T, 1);
            
            if dim == 3 && nodesPerElem >= 4
                % 3D 四面体体积计算
                % V = |det(v1, v2, v3)| / 6
                % 向量化计算: V = |dot(cross(v1, v2), v3)| / 6
                
                % 提取四个顶点的坐标
                p1 = obj.P(:, obj.T(1,:));
                p2 = obj.P(:, obj.T(2,:));
                p3 = obj.P(:, obj.T(3,:));
                p4 = obj.P(:, obj.T(4,:));
                
                v1 = p2 - p1;
                v2 = p3 - p1;
                v3 = p4 - p1;
                
                % cross(v1, v2, 1) 表示沿第一维度(行)做叉乘
                cp = cross(v1, v2, 1);
                
                % 点乘后取绝对值，除以 6
                obj.CellVolumes = abs(sum(cp .* v3, 1)) / 6.0;
                
            elseif dim == 2 && nodesPerElem >= 3
                % 2D 三角形面积计算 (作为体积处理)
                p1 = obj.P(:, obj.T(1,:));
                p2 = obj.P(:, obj.T(2,:));
                p3 = obj.P(:, obj.T(3,:));
                
                v1 = p2 - p1;
                v2 = p3 - p1;
                
                % 2D 叉乘模长: |x1*y2 - x2*y1|
                obj.CellVolumes = 0.5 * abs(v1(1,:).*v2(2,:) - v1(2,:).*v2(1,:));
            else
                warning('Mesh:Volume', 'Volume calculation not implemented for this element type. CellVolumes will be empty.');
                obj.CellVolumes = [];
            end
        end
        
        function generateEdges(obj)
            if ~isempty(obj.Edges), return; end
            
            fprintf('[Mesh] Generating Edge topology...\n');
            local_edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
            idx1 = local_edge_defs(:, 1);
            idx2 = local_edge_defs(:, 2);
            
            start_nodes = obj.T(idx1, :); 
            end_nodes   = obj.T(idx2, :); 
            
            raw_edges = [start_nodes(:), end_nodes(:)];
            raw_edges_sorted = sort(raw_edges, 2);
            
            [unique_edges, ~, ic] = unique(raw_edges_sorted, 'rows');
            
            obj.Edges = unique_edges';
            numEdges = size(obj.Edges, 2);
            obj.T2E = reshape(ic, 6, obj.NumElements);
            
            is_same = (raw_edges(:,1) == raw_edges_sorted(:,1));
            obj.T2E_Sign = reshape(int8(is_same)*2 - 1, 6, obj.NumElements);
            
            fprintf('[Mesh] Topology built: %d Edges found.\n', numEdges);
        end
        
        function generateFaces(obj)
            % GENERATEFACES 生成四面体网格的面拓扑
            if ~isempty(obj.Faces) && size(obj.Faces, 2) > size(obj.T, 2), return; end
            
            fprintf('[Mesh] Generating Face topology...\n');
            local_face_defs = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
            
            F1 = obj.T(local_face_defs(1,:), :)';
            F2 = obj.T(local_face_defs(2,:), :)';
            F3 = obj.T(local_face_defs(3,:), :)';
            F4 = obj.T(local_face_defs(4,:), :)';
            
            all_raw_faces = [F1; F2; F3; F4];
            all_sorted_faces = sort(all_raw_faces, 2);
            unique_faces = unique(all_sorted_faces, 'rows');
            
            obj.Faces = unique_faces';
            fprintf('[Mesh] Faces generated: %d unique faces found.\n', size(obj.Faces, 2));
        end
        
        function stats(obj)
            fprintf('--- Mesh Statistics ---\n');
            fprintf('Nodes: %d\n', obj.NumNodes);
            fprintf('Elements: %d\n', obj.NumElements);
            if ~isempty(obj.Edges)
                fprintf('Edges: %d\n', size(obj.Edges, 2));
            end
            if ~isempty(obj.Faces)
                fprintf('Faces: %d\n', size(obj.Faces, 2));
            end
            if ~isempty(obj.CellVolumes)
                fprintf('Total Volume: %.4e m^3\n', sum(obj.CellVolumes));
            end
        end
    end
end