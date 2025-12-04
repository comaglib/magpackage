classdef Mesh < handle
    % MESH 有限元网格基础类 (v2.0 - Unit Support)
    % 
    % 更新:
    %   1. load 方法现在支持 unit_tag 参数 (例如 'mm')，用于坐标缩放。
    
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
                % 注意: read_msh 可能需要更新以支持 unit_tag，这里暂只处理 mphtxt
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
        end
    end
    
    methods
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
            % GENERATEFACES 生成四面体网格的面拓扑 (Triangular Faces)
            %
            % 功能:
            %   从四面体连接关系 (obj.T) 中提取所有唯一的三角形面，
            %   并存储在 obj.Faces 中。
            %   (如果 obj.Faces 已存在且包含数据，该函数可能会覆盖它以确保包含内部面)
            
            % 如果需要强制重新生成，可注释掉下面这行
            if ~isempty(obj.Faces) && size(obj.Faces, 2) > size(obj.T, 2), return; end
            
            fprintf('[Mesh] Generating Face topology...\n');
            
            % 四面体的4个面 (局部节点索引定义)
            % 组合: [1 2 3], [1 2 4], [1 3 4], [2 3 4]
            local_face_defs = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
            
            % 预分配/构造大矩阵
            % obj.T 是 [4 x N_elems]
            
            % 提取所有单元的4个面
            % 为了向量化操作，我们将所有面堆叠起来
            F1 = obj.T(local_face_defs(1,:), :)'; % [N_elems x 3]
            F2 = obj.T(local_face_defs(2,:), :)';
            F3 = obj.T(local_face_defs(3,:), :)';
            F4 = obj.T(local_face_defs(4,:), :)';
            
            % 合并所有面: [4*N_elems x 3]
            all_raw_faces = [F1; F2; F3; F4];
            
            % 对每个面的节点索引进行排序 (行内排序)，以确保唯一性判断不受节点顺序影响
            % 例如 [3 1 2] 和 [1 2 3] 被视为同一个面
            all_sorted_faces = sort(all_raw_faces, 2);
            
            % 使用 unique 函数查找唯一面
            % 'rows' 表示按行比较
            unique_faces = unique(all_sorted_faces, 'rows');
            
            % 转置并存储结果: [3 x N_unique_faces]
            obj.Faces = unique_faces';
            
            % 这里的 NumFaces 是所有面（包括内部面和边界面）的总数
            fprintf('[Mesh] Faces generated: %d unique faces found.\n', size(obj.Faces, 2));
            
            % (可选) 如果您以后需要单元到面的映射 (T2F)，可以在这里扩展:
            % [unique_faces, ~, ic] = unique(...);
            % obj.T2F = reshape(ic, [], 4)';
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
        end
    end
end