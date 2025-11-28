classdef Mesh < handle
    % MESH 有限元网格基础类
    % 管理节点(P)、单元(T)以及拓扑关系(Edges, Faces)
    
    properties
        P           % 节点坐标 [3 x N_nodes]
        T           % 单元连接关系 [N_nodes_per_elem x N_elems] (Local Node IDs)
        RegionTags  % 单元所属区域 ID [1 x N_elems]
        
        % 拓扑实体 (Lazy Loading: 只有在需要时才计算)
        Edges       % [2 x N_edges] 存储每条棱的两个端点 (全局ID)
        Faces       % [3 x N_faces] (暂预留)
        
        % 拓扑映射
        T2E         % [6 x N_elems] 单元到棱的映射 (Edge IDs)
        T2E_Sign    % [6 x N_elems] 局部棱与全局棱的方向一致性 (+1/-1)
        
        NumNodes
        NumElements
    end
    
    methods (Static)
        function obj = load(filename)
            % 工厂方法：根据后缀自动选择读取器
            obj = Mesh();
            [~, ~, ext] = fileparts(filename);
            
            fprintf('[Mesh] Loading %s ...\n', filename);
            if strcmpi(ext, '.msh')
                raw = read_msh(filename); % 调用原有的读取函数
            elseif strcmpi(ext, '.mphtxt')
                raw = read_mphtxt(filename);
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
            % 生成棱拓扑信息 (移植自原 build_topology.m)
            % 如果已经生成过，则跳过
            if ~isempty(obj.Edges), return; end
            
            fprintf('[Mesh] Generating Edge topology...\n');
            
            % 1. 定义局部棱 (Tetrahedron: 6 edges)
            % 节点对: [1-2, 1-3, 1-4, 2-3, 2-4, 3-4]
            local_edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
            
            % 2. 向量化构建所有潜在棱
            idx1 = local_edge_defs(:, 1);
            idx2 = local_edge_defs(:, 2);
            
            start_nodes = obj.T(idx1, :); % 6 x Ne
            end_nodes   = obj.T(idx2, :); % 6 x Ne
            
            % 3. 排序并去重
            % 确保每条棱都是 min_id -> max_id
            raw_edges = [start_nodes(:), end_nodes(:)];
            raw_edges_sorted = sort(raw_edges, 2);
            
            [unique_edges, ~, ic] = unique(raw_edges_sorted, 'rows');
            
            % 4. 存储全局棱
            obj.Edges = unique_edges';
            numEdges = size(obj.Edges, 2);
            
            % 5. 构建映射表 (重塑 ic)
            obj.T2E = reshape(ic, 6, obj.NumElements);
            
            % 6. 计算方向符号
            % 原方向(raw)与排序后方向(sorted)一致则为+1，否则-1
            is_same = (raw_edges(:,1) == raw_edges_sorted(:,1));
            obj.T2E_Sign = reshape(int8(is_same)*2 - 1, 6, obj.NumElements);
            
            fprintf('[Mesh] Topology built: %d Edges found.\n', numEdges);
        end
        
        function stats(obj)
            fprintf('--- Mesh Statistics ---\n');
            fprintf('Nodes: %d\n', obj.NumNodes);
            fprintf('Elements: %d\n', obj.NumElements);
            if ~isempty(obj.Edges)
                fprintf('Edges: %d\n', size(obj.Edges, 2));
            else
                fprintf('Edges: Not generated yet.\n');
            end
        end
    end
end