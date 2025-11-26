function DofData = create_dof_map(Model)
% CREATE_DOF_MAP 生成优化的自由度映射 (A-V Formulation)
% 
% 策略:
% 1. A-DoFs: 定义在所有棱上 (1 .. nEdges)
% 2. V-DoFs: 仅定义在导体区域的节点上 (nEdges+1 .. nEdges+nActiveNodes)
%
% 输出 DofData 结构体:
%   .NumEdges      - 棱总数
%   .NumNodes      - 节点总数
%   .NumActiveNodes- 导电区节点数
%   .TotalDoFs     - 系统总自由度
%   .Node2DoF      - [N_nodes x 1] 映射表 (NodeID -> Global DoF ID)
%                    如果节点非导电，值为 0

    fprintf('正在构建自由度映射 (Active V-Nodes Only)...\n');
    
    Mesh = Model.Mesh;
    T = Mesh.T;
    numNodes = size(Mesh.P, 2);
    numEdges = size(Mesh.Edges, 2);
    numElems = size(T, 2);
    
    % 1. 识别导电节点
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    IsNodeActive = false(numNodes, 1);
    
    % 遍历所有单元，如果单元导电，其所有节点标记为 Active
    for i = 1:numElems
        mat_id = MatMap(i);
        if mat_id > length(MatLib), continue; end
        
        if MatLib(mat_id).Sigma > 1e-12
            nodes = T(:, i);
            IsNodeActive(nodes) = true;
        end
    end
    
    % 2. 构建映射表
    active_nodes_indices = find(IsNodeActive);
    numActive = length(active_nodes_indices);
    
    Node2DoF = zeros(numNodes, 1);
    
    % V 的自由度从 numEdges + 1 开始
    current_dof = numEdges + 1;
    
    for k = 1:numActive
        nid = active_nodes_indices(k);
        Node2DoF(nid) = current_dof;
        current_dof = current_dof + 1;
    end
    
    % 3. 打包输出
    DofData.NumEdges = numEdges;
    DofData.NumNodes = numNodes;
    DofData.NumActiveNodes = numActive;
    DofData.TotalDoFs = numEdges + numActive;
    DofData.Node2DoF = Node2DoF;
    
    fprintf('  - 系统总自由度: %d (Edges: %d, Active Nodes: %d)\n', ...
        DofData.TotalDoFs, numEdges, numActive);
end