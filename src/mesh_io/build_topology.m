function Mesh = build_topology(MeshRaw)
% BUILD_TOPOLOGY 构建由节点、单元和棱组成的完整拓扑结构 (SoA)
% 输入:
%   MeshRaw - 包含 P, T, RegionTags 的简单结构体
% 输出:
%   Mesh - 完整的 Model.Mesh 结构体，包含 Edges, T2E, T2E_Sign

    fprintf('正在构建网格拓扑 (Edges & T2E)...\n');
    
    Mesh.P = MeshRaw.P;
    Mesh.T = MeshRaw.T;
    Mesh.RegionTags = MeshRaw.RegionTags;
    
    numElems = size(Mesh.T, 2);
    
    % 1. 定义四面体的6条局部棱 (基于局部节点索引 1..4)
    % Nédélec 单元标准顺序: 
    % e1:1-2, e2:1-3, e3:1-4, e4:2-3, e5:2-4, e6:3-4
    local_edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
    
    % 2. 生成所有单元的所有棱 (未去重)
    % 维度: 2 x (6 * numElems)
    % 策略: 先构造大矩阵，再 unique
    
    % 展开 T 矩阵以便向量化操作
    % T 的行是节点索引
    all_edges_start = Mesh.T(local_edge_defs(:,1), :); % 6 x Ne
    all_edges_end   = Mesh.T(local_edge_defs(:,2), :); % 6 x Ne
    
    all_edges_start = all_edges_start(:); % 6Ne x 1
    all_edges_end   = all_edges_end(:);   % 6Ne x 1
    
    % 3. 规范化棱的方向用于去重
    % 全局棱总是定义为: min_idx -> max_idx
    raw_edges = [all_edges_start, all_edges_end];
    raw_edges_sorted = sort(raw_edges, 2); % 每一行按从小到大排序
    
    % 4. 寻找唯一棱 (Unique)
    % [C, ia, ic] = unique(...)
    % C: 唯一棱列表 (Edges)
    % ic: 原始棱在 C 中的索引，这将直接构成 T2E
    fprintf('  - 正在去重并生成全局棱列表...\n');
    [unique_edges, ~, ic] = unique(raw_edges_sorted, 'rows');
    
    Mesh.Edges = unique_edges'; % 转置为 2 x N_edges
    numEdges = size(Mesh.Edges, 2);
    fprintf('  - 生成全局棱数: %d\n', numEdges);
    
    % 5. 构建 T2E (Element to Edge Map)
    % ic 的维度是 6Ne x 1，重塑回 6 x Ne
    Mesh.T2E = reshape(ic, 6, numElems);
    
    % 6. 计算方向符号 T2E_Sign
    % 规则: 
    % 局部棱方向由 Mesh.T 中的节点顺序决定 (Node_local_a -> Node_local_b)
    % 全局棱方向由 sort 决定 (Node_min -> Node_max)
    % 如果 raw_edges(k,1) == raw_edges_sorted(k,1)，说明原方向就是从小到大，Sign = +1
    % 否则 Sign = -1
    
    fprintf('  - 计算棱方向符号...\n');
    is_same_dir = (raw_edges(:, 1) == raw_edges_sorted(:, 1));
    signs_vec = int8(is_same_dir) * 2 - 1; % true->1, false->-1
    
    Mesh.T2E_Sign = reshape(signs_vec, 6, numElems);
    
    % 7. 预计算单元体积 (用于诊断和积分权重)
    fprintf('  - 计算单元体积...\n');
    % 使用简单的行列式计算体积
    % Vol = det([1 1 1 1; x; y; z]) / 6
    P = Mesh.P;
    T = Mesh.T;
    
    % 向量化计算体积
    % V1 = P(:, T(1,:)); V2 = P(:, T(2,:)); ...
    % d1 = V2-V1; d2 = V3-V1; d3 = V4-V1;
    % vol = dot(cross(d1, d2), d3) / 6
    v1 = P(:, T(1,:));
    v2 = P(:, T(2,:));
    v3 = P(:, T(3,:));
    v4 = P(:, T(4,:));
    
    d1 = v2 - v1;
    d2 = v3 - v1;
    d3 = v4 - v1;
    
    % 叉乘 (手动展开更快)
    c1 = d1(2,:).*d2(3,:) - d1(3,:).*d2(2,:);
    c2 = d1(3,:).*d2(1,:) - d1(1,:).*d2(3,:);
    c3 = d1(1,:).*d2(2,:) - d1(2,:).*d2(1,:);
    
    Mesh.Volumes = (c1.*d3(1,:) + c2.*d3(2,:) + c3.*d3(3,:)) / 6.0;
    
    % 检查负体积 (意味着节点编号不符合右手定则，可能需要纠正 T)
    neg_vols = find(Mesh.Volumes < 0);
    if ~isempty(neg_vols)
        warning('检测到 %d 个单元体积为负，建议检查网格节点顺序。代码将自动取绝对值。', length(neg_vols));
        Mesh.Volumes = abs(Mesh.Volumes);
    end
    
    fprintf('网格拓扑构建完成。\n');
end