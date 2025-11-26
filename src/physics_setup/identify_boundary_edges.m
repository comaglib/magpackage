function fixed_edge_indices = identify_boundary_edges(Mesh, MeshRaw, target_surface_tags)
% IDENTIFY_BOUNDARY_EDGES 查找属于指定物理表面的全局棱索引
% 
% 输入:
%   Mesh       - 已构建拓扑的网格 (必须包含 Mesh.Edges)
%   MeshRaw    - 原始网格数据 (包含 Faces 和 FaceTags)
%   target_tags- 需要查找的物理表面 Tag 列表 (数组)
%
% 输出:
%   fixed_edge_indices - 符合条件的全局棱索引列向量

    if isempty(MeshRaw.Faces)
        warning('MeshRaw 不包含表面单元信息。无法识别边界条件。');
        fixed_edge_indices = [];
        return;
    end

    % 1. 筛选目标表面单元
    % 找出 FaceTags 属于 target_surface_tags 的列索引
    mask = ismember(MeshRaw.FaceTags, target_surface_tags);
    target_faces = MeshRaw.Faces(:, mask); % 3 x N_faces
    
    if isempty(target_faces)
        warning('未找到指定 Tag (%s) 的表面单元。', num2str(target_surface_tags));
        fixed_edge_indices = [];
        return;
    end
    
    % 2. 提取这些表面单元的所有棱
    % 一个三角形有 3 条棱: (1-2), (2-3), (3-1)
    nFaces = size(target_faces, 2);
    
    % 构造棱列表 (起点, 终点)
    e1 = target_faces([1, 2], :);
    e2 = target_faces([2, 3], :);
    e3 = target_faces([3, 1], :);
    
    boundary_edges_raw = [e1, e2, e3]; % 2 x (3*N_faces)
    
    % 3. 规范化方向 (min -> max) 以便匹配
    boundary_edges_sorted = sort(boundary_edges_raw, 1)'; % 转置为 N x 2
    
    % 去重 (因为相邻三角形共用棱)
    boundary_edges_unique = unique(boundary_edges_sorted, 'rows');
    
    % 4. 在全局棱列表中查找索引
    % Mesh.Edges 也是 (2 x N_global)，先转置
    global_edges = Mesh.Edges'; % N_global x 2
    
    % 使用 ismember 查找
    % [Lia, Locb] = ismember(A, B, 'rows')
    % Lia 是逻辑索引，表示 A 中的行是否存在于 B 中
    % Locb 是 B 中的索引 (我们需要这个，但 ismember 的 Locb 对于重复项行为需注意)
    % 实际上我们只需要知道 B 中的哪些行对应 A。
    % 更快的方法: ismember(global, boundary) -> logic index
    
    [is_boundary, ~] = ismember(global_edges, boundary_edges_unique, 'rows');
    
    fixed_edge_indices = find(is_boundary);
    
    fprintf('  - 识别到 %d 条边界棱 (Surface Tags: %s)\n', ...
        length(fixed_edge_indices), num2str(target_surface_tags));
end