function outer_tags = detect_outer_boundary_tags(MeshRaw)
% DETECT_OUTER_BOUNDARY_TAGS 自动识别网格最外层边界的 Face Tags
%
% 原理: 计算网格包围盒，采样每个 Tag 对应的面，检查节点是否位于包围盒表面。
%
% 存放位置: src/physics_setup/detect_outer_boundary_tags.m

    if ~isfield(MeshRaw, 'Faces') || isempty(MeshRaw.Faces)
        warning('detect_outer_boundary_tags:NoFaces', '网格数据中不包含 Faces 信息，无法检测边界。');
        outer_tags = [];
        return;
    end

    P = MeshRaw.P;
    
    % 1. 计算全局包围盒
    min_box = min(P, [], 2);
    max_box = max(P, [], 2);
    
    % 设定容差 (相对于网格尺寸)
    L_char = max(max_box - min_box);
    tol = L_char * 1e-3; 
    
    unique_tags = unique(MeshRaw.FaceTags);
    outer_tags = [];
    
    for t = unique_tags
        % 获取该 Tag 的所有面索引
        f_ids = find(MeshRaw.FaceTags == t);
        
        % 采样优化: 为了速度，只检查前 50 个面
        % 如果这 50 个面都在边界上，我们认为整个 Tag 都是边界
        n_sample = min(length(f_ids), 50);
        sample_faces = MeshRaw.Faces(:, f_ids(1:n_sample));
        
        % 获取采样面的所有节点
        pts = P(:, sample_faces(:));
        
        % 检查节点是否落在 6 个包围盒平面上
        % (X=Xmin OR X=Xmax OR Y=Ymin ...)
        is_on_bd = ...
            abs(pts(1,:) - min_box(1)) < tol | abs(pts(1,:) - max_box(1)) < tol | ...
            abs(pts(2,:) - min_box(2)) < tol | abs(pts(2,:) - max_box(2)) < tol | ...
            abs(pts(3,:) - min_box(3)) < tol | abs(pts(3,:) - max_box(3)) < tol;
        
        % 如果所有采样节点都在边界上，则判定为外边界 Tag
        if all(is_on_bd)
            outer_tags = [outer_tags, t]; %#ok<AGROW>
        end
    end
end