function MeshRaw = create_test_mesh_direct(filename)
    fprintf('  [MeshGen] Generating mesh via Delaunay (Direct Memory)...\n');
    
    % 1. 点云生成
    % 铁心区域 [-0.15, 0.15]
    [X1, Y1, Z1] = ndgrid(linspace(-0.15, 0.15, 5));
    P_inner = [X1(:), Y1(:), Z1(:)];
    
    % 空气区域 [-0.6, 0.6] (确保完全包裹线圈)
    [X2, Y2, Z2] = ndgrid(linspace(-0.6, 0.6, 7));
    P_outer = [X2(:), Y2(:), Z2(:)];
    
    P = unique([P_inner; P_outer], 'rows');
    nNodes = size(P, 1);
    
    % 2. Delaunay 三角剖分
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList; 
    
    % 3. 手性修正 (Chirality Fix)
    % 确保所有四面体体积为正 (符合右手定则)
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    vols = sum(cross(p2-p1, p3-p1, 2) .* (p4-p1), 2);
    
    neg_mask = vols < 0;
    if any(neg_mask)
        % 交换节点 1 和 2 以翻转手性
        T(neg_mask, [1 2]) = T(neg_mask, [2 1]);
    end
    
    % 4. 构建数据结构
    P_io = P'; 
    T_io = T'; 
    nElems = size(T_io, 2);
    
    % 材料区域标记
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    in_iron = abs(Centers(:,1)) < 0.16 & abs(Centers(:,2)) < 0.16 & abs(Centers(:,3)) < 0.16;
    
    ElemTags = ones(1, nElems); 
    ElemTags(in_iron) = 2; % Tag 2 = Iron
    
    % 提取边界表面
    TR = triangulation(T, P);
    FB = freeBoundary(TR); 
    BoundaryFaces = FB'; 
    nFaces = size(BoundaryFaces, 2);
    FaceTags = repmat(100, 1, nFaces); % Tag 100 = Boundary
    
    % 5. 返回 Mesh 结构体
    MeshRaw.P = P_io;
    MeshRaw.T = T_io;
    MeshRaw.RegionTags = ElemTags;
    MeshRaw.Faces = BoundaryFaces;
    MeshRaw.FaceTags = FaceTags;
    
    % 6. 写入文件 (用于后续检查或可视化)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    % 格式: ID X Y Z
    NodeData = [(1:nNodes); P_io];
    fprintf(fid, '%d %.6f %.6f %.6f\n', NodeData);
    fprintf(fid, '$EndNodes\n');
    
    nTotal = nFaces + nElems;
    idx_air = (ElemTags == 1);
    idx_iron = (ElemTags == 2);
    nAir = sum(idx_air);
    nIron = sum(idx_iron);
    
    nBlocks = 1 + (nAir>0) + (nIron>0);
    fprintf(fid, '$Elements\n%d %d 1 %d\n', nBlocks, nTotal, nTotal);
    
    % Surface Block
    fprintf(fid, '2 100 2 %d\n', nFaces);
    TriData = [(1:nFaces); BoundaryFaces];
    fprintf(fid, '%d %d %d %d\n', TriData);
    
    id = nFaces + 1;
    % Air Block
    if nAir > 0
        fprintf(fid, '3 1 4 %d\n', nAir);
        AirData = [(id : id+nAir-1); T_io(:, idx_air)];
        fprintf(fid, '%d %d %d %d %d\n', AirData);
        id = id + nAir;
    end
    
    % Iron Block
    if nIron > 0
        fprintf(fid, '3 2 4 %d\n', nIron);
        IronData = [(id : id+nIron-1); T_io(:, idx_iron)];
        fprintf(fid, '%d %d %d %d %d\n', IronData);
    end
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end