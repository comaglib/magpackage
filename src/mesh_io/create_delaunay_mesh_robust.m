function create_delaunay_mesh_robust(filename)
% CREATE_DELAUNAY_MESH_ROBUST 生成鲁棒的四面体网格文件 (.msh)
% 
% 功能:
%   生成包含铁心(Inner)和空气(Outer)的测试网格。
%   强制检查并修正四面体节点顺序，确保体积为正。
%   同时支持返回数据结构和写入文件。
% 
% 输入:
%   filename - 输出 .msh 文件名 (可选)
% 
% 输出:
%   MeshRaw  - 包含 .P, .T, .RegionTags, .Faces, .FaceTags 的结构体
% 
% 修复:
% 1. 使用 %.16g 高精度写入节点坐标，防止因 ASCII 截断导致的单元退化(负体积)。
% 2. 写入前执行严格的手性检查和修正。

    fprintf('  [MeshGen] Generating mesh via Delaunay (High Precision)...\n');
    
    % 1. 点云生成
    % 铁心 [-0.15, 0.15]
    [X1, Y1, Z1] = ndgrid(linspace(-0.15, 0.15, 5));
    P_inner = [X1(:), Y1(:), Z1(:)];
    
    % 空气 [-0.5, 0.5]
    [X2, Y2, Z2] = ndgrid(linspace(-0.5, 0.5, 6));
    P_outer = [X2(:), Y2(:), Z2(:)];
    
    P = unique([P_inner; P_outer], 'rows');
    nNodes = size(P, 1);
    
    % 2. Delaunay 剖分
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList; 
    
    % 3. 手性修正 (Chirality Fix)
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    % 混合积 (p2-p1) . ((p3-p1) x (p4-p1))
    vols = sum(cross(p2-p1, p3-p1, 2) .* (p4-p1), 2);
    
    neg_mask = vols < 0;
    if any(neg_mask)
        fprintf('    Fixing %d inverted elements...\n', sum(neg_mask));
        T(neg_mask, [1 2]) = T(neg_mask, [2 1]);
    end
    
    % 再次检查
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    vols_check = sum(cross(p2-p1, p3-p1, 2) .* (p4-p1), 2);
    if any(vols_check < 0)
        error('Fatal: Failed to fix mesh chirality.');
    end
    
    % 4. 准备数据结构 (转置为 MATLAB 列优先格式)
    P_io = P'; 
    T_io = T'; 
    nElems = size(T_io, 2);
    
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    in_iron = abs(Centers(:,1)) < 0.16 & abs(Centers(:,2)) < 0.16 & abs(Centers(:,3)) < 0.16;
    
    ElemTags = ones(1, nElems); 
    ElemTags(in_iron) = 2;
    
    TR = triangulation(T, P);
    FB = freeBoundary(TR); 
    BoundaryFaces = FB'; 
    nFaces = size(BoundaryFaces, 2);
    
    % 5. 写入文件
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    
    % Nodes
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    % 构造矩阵: [ID; X; Y; Z]
    NodeData = [(1:nNodes); P_io];
    % [关键修复] 使用 %.16g 保证精度
    fprintf(fid, '%d %.16g %.16g %.16g\n', NodeData);
    fprintf(fid, '$EndNodes\n');
    
    % Elements
    nTotal = nFaces + nElems;
    idx_air = (ElemTags == 1);
    idx_iron = (ElemTags == 2);
    nAir = sum(idx_air);
    nIron = sum(idx_iron);
    
    nBlocks = 1 + (nAir>0) + (nIron>0);
    
    fprintf(fid, '$Elements\n%d %d 1 %d\n', nBlocks, nTotal, nTotal);
    
    % Surface
    fprintf(fid, '2 100 2 %d\n', nFaces);
    % Tri Format: ID n1 n2 n3
    TriData = [(1:nFaces); BoundaryFaces];
    fprintf(fid, '%d %d %d %d\n', TriData);
    
    id = nFaces + 1;
    % Air
    if nAir > 0
        fprintf(fid, '3 1 4 %d\n', nAir);
        % Tet Format: ID n1 n2 n3 n4
        AirData = [(id : id+nAir-1); T_io(:, idx_air)];
        fprintf(fid, '%d %d %d %d %d\n', AirData);
        id = id + nAir;
    end
    
    % Iron
    if nIron > 0
        fprintf(fid, '3 2 4 %d\n', nIron);
        IronData = [(id : id+nIron-1); T_io(:, idx_iron)];
        fprintf(fid, '%d %d %d %d %d\n', IronData);
    end
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end