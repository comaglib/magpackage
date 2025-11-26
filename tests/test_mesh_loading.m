function test_mesh_loading()
% TEST_MESH_LOADING 验证网格读取与拓扑构建模块

    addpath(genpath('src'));
    
    fprintf('Starting Mesh Loading Test...\n');
    
    % 1. 创建临时文件
    msh_file = 'test_cube_v2.msh';
    create_dummy_msh(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    % 2. 读取
    Raw = read_msh(msh_file);
    
    % 调试输出
    if isempty(Raw.T)
        error('Failed to read elements from dummy file.');
    end
    
    % 3. 构建拓扑
    Mesh = build_topology(Raw);
    
    nV = size(Mesh.P, 2);
    nE = size(Mesh.Edges, 2);
    nT = size(Mesh.T, 2);
    
    fprintf('Nodes: %d, Edges: %d, Elements: %d\n', nV, nE, nT);
    
    % 检查: 一个标准四面体应该有 6 条棱
    assert(nT == 1, 'Should have 1 element');
    assert(nV == 4, 'Should have 4 nodes');
    if nE ~= 6
        disp('Mesh.T:');
        disp(Mesh.T);
        error('Edge count mismatch! Expected 6 for a single tet, got %d. Check read_msh.', nE);
    end
    
    % 4. 验证 T2E 符号
    elem_idx = 1;
    % 检查第1条局部棱: Node1 -> Node2
    edge_local_idx = 1; 
    
    n1 = Mesh.T(1, elem_idx);
    n2 = Mesh.T(2, elem_idx);
    
    sign_recorded = Mesh.T2E_Sign(edge_local_idx, elem_idx);
    
    % 逻辑判断
    if n1 < n2
        expected_sign = 1; % Global def is n1->n2 (Same)
    else
        expected_sign = -1; % Global def is n2->n1 (Opposite)
    end
    
    if sign_recorded == expected_sign
        fprintf('[PASS] Edge Orientation Sign Check (Sign=%d)\n', sign_recorded);
    else
        error('[FAIL] Sign Mismatch! n1=%d, n2=%d, Rec=%d, Exp=%d', ...
            n1, n2, sign_recorded, expected_sign);
    end
    
    fprintf('Mesh Test Passed!\n');
end

function create_dummy_msh(filename)
    % 创建一个有效的 Gmsh 4.1 ASCII 文件
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    
    % Nodes Block
    fprintf(fid, '$Nodes\n');
    fprintf(fid, '1 4 1 4\n'); % 1 block, 4 nodes total, min 1, max 4
    fprintf(fid, '2 1 0 4\n'); % Dim=2, Tag=1, Param=0, Nodes=4
    fprintf(fid, '1 2 3 4\n'); % Node Tags
    fprintf(fid, '0 0 0\n1 0 0\n0 1 0\n0 0 1\n'); % Coordinates
    fprintf(fid, '$EndNodes\n');
    
    % Elements Block
    fprintf(fid, '$Elements\n');
    fprintf(fid, '1 1 1 1\n'); % 1 block, 1 elem total, min 1, max 1
    fprintf(fid, '2 1 4 1\n'); % Dim=2, Tag=1, Type=4(Tet), Count=1
    % Element Data: ID Node1 Node2 Node3 Node4
    fprintf(fid, '1 1 2 3 4\n'); 
    fprintf(fid, '$EndElements\n');
    
    fclose(fid);
end