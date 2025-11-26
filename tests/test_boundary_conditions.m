function test_boundary_conditions()
% TEST_BOUNDARY_CONDITIONS 验证 Dirichlet 边界条件的施加与求解

    addpath(genpath('src'));
    fprintf('Starting Boundary Condition Test...\n');
    
    % 1. 生成带 Surface Tag 的测试网格
    msh_file = 'test_bc.msh';
    create_bc_test_msh(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    % 2. 读取与构建
    Raw = read_msh(msh_file);
    % 验证是否读到了 Faces
    assert(~isempty(Raw.Faces), 'Failed to read surface triangles.');
    
    Model.Mesh = build_topology(Raw);
    
    % 3. 组装空载荷系统 (Kx = 0)
    % 材料
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    Model.Materials.ActiveMap = ones(1, size(Model.Mesh.T, 2));
    
    % 启动 Pool
    if isempty(gcp('nocreate')), parpool('local'); end
    
    K = assemble_magnetic_stiffness(Model);
    nEdges = size(K, 1);
    b = zeros(nEdges, 1);
    
    % 4. 识别边界
    % 在 create_bc_test_msh 中，我们定义 Surface Tag 100 为底面
    bc_tag = 100;
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, bc_tag);
    
    assert(~isempty(fixed_edges), 'Failed to identify boundary edges.');
    
    % 5. 施加非零边界条件 (比如 A = 5.0)
    bc_val = 5.0;
    [K_bc, b_bc] = apply_dirichlet_bc(K, b, fixed_edges, bc_val);
    
    % 6. 求解
    Model.Solver.Linear.Interface = 'MUMPS'; % 此时应已配置好
    Model.Solver.Linear.Symmetric = true;
    x = linear_solve(K_bc, b_bc, Model);
    
    % 7. 验证
    % 检查边界棱的值是否等于 5.0
    x_boundary = x(fixed_edges);
    err = norm(x_boundary - bc_val, inf);
    fprintf('Boundary Value Error: %e\n', err);
    
    if err < 1e-10
        fprintf('[PASS] Dirichlet BC enforced correctly.\n');
    else
        error('[FAIL] Boundary values do not match target.');
    end
    
    % 检查内部值 (对于拉普拉斯方程，内部应该是边界的插值，此处全域常数 5.0 也是解)
    % 因为 curl(curl A) = 0, A=const 导致 B=0, curl B = 0. 满足方程.
    % 所以全域应该都是 5.0 (如果连通)
    err_internal = norm(x - bc_val, inf);
    fprintf('Internal Field Uniformity Error: %e\n', err_internal);
    if err_internal < 1e-10
        fprintf('[PASS] Trivial solution (Constant field) verified.\n');
    end
end

function create_bc_test_msh(filename)
    % 创建一个单个四面体，并定义其中一个面为物理表面 100
    % Nodes: 1(0,0,0), 2(1,0,0), 3(0,1,0), 4(0,0,1)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    
    fprintf(fid, '$Nodes\n1 4 1 4\n2 1 0 4\n1 2 3 4\n');
    fprintf(fid, '0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n');
    
    fprintf(fid, '$Elements\n2 2 1 2\n'); 
    % Block 1: Surface 100 (Type 2), 1 elem: Face (1,2,3) (z=0平面)
    fprintf(fid, '2 100 2 1\n'); 
    fprintf(fid, '1 1 2 3\n');
    
    % Block 2: Volume 1 (Type 4), 1 elem: Tet (1,2,3,4)
    fprintf(fid, '3 1 4 1\n');
    fprintf(fid, '2 1 2 3 4\n');
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end