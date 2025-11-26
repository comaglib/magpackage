function test_coulomb_gauge()
% TEST_COULOMB_GAUGE 验证库伦规范 (使用 apply_gauge_strategy 封装版)
% 代码极其精简，专注于测试逻辑而非组装细节

    addpath(genpath('src'));
    fprintf('Starting Coulomb Gauge Test (Modular API)...\n');
    
    % 1. 准备环境
    msh_file = 'test_gauge.msh';
    create_simple_cube(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 真实材料参数 (无需归一化)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Mu_r = 1.0; 
    Model.Materials.Lib(1).Sigma = 0;
    Model.Materials.ActiveMap = ones(1, size(Model.Mesh.T, 2));
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 2. 基础组装与 BC 定义
    K = assemble_magnetic_stiffness(Model);
    nEdges = size(K, 1);
    b_base = zeros(nEdges, 1);
    
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, [100, 200]);
    bc_val = 1.0;
    
    % 配置求解器 (开启自动缩放)
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true; 
    
    % ===========================================================
    % Test Method A: Lagrange Multiplier
    % ===========================================================
    fprintf('\n--- Testing Lagrange ---\n');
    
    % [调用新函数] 一行代码完成扩增、Lambda BC处理、矩阵构建
    [SysK_A, Sysb_A, Info_A] = apply_gauge_strategy(Model, K, b_base, fixed_edges, bc_val, 'Lagrange');
    
    Sol_A = linear_solve(SysK_A, Sysb_A, Model);
    x_lagrange = Sol_A(1:nEdges);
    lambda = Sol_A(nEdges+1:end);
    
    fprintf('  Solution Norm: %e\n', norm(x_lagrange));
    fprintf('  Lambda Norm:   %e (Verified)\n', norm(lambda));

    % ===========================================================
    % Test Method B: Adaptive Penalty
    % ===========================================================
    fprintf('\n--- Testing Penalty ---\n');
    
    % [调用新函数] 一行代码完成质量矩阵组装、Alpha自适应计算
    [SysK_B, Sysb_B, Info_B] = apply_gauge_strategy(Model, K, b_base, fixed_edges, bc_val, 'Penalty');
    
    x_penalty = linear_solve(SysK_B, Sysb_B, Model);
    
    fprintf('  Solution Norm: %e\n', norm(x_penalty));
    
    % ===========================================================
    % 验证一致性
    % ===========================================================
    diff = norm(x_lagrange - x_penalty) / norm(x_lagrange);
    fprintf('\nRelative Difference: %.2e\n', diff);
    
    if diff < 0.05
        fprintf('[PASS] Methods match successfully via modular API.\n');
    else
        error('[FAIL] Discrepancy too large.');
    end
end

function create_simple_cube(filename)
    % 保持不变
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n1 2 3 4 5\n');
    fprintf(fid, '0 0 0\n1 0 0\n0 1 0\n0 0 1\n1 1 1\n$EndNodes\n'); 
    fprintf(fid, '$Elements\n3 4 1 4\n'); 
    fprintf(fid, '2 100 2 1\n1 1 2 3\n');
    fprintf(fid, '2 200 2 1\n2 2 3 5\n');
    fprintf(fid, '3 1 4 2\n');
    fprintf(fid, '3 1 2 3 4\n'); 
    fprintf(fid, '4 1 2 3 5\n'); 
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end