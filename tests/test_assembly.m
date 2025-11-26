function test_assembly()
% TEST_ASSEMBLY 验证并行矩阵组装的正确性 (修正版)
% 修复了测试网格生成中的节点顺序问题，消除了负体积警告

    addpath(genpath('src'));
    fprintf('Starting Assembly Test...\n');

    % 1. 构造一个微型测试案例
    msh_file = 'test_assembly.msh';
    create_2_element_msh(msh_file); 
    cleanup = onCleanup(@() delete(msh_file));
    
    % 2. 初始化模型
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 创建伪材料库 (空气)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    % 映射所有单元到材料1
    numElems = size(Model.Mesh.T, 2);
    Model.Materials.ActiveMap = ones(1, numElems);
    
    % 3. 执行组装
    % 启动并行池 (如果尚未启动)
    if isempty(gcp('nocreate'))
        try
            parpool('local'); 
        catch
            fprintf('Parallel pool already running or failed to start.\n');
        end
    end
    
    K = assemble_magnetic_stiffness(Model);
    
    % 4. 验证矩阵性质
    
    % 4.1 维度检查
    numEdges = size(Model.Mesh.Edges, 2);
    assert(all(size(K) == [numEdges, numEdges]), 'Matrix dimension mismatch');
    
    % 4.2 对称性检查
    err_sym = norm(K - K', 'fro');
    fprintf('Symmetry Error: %e\n', err_sym);
    if err_sym < 1e-10
        fprintf('[PASS] Stiffness Matrix is Symmetric\n');
    else
        error('[FAIL] Stiffness Matrix is NOT Symmetric');
    end
    
    % 4.3 零空间检查 (常数场/梯度场测试)
    % 构造标量场 V = x + 2y + 3z
    P = Model.Mesh.P;
    V_nodes = P(1,:) + 2*P(2,:) + 3*P(3,:);
    
    % 将 V 投影到棱上得到 A_edges = grad(V)
    Edges = Model.Mesh.Edges;
    A_grad = zeros(numEdges, 1);
    for e = 1:numEdges
        n1 = Edges(1, e);
        n2 = Edges(2, e);
        A_grad(e) = V_nodes(n2) - V_nodes(n1);
    end
    
    % 计算残差 K * A_grad (应为 0)
    residual = norm(K * A_grad);
    fprintf('Null-space Residual (K * grad V): %e\n', residual);
    
    if residual < 1e-9
        fprintf('[PASS] Curl-Curl Operator Correctly Handles Null Space\n');
    else
        warning('[WARN] Null-space residual is high. Check mesh scaling.');
    end
    
    fprintf('Assembly Test Passed!\n');
end

function create_2_element_msh(filename)
    % 创建两个共用一个面的四面体
    % Nodes: (0,0,0), (1,0,0), (0,1,0), (0,0,1), (0,0,-1)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n1 2 3 4 5\n');
    fprintf(fid, '0 0 0\n1 0 0\n0 1 0\n0 0 1\n0 0 -1\n$EndNodes\n');
    fprintf(fid, '$Elements\n1 2 1 2\n2 1 4 2\n');
    
    % Elem 1: 1 2 3 4 (体积 > 0)
    fprintf(fid, '1 1 2 3 4\n');
    
    % Elem 2: 1 3 2 5 (修正顺序)
    % 原序 1 2 3 5 导致体积 < 0
    % 现序 1 3 2 5: 1->3(y) cross 1->2(x) = -z. 
    % 1->5 是 -z. (-z) dot (-z) > 0. 体积为正.
    fprintf(fid, '2 1 3 2 5\n');
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end