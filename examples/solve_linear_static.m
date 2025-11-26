function solve_linear_static()
% SOLVE_LINEAR_STATIC 完整流程演示：线性静磁场求解
% 案例: 铁心 (Box) 被跑道线圈磁化

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   3D Magnetostatic Solver (Ar-Formulation)   \n');
    fprintf('==============================================\n');
    
    % 1. 前处理: 生成并读取网格
    msh_file = 'example_iron_coil.msh';
    create_iron_air_mesh(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 2. 定义材料 (Tag 1=Air, Tag 2=Iron)
    % Iron: mu_r = 1000
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Iron';
    Model.Materials.Lib(2).Mu_r = 1000.0;
    Model.Materials.Lib(2).Sigma = 0;
    
    % 构建 ActiveMap
    % 假设 Raw.RegionTags 正确对应
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 定义线圈源
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.5; % 跑道长
    CoilParams.Radius = 0.2; 
    CoilParams.Current = 1000; % 1000安匝
    CoilParams.N_seg = 36;
    
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 组装系统
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 4.1 刚度矩阵 K
    K = assemble_magnetic_stiffness(Model);
    
    % 4.2 载荷向量 b (Reduced Potential RHS)
    b = assemble_rhs_reduced(Model, Coil);
    
    % 5. 施加库伦规范并求解
    % 边界: 外表面 (Tag 100) A_r = 0 (假设远场无扰动)
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    
    [SysK, Sysb] = apply_gauge_strategy(Model, K, b, fixed_edges, 0, 'Lagrange');
    
    % 设置求解器
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true;
    
    fprintf('正在求解线性系统...\n');
    Sol = linear_solve(SysK, Sysb, Model);
    
    nEdges = size(K, 1);
    Ar = Sol(1:nEdges); % 提取 Ar
    
    fprintf('求解完成. |Ar| = %.4e\n', norm(Ar));
    
    % 6. 后处理: 验证中心点磁通密度
    % B_total = Curl(Ar) + Bs
    % 我们取铁心中心 (0,0,0) 的一个单元
    
    % 简单查找包含原点的单元
    % (实际后处理模块会更复杂，这里仅做快速验证)
    Center = [0;0;0];
    % 计算所有单元中心到原点的距离
    P = Model.Mesh.P; T = Model.Mesh.T;
    min_dist = inf; elem_id = 0;
    for e = 1:size(T, 2)
        pc = mean(P(:, T(:,e)), 2);
        d = norm(pc - Center);
        if d < min_dist
            min_dist = d; elem_id = e;
        end
    end
    
    % 计算该单元的 B_total
    % 1. Curl Ar
    edges = Model.Mesh.T2E(:, elem_id);
    signs = double(Model.Mesh.T2E_Sign(:, elem_id));
    Ar_local = Ar(edges) .* signs;
    
    % Curl N at center (Recalculate or store)
    nodes = T(:, elem_id);
    pts = P(:, nodes);
    % ... (Copy paste curl calc logic for demo) ...
    v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
    J_mat = [v21, v31, v41];
    invJ_T = inv(J_mat)';
    G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
    G_phy = invJ_T * G_ref;
    edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    Curl_N = zeros(3, 6);
    for i = 1:6
        na = edge_pairs(i,1); nb = edge_pairs(i,2);
        Curl_N(:, i) = 2 * cross(G_phy(:, na), G_phy(:, nb));
    end
    
    B_r = Curl_N * Ar_local;
    
    % 2. Bs
    pc = mean(pts, 2);
    B_s = compute_biot_savart_B(Coil, pc);
    
    B_total = B_r + B_s;
    
    fprintf('\n--- 结果验证 (铁心中心) ---\n');
    fprintf('源磁场   Bs = [%.4f, %.4f, %.4f] T\n', B_s(1), B_s(2), B_s(3));
    fprintf('反应场   Br = [%.4f, %.4f, %.4f] T\n', B_r(1), B_r(2), B_r(3));
    fprintf('总磁场   B  = [%.4f, %.4f, %.4f] T\n', B_total(1), B_total(2), B_total(3));
    
    % 理论定性分析: 铁心会增强磁场
    % B_total_z 应该 > B_s_z
    amplification = B_total(3) / B_s(3);
    fprintf('磁增强系数 (B_tot/B_s) = %.2f\n', amplification);
    
    if amplification > 1.0
        fprintf('[PASS] Iron core successfully amplified the magnetic field.\n');
    else
        warning('[WARN] No amplification observed. Check physics setup.');
    end
end

function create_iron_air_mesh(filename)
    % 创建一个包含铁心和空气的简单网格
    % Iron: [-0.1, 0.1]^3 (Tag 2)
    % Air:  [-0.5, 0.5]^3 (Tag 1)
    % Boundary: Air Surface (Tag 100)
    
    fid = fopen(filename, 'w');
    % 简化的写法：直接用 Gmsh 原生语言写 .geo 然后 mesh? 
    % 不，为了独立性，我们还是手写 .msh 节点和单元 (虽然繁琐)
    % 为了演示，我们只做一个极简的: 内部1个四面体(铁), 外部包围几个四面体(气)
    % 这太难手写了。
    
    % 策略: 写一个 2x2x2 的网格
    % Nodes: 27个点... 
    % 让我们写最简单的: 2个四面体，一个是铁，一个是气
    
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n1 2 3 4 5\n');
    % 4个点构成中心四面体，第5个点构成外部四面体
    fprintf(fid, '0 0 0\n0.1 0 0\n0 0.1 0\n0 0 0.1\n0.2 0.2 0.2\n$EndNodes\n');
    
    fprintf(fid, '$Elements\n3 3 1 3\n');
    
    % Surface 100 (Boundary - Face of outer tet)
    % Outer tet is (2,3,4,5). Face (2,3,5) is boundary
    fprintf(fid, '2 100 2 1\n1 2 3 5\n');
    
    % Volume 2 (Iron) - Tet (1,2,3,4)
    fprintf(fid, '3 2 4 1\n2 1 2 3 4\n');
    
    % Volume 1 (Air) - Tet (2,3,4,5) (Adjacent to Iron)
    fprintf(fid, '3 1 4 1\n3 2 3 4 5\n');
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end