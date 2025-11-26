function solve_linear_static()
% SOLVE_LINEAR_STATIC 完整流程演示：线性静磁场求解 (修正版)
% 
% 修复:
% 1. 更新 create_iron_air_mesh 以匹配最新的 read_msh 格式 (ID X Y Z)。
% 2. 解决了因网格读取错位导致的矩阵奇异和 NaN 问题。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   3D Magnetostatic Solver (Ar-Formulation)   \n');
    fprintf('==============================================\n');
    
    % 1. 前处理
    msh_file = 'example_iron_coil.msh';
    create_iron_air_mesh(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 2. 材料
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Iron';
    Model.Materials.Lib(2).Mu_r = 1000.0;
    Model.Materials.Lib(2).Sigma = 0;
    
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.5; 
    CoilParams.Radius = 0.2; 
    CoilParams.Current = 1000; 
    CoilParams.N_seg = 36;
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 组装
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    K = assemble_magnetic_stiffness(Model);
    b = assemble_rhs_reduced(Model, Coil);
    
    % 5. 求解
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    if isempty(fixed_edges), warning('No boundary edges found.'); end
    
    % 自动处理 Lagrange 规范 (包含 Lambda BC)
    [SysK, Sysb] = apply_gauge_strategy(Model, K, b, fixed_edges, 0, 'Lagrange');
    
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true;
    Model.Solver.Linear.MumpsICNTL.i8 = 77; % Auto scaling
    
    fprintf('正在求解线性系统...\n');
    Sol = linear_solve(SysK, Sysb, Model);
    
    nEdges = size(K, 1);
    Ar = Sol(1:nEdges); 
    
    fprintf('求解完成. |Ar| = %.4e\n', norm(Ar));
    
    % 6. 后处理
    iron_elems = find(Model.Mesh.RegionTags == 2);
    center_elem = iron_elems(1);
    
    P = Model.Mesh.P; T = Model.Mesh.T;
    nodes = T(:, center_elem);
    pts = P(:, nodes);
    center_pt = mean(pts, 2);
    
    edges = Model.Mesh.T2E(:, center_elem);
    signs = double(Model.Mesh.T2E_Sign(:, center_elem));
    a_vec = Ar(edges) .* signs;
    
    B_r = calc_element_B(pts, a_vec);
    B_s = compute_biot_savart_B(Coil, center_pt);
    B_total = B_r + B_s;
    
    fprintf('\n--- 结果验证 (铁心中心) ---\n');
    fprintf('源磁场   Bs = [%.4f, %.4f, %.4f] T\n', B_s(1), B_s(2), B_s(3));
    fprintf('反应场   Br = [%.4f, %.4f, %.4f] T\n', B_r(1), B_r(2), B_r(3));
    fprintf('总磁场   B  = [%.4f, %.4f, %.4f] T\n', B_total(1), B_total(2), B_total(3));
    
    amplification = B_total(3) / B_s(3);
    fprintf('磁增强系数 (B_tot/B_s) = %.2f\n', amplification);
    
    if amplification > 1.0
        fprintf('[PASS] Iron core successfully amplified the magnetic field.\n');
    else
        warning('[WARN] No amplification observed. Check physics setup.');
    end
end

function create_iron_air_mesh(filename)
    % [修正] 生成格式兼容 read_msh 的网格 (ID X Y Z)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n');
    
    % 节点列表: ID X Y Z
    fprintf(fid, '1 0 0 0\n');
    fprintf(fid, '2 0.1 0 0\n');
    fprintf(fid, '3 0 0.1 0\n');
    fprintf(fid, '4 0 0 0.1\n');
    fprintf(fid, '5 0.2 0.2 0.2\n');
    
    fprintf(fid, '$EndNodes\n');
    
    fprintf(fid, '$Elements\n3 3 1 3\n');
    
    % Surface 100 (Boundary)
    fprintf(fid, '2 100 2 1\n1 2 3 5\n');
    
    % Volume 2 (Iron)
    fprintf(fid, '3 2 4 1\n2 1 2 3 4\n');
    
    % Volume 1 (Air)
    fprintf(fid, '3 1 4 1\n3 2 3 4 5\n');
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end