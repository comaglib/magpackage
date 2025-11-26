function solve_nonlinear_static()
% SOLVE_NONLINEAR_STATIC 完整流程演示：非线性静磁场求解 (修正版)

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   3D Nonlinear Magnetostatic Solver          \n');
    fprintf('==============================================\n');
    
    % 1. 前处理
    msh_file = 'example_nonlinear.msh';
    create_iron_air_mesh(msh_file); 
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 2. 材料定义
    % Material 1: Air
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    % Material 2: Nonlinear Steel
    Model.Materials.Lib(2).Name = 'Steel-1010';
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Sigma = 0;
    
    % [修复推荐] 显式设置初始 Mu_r，虽然 preprocess 也会自动补全
    % 根据 H=50, B=0.4 -> mu = 0.4/50 = 0.008 -> mu_r = 0.008/mu0 ~ 6366
    Model.Materials.Lib(2).Mu_r = 6000; 
    
    H_data = [0, 50, 100, 200, 300, 500, 1000, 2000, 5000, 10000, 50000];
    B_data = [0, 0.4, 0.7, 1.1, 1.3, 1.45, 1.55, 1.65, 1.8, 1.9, 2.1];
    
    Model.Materials.Lib(2).BH_Curve.H = H_data;
    Model.Materials.Lib(2).BH_Curve.B = B_data;
    
    % 预处理 (现在会自动填补 Mu_r，如果上面没写的话)
    Model = preprocess_materials(Model);
    
    % 映射材料
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 定义激励源
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.5; 
    CoilParams.Radius = 0.2; 
    CoilParams.Current = 50000; % 50 kA
    CoilParams.N_seg = 36;
    
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 组装右端项 (RHS)
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    fprintf('组装初始右端项...\n');
    b_ext = assemble_rhs_reduced(Model, Coil);
    
    % 5. 求解器配置
    Model.Solver.Nonlinear.Method = 'NewtonRaphson';
    Model.Solver.Nonlinear.MaxIter = 25;
    Model.Solver.Nonlinear.Tolerance = 1e-5;
    
    Model.Solver.Gauge.Method = 'Lagrange';
    
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true;
    Model.Solver.Linear.MumpsICNTL.i8 = 77; % Auto scaling
    
    % 边界条件
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    Model.Runtime.FixedEdges = fixed_edges;
    Model.Runtime.FixedVals = zeros(size(fixed_edges));
    
    % 6. 执行求解
    nEdges = size(Model.Mesh.Edges, 2);
    nNodes = size(Model.Mesh.P, 2);
    x0 = zeros(nEdges + nNodes, 1);
    
    [x_sol, info] = solve_nr(Model, b_ext, x0);
    
    Ar = x_sol(1:nEdges);
    
    % 7. 结果验证
    fprintf('\n--- 结果分析 ---\n');
    fprintf('迭代步数: %d\n', info.Iterations);
    fprintf('最终残差: %.4e\n', info.Residuals(end));
    
    % 验证饱和度
    iron_elems = find(Model.Mesh.RegionTags == 2);
    center_elem = iron_elems(1);
    
    P = Model.Mesh.P; T = Model.Mesh.T;
    nodes = T(:, center_elem);
    pts = P(:, nodes);
    center_pt = mean(pts, 2);
    
    B_s = compute_biot_savart_B(Coil, center_pt);
    
    edges = Model.Mesh.T2E(:, center_elem);
    signs = double(Model.Mesh.T2E_Sign(:, center_elem));
    a_vec = Ar(edges) .* signs;
    B_r = calc_element_B(pts, a_vec);
    
    B_total = B_r + B_s;
    B_mag = norm(B_total);
    
    fprintf('中心总磁场 B = %.4f T\n', B_mag);
    
    mat_info = Model.Materials.Lib(2);
    [nu_sat, ~] = eval_material_nu(B_mag^2, mat_info);
    mu0 = 4*pi*1e-7;
    mu_r_eff = 1 / (nu_sat * mu0);
    
    fprintf('有效相对磁导率 mu_r = %.2f (初始 ~%.0f)\n', mu_r_eff, Model.Materials.Lib(2).Mu_r);
    
    if mu_r_eff < Model.Materials.Lib(2).Mu_r * 0.8
        fprintf('[PASS] Material is saturating.\n');
    else
        fprintf('[INFO] Material is in linear region.\n');
    end
end

function create_iron_air_mesh(filename)
    % 生成简单的 Box-in-Box 网格
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n1 2 3 4 5\n');
    fprintf(fid, '0 0 0\n0.1 0 0\n0 0.1 0\n0 0 0.1\n0.2 0.2 0.2\n$EndNodes\n');
    fprintf(fid, '$Elements\n3 3 1 3\n');
    fprintf(fid, '2 100 2 1\n1 2 3 5\n'); % Surface
    fprintf(fid, '3 2 4 1\n2 1 2 3 4\n'); % Iron
    fprintf(fid, '3 1 4 1\n3 2 3 4 5\n'); % Air
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end