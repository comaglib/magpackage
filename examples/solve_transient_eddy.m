function solve_transient_eddy()
% SOLVE_TRANSIENT_EDDY 测试瞬态涡流求解器 (终极修正版)
% 
% 修复:
% 1. 网格生成器现在包含严格的手性验证 (Double Check)，确保 0 负体积。
% 2. 物理参数调整以确保稳定性。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Transient Eddy Current Test (Step Response)\n');
    fprintf('==============================================\n');
    
    % 1. 生成网格 (严格保证正体积)
    msh_file = 'test_eddy_delaunay.msh';
    create_delaunay_mesh_robust(msh_file);
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    
    % 再次验证 (确保读取后也是正的)
    Model.Mesh = build_topology(Raw);
    if any(Model.Mesh.Volumes < 0)
        error('Fatal: Negative volumes detected after read. Mesh I/O pipeline is broken.');
    end
    
    % 2. 材料定义
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'TestConductor';
    Model.Materials.Lib(2).Type = 'Linear';
    Model.Materials.Lib(2).Mu_r = 1.0;
    Model.Materials.Lib(2).Sigma = 1e6; 
    
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.4; 
    CoilParams.Radius = 0.4; 
    CoilParams.Current = 0; 
    CoilParams.N_seg = 72;
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 时间设置
    TargetCurrent = 1e5; 
    TimeParams.I_func = @(t) TargetCurrent * (t > 0);
    TimeParams.dt = 0.005; 
    TimeParams.N_steps = 15;
    
    % 5. 边界条件
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    if isempty(fixed_edges), warning('Boundary failed.'); end
    Model.Runtime.FixedEdges = fixed_edges;
    
    % 6. 求解
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    [Sol, ~] = solve_transient_linear(Model, Coil, TimeParams);
    
    % 7. 后处理
    fprintf('\n--- 后处理验证 ---\n');
    Coil.I(:) = TargetCurrent;
    
    P = Model.Mesh.P; T = Model.Mesh.T;
    
    function B_mag = get_B_mag(pt)
        centers = squeeze(mean(reshape(P(:, T), 3, 4, []), 2));
        [~, e_id] = min(vecnorm(centers - pt, 2, 1));
        
        edges = Model.Mesh.T2E(:, e_id);
        signs = double(Model.Mesh.T2E_Sign(:, e_id));
        a_vec = Sol.A(edges) .* signs;
        
        nodes = T(:, e_id);
        elem_pts = P(:, nodes);
        B_r = calc_element_B(elem_pts, a_vec);
        
        B_s = compute_biot_savart_B(Coil, centers(:, e_id));
        B_mag = norm(B_r + B_s);
    end

    B_center = get_B_mag([0;0;0]);
    B_surface = get_B_mag([0.3; 0; 0]); 
    
    fprintf('中心磁场 B_center  = %.4e T\n', B_center);
    fprintf('表面磁场 B_surface = %.4e T\n', B_surface);
    
    ratio = B_center / B_surface;
    fprintf('屏蔽比 (Center/Surface) = %.2f\n', ratio);
    
    if ratio < 0.9
        fprintf('[PASS] Shielding effect observed.\n');
    else
        fprintf('[WARN] Weak shielding. Current ratio suggests fast diffusion.\n');
    end
end
