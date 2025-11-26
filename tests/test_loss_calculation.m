function test_loss_calculation()
% TEST_LOSS_CALCULATION 验证损耗和能量计算 (修正版)
% 
% 修复: 修正了 compute_ohmic_loss 的调用参数顺序错误。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Loss & Energy Calculation Test             \n');
    fprintf('==============================================\n');
    
    % 1. 准备环境
    msh_file = 'test_loss.msh';
    create_delaunay_mesh_robust(msh_file); 
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 材料
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Conductor';
    Model.Materials.Lib(2).Type = 'Linear';
    Model.Materials.Lib(2).Mu_r = 1.0;
    Model.Materials.Lib(2).Sigma = 1e6; 
    
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.4; 
    CoilParams.Radius = 0.4; 
    CoilParams.Current = 0; 
    CoilParams.N_seg = 72;
    Coil = create_racetrack_coil(CoilParams);
    
    % 时间设置
    TargetI = 1e5;
    TimeParams.I_func = @(t) TargetI * (t > 0);
    TimeParams.dt = 0.005;
    TimeParams.N_steps = 10;
    
    Model.Runtime.FixedEdges = identify_boundary_edges(Model.Mesh, Raw, 100);
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 2. 求解循环 (手动展开)
    fprintf('  [Init] Assembling matrices...\n');
    
    DofData = create_dof_map(Model);
    K = assemble_magnetic_stiffness(Model);
    M = assemble_mass_matrix(Model, 'ConductorOnly');
    G = assemble_scalar_laplacian(Model, DofData);
    C = assemble_coupling_matrix(Model, 'Physical', DofData);
    
    M_reg = assemble_mass_matrix(Model, 'All');
    eps_A = 1e-8 * mean(abs(nonzeros(K)));
    
    dt = TimeParams.dt;
    
    K_eff = K + (1/dt)*M + eps_A*M_reg;
    Ct = C';
    G_eff = G * dt;
    
    SysK = [K_eff, C; Ct, G_eff];
    
    CoilUnit = Coil; CoilUnit.I(:)=1;
    As_unit = project_source_A_on_edges(Model, CoilUnit);
    
    A_prev = zeros(DofData.NumEdges, 1);
    I_prev = 0;
    
    fixed_vals_A = zeros(length(Model.Runtime.FixedEdges), 1);
    Model.Solver.Linear.Interface='MUMPS'; 
    Model.Solver.Linear.Symmetric=true;
    if ~isfield(Model.Solver.Linear, 'MumpsICNTL')
        Model.Solver.Linear.MumpsICNTL = struct();
    end
    Model.Solver.Linear.MumpsICNTL.i8 = 77;
    
    LossHistory = [];
    EnergyHistory = [];
    
    fprintf('  [Test] Stepping and calculating metrics...\n');
    
    for n = 1:TimeParams.N_steps
        t = n*dt; 
        I_now = TargetI; 
        
        % RHS
        Coil.I(:) = I_now;
        b_mag = assemble_rhs_reduced(Model, Coil);
        
        dI = I_now - I_prev;
        b_ind = - M * (As_unit * dI/dt);
        b_ext = b_mag + b_ind;
        
        rhs_A = b_ext + (1/dt)*M*A_prev;
        rhs_V = Ct * A_prev;
        
        [K_bc, R_bc] = apply_dirichlet_bc(SysK, [rhs_A; rhs_V], Model.Runtime.FixedEdges, fixed_vals_A);
        
        x = linear_solve(K_bc, R_bc, Model);
        
        A_curr = x(1:DofData.NumEdges);
        V_curr = x(DofData.NumEdges+1:end);
        
        % --- 计算指标 ---
        % [核心修正] 参数顺序: Model, A, V, A_prev, dt
        [P_ohm, ~] = compute_ohmic_loss(Model, A_curr, V_curr, A_prev, dt);
        
        W_mag = compute_magnetic_energy(Model, A_curr, Coil);
        
        LossHistory(n) = P_ohm; %#ok<AGROW>
        EnergyHistory(n) = W_mag; %#ok<AGROW>
        
        fprintf('    Step %d: Loss=%.2e W, Energy=%.2e J\n', n, P_ohm, W_mag);
        
        A_prev = A_curr;
        I_prev = I_now;
    end
    
    % 3. 验证
    loss_ratio = LossHistory(end) / LossHistory(1);
    is_decaying = loss_ratio < 0.01; 
    
    energy_change = abs(EnergyHistory(end) - EnergyHistory(1)) / EnergyHistory(1);
    is_stable = energy_change < 0.05; 
    
    fprintf('\n--- Validation ---\n');
    fprintf('Loss Decay Ratio: %.2e (Pass if < 0.01)\n', loss_ratio);
    fprintf('Energy Change:    %.2f%% (Pass if < 5%%)\n', energy_change*100);
    
    if is_decaying && is_stable
        fprintf('[PASS] Physics trends are reasonable.\n');
    else
        if ~is_decaying
            warning('Loss did not decay significantly.');
        end
        if ~is_stable
            warning('Energy is unstable.');
        end
    end
end

function create_delaunay_mesh_robust(filename)
    % 生成鲁棒的四面体网格
    % fprintf('  [MeshGen] Generating mesh via Delaunay (Robust)...\n');
    
    [X1, Y1, Z1] = ndgrid(linspace(-0.15, 0.15, 6));
    P_inner = [X1(:), Y1(:), Z1(:)];
    [X2, Y2, Z2] = ndgrid(linspace(-0.4, 0.4, 5));
    P_outer = [X2(:), Y2(:), Z2(:)];
    
    P = unique([P_inner; P_outer], 'rows');
    nNodes = size(P, 1);
    
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList; 
    
    % 手性修正
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    v1 = p2 - p1; v2 = p3 - p1; v3 = p4 - p1;
    vols = sum(cross(v1, v2, 2) .* v3, 2);
    if any(vols < 0), T(vols < 0, [1 2]) = T(vols < 0, [2 1]); end
    
    P_io = P'; T_io = T';
    nElems = size(T_io, 2);
    
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    in_cond = abs(Centers(:,1)) < 0.16 & abs(Centers(:,2)) < 0.16 & abs(Centers(:,3)) < 0.16;
    ElemTags = ones(1, nElems);
    ElemTags(in_cond) = 2;
    
    TR = triangulation(T, P);
    FB = freeBoundary(TR);
    BoundaryFaces = FB';
    nFaces = size(BoundaryFaces, 2);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    fprintf(fid, '%d %.6f %.6f %.6f\n', [(1:nNodes); P_io]);
    fprintf(fid, '$EndNodes\n');
    
    nTotal = nFaces + nElems;
    idx_air = (ElemTags == 1); idx_cond = (ElemTags == 2);
    nAir = sum(idx_air); nCond = sum(idx_cond);
    nBlocks = 1 + (nAir>0) + (nCond>0);
    
    fprintf(fid, '$Elements\n%d %d 1 %d\n', nBlocks, nTotal, nTotal);
    fprintf(fid, '2 100 2 %d\n', nFaces);
    fprintf(fid, '%d %d %d %d\n', [(1:nFaces); BoundaryFaces]);
    
    id = nFaces + 1;
    if nAir > 0
        fprintf(fid, '3 1 4 %d\n', nAir);
        fprintf(fid, '%d %d %d %d %d\n', [(id:id+nAir-1); T_io(:, idx_air)]);
        id = id + nAir;
    end
    if nCond > 0
        fprintf(fid, '3 2 4 %d\n', nCond);
        fprintf(fid, '%d %d %d %d %d\n', [(id:id+nCond-1); T_io(:, idx_cond)]);
    end
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end