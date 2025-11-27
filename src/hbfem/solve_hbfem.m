function [Solution, Info] = solve_hbfem(Model, Coil, HBFEMParams)
% SOLVE_HBFEM 谐波平衡法求解器 (HBFEM) - 自适应高性能版
%
% 更新:
% 1. 自适应计算 TimeSteps，确保满足奈奎斯特采样定理。
% 2. 优化内存与广播变量。

    fprintf('==============================================\n');
    fprintf('   Harmonic Balance FEM Solver (Adaptive)     \n');
    fprintf('   Method: Picard Iteration w/ AFT            \n');
    fprintf('==============================================\n');

    base_freq = HBFEMParams.Frequency;
    harmonics = HBFEMParams.Harmonics; 
    
    % --- [Update 1] 自适应 TimeSteps 设置 ---
    max_h = max(harmonics);
    min_required = 2 * max_h + 1;
    
    if ~isfield(HBFEMParams, 'TimeSteps') || isempty(HBFEMParams.TimeSteps) || HBFEMParams.TimeSteps < min_required
        % 自动选择: 取大于最小需求的下一个 2 的幂次，且至少为 32
        % (2的幂次对 FFT 效率最高)
        n_time = max(32, 2^nextpow2(min_required + 2));
        fprintf('  [Auto-Setup] TimeSteps adjusted to %d (Min req: %d)\n', n_time, min_required);
    else
        n_time = HBFEMParams.TimeSteps;
    end
    % ----------------------------------------
    
    numHarm = length(harmonics);
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 1. 初始化 AFT 模块
    AFT = aft_module();
    Info_A = AFT.prepare_indices(harmonics, n_time);
    
    % Info_Nu 用于磁阻率 nu (偶次谐波)
    harmonics_nu = 0 : 2 : (2*max_h); 
    Info_Nu = AFT.prepare_indices(harmonics_nu, n_time);
    numHarmNu = length(harmonics_nu);
    
    % 2. 恒定矩阵
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    M_reg = assemble_mass_matrix(Model, 'All');
    
    % 3. 源电流设置
    I_vec = zeros(numHarm, 1);
    if isfield(Coil, 'HarmonicCurrents')
        I_vec = Coil.HarmonicCurrents;
    else
        if isfield(Coil, 'I') && ~isempty(Coil.I)
            I_vec(1) = Coil.I(1); 
        else
            I_vec(1) = 1.0; 
        end
    end
    
    if isfield(Coil, 'Turns'), N_turns = Coil.Turns; else, N_turns = 1.0; end
    CoilUnit = Coil; CoilUnit.I(:) = 1.0 * N_turns;
    As_unit = project_source_A_on_edges(Model, CoilUnit);
    
    % 4. 边界条件预处理 (向量化)
    if isfield(Model.Runtime,'FixedEdges'), fe=Model.Runtime.FixedEdges; else, fe=[]; end
    fixed_dofs_all = [];
    if ~isempty(fe)
        offsets = (0:numHarm-1) * numEdges; 
        fe_mat = repmat(fe, 1, numHarm) + repmat(offsets, length(fe), 1);
        fixed_dofs_all = fe_mat(:);
    end
    fixed_vals_all = zeros(size(fixed_dofs_all));
    
    % 5. 初始化解向量
    X_curr = zeros(numEdges * numHarm, 1); 
    
    % 初始磁阻率 (仅 DC 分量非零)
    numElems = size(Model.Mesh.T, 2);
    Nu_fft_curr = zeros(numElems, numHarmNu); 
    
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    mu0 = 4*pi*1e-7;
    
    dc_col_idx = 1; % Info_Nu 中 0Hz 在第一列
    for i = 1:numElems
        if MatMap(i) <= length(MatLib)
            mat = MatLib(MatMap(i));
            Nu_fft_curr(i, dc_col_idx) = 1.0 / (mat.Mu_r * mu0);
        else
             Nu_fft_curr(i, dc_col_idx) = 1.0 / mu0;
        end
    end
    
    MaxIter = 20; Tol = 1e-4;
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true; 
    if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
    Model.Solver.Linear.MumpsICNTL.i8 = 77;
    
    alpha = 1.0; err_prev = inf;
    
    fprintf('  [HBFEM] Starting Iterations...\n');
    
    for iter = 1:MaxIter
        % (A) 组装 LHS
        K_toep = assemble_harmonic_stiffness(Model, Nu_fft_curr, harmonics, Info_Nu);
        
        % 添加涡流项
        scale_K = mean(abs(nonzeros(K_toep)));
        eps_val = 1e-6 * scale_K;
        
        for k = 1:numHarm
            h_order = harmonics(k);
            w_k = h_order * 2 * pi * base_freq;
            idx = (k-1)*numEdges + 1 : k*numEdges;
            M_term = 1i * w_k * M_sigma + eps_val * M_reg;
            K_toep(idx, idx) = K_toep(idx, idx) + M_term;
        end
        
        % (B) AFT 更新 & RHS
        [Nu_fft_new, RHS_mag_vec, max_B] = ...
            update_hbfem_properties(Model, X_curr, Coil, AFT, Info_A, Info_Nu, I_vec, base_freq);
            
        % (C) 组装总 RHS
        RHS_total = RHS_mag_vec;
        
        for k = 1:numHarm
            h_order = harmonics(k);
            w_k = h_order * 2 * pi * base_freq;
            idx = (k-1)*numEdges + 1 : k*numEdges;
            if abs(I_vec(k)) > 1e-9
                b_ind_k = -1i * w_k * M_sigma * (As_unit * I_vec(k));
                RHS_total(idx) = RHS_total(idx) + b_ind_k;
            end
        end
        
        % (D) 求解
        [K_solve, R_solve] = apply_dirichlet_bc(K_toep, RHS_total, fixed_dofs_all, fixed_vals_all);
        
        try
            X_new = linear_solve(K_solve, R_solve, Model);
        catch
            X_new = K_solve \ R_solve;
        end
        
        % (E) 收敛检查
        err = norm(X_new - X_curr) / (norm(X_new) + 1e-10);
        fprintf('    Iter %2d: Max|B|=%.2fT, Err=%.2e (Relax=%.2f)\n', iter, max_B, err, alpha);
        
        if err < Tol
            fprintf('  [Converged]\n');
            X_curr = X_new;
            break;
        end
        
        if iter > 1
            if err > err_prev, alpha = max(alpha*0.5, 0.1); else, alpha = min(alpha*1.1, 1.0); end
        end
        
        X_curr = (1-alpha)*X_curr + alpha*X_new;
        Nu_fft_curr = (1-alpha)*Nu_fft_curr + alpha*Nu_fft_new; 
        err_prev = err;
    end
    
    Solution.Harmonics = reshape(X_curr, numEdges, numHarm);
    Info.Iterations = iter;
    Info.TimeSteps = n_time; % 返回实际使用的步数
end

function [Nu_fft, RHS_vec, max_B] = update_hbfem_properties(Model, X_curr, Coil, AFT, Info_A, Info_Nu, I_vec, freq)
    
    Mesh = Model.Mesh;
    numEdges = size(Mesh.Edges, 2);
    numElems = size(Mesh.T, 2);
    numHarmA = Info_A.NumHarmonics;
    numHarmNu = Info_Nu.NumHarmonics;
    Nt = Info_A.Nt; 
    
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw=Mesh.T; T2E_raw=Mesh.T2E; Signs_raw=double(Mesh.T2E_Sign); MatMap_raw=Model.Materials.ActiveMap;
    for k=1:numChunks
        idx=(k-1)*chunkSize+1 : min(k*chunkSize, numElems);
        ElementBlocks{k}.T=T_raw(:,idx); ElementBlocks{k}.T2E=T2E_raw(:,idx);
        ElementBlocks{k}.Signs=Signs_raw(:,idx); ElementBlocks{k}.MatIDs=MatMap_raw(idx);
        ElementBlocks{k}.Count=length(idx);
    end
    
    A_harm = reshape(X_curr, numEdges, numHarmA);
    
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(A_harm);
    C_Coil = parallel.pool.Constant(Coil);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    C_AFT = parallel.pool.Constant(AFT);
    
    Nu_cell = cell(numChunks, 1);
    RHS_cell = cell(numChunks, 1);
    B_cell = cell(numChunks, 1);
    
    I_t = zeros(1, Nt);
    for k = 1:numHarmA
        coeffs = zeros(1, numHarmA);
        coeffs(k) = I_vec(k);
        I_t = I_t + AFT.freq2time(coeffs, Info_A);
    end
    C_It = parallel.pool.Constant(I_t);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        cnt = Block.Count;
        loc_P = C_P.Value; loc_A = C_A.Value; loc_Coil = C_Coil.Value; 
        loc_MatLib = C_MatLib.Value; loc_It = C_It.Value; loc_AFT = C_AFT.Value;
        
        nu_fft_local = zeros(cnt, numHarmNu); 
        rhs_elem_data = zeros(6, numHarmA, cnt);
        b_max_local = 0;
        
        Centers = zeros(3, cnt);
        for e=1:cnt, nodes=Block.T(:,e); Centers(:,e)=mean(loc_P(:,nodes),2); end
        Bs_unit = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e = 1:cnt
            nodes = Block.T(:,e); pts = loc_P(:, nodes); edges = Block.T2E(:,e); s = Block.Signs(:,e);
            
            a_harm_elem = loc_A(edges, :) .* s; 
            a_t = loc_AFT.freq2time(a_harm_elem, Info_A); 
            
            Br_t = calc_element_B_time(pts, a_t);
            Bs_t = Bs_unit(:, e) * loc_It;
            B_tot_t = Br_t + Bs_t;
            
            B_mag_t = sqrt(sum(B_tot_t.^2, 1));
            if max(B_mag_t) > b_max_local, b_max_local = max(B_mag_t); end
            
            mat_id = Block.MatIDs(e); mat_info = loc_MatLib(mat_id);
            [nu_t, ~] = eval_material_nu(B_mag_t.^2, mat_info);
            
            nu_fft_local(e, :) = loc_AFT.time2freq(nu_t, Info_Nu);
            
            mu0 = 4*pi*1e-7; nu0 = 1/mu0;
            M_t = (nu0 - nu_t) .* Bs_t; 
            M_harm = loc_AFT.time2freq(M_t, Info_A);
            
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            J = [v21, v31, v41]; detJ = det(J); invJ_T = inv(J)';
            G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1]; G_phy = invJ_T * G_ref;
            
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            Curl_N = zeros(3, 6);
            for i_edge = 1:6
                na = edge_pairs(i_edge,1); nb = edge_pairs(i_edge,2);
                Curl_N(:, i_edge) = 2 * cross(G_phy(:, na), G_phy(:, nb));
            end
            Vol = abs(detJ) / 6.0;
            
            rhs_val = Vol * (Curl_N' * M_harm); 
            rhs_val = rhs_val .* s; 
            rhs_elem_data(:, :, e) = rhs_val;
        end
        Nu_cell{k} = nu_fft_local;
        RHS_cell{k} = rhs_elem_data;
        B_cell{k} = b_max_local;
    end
    
    Nu_fft = vertcat(Nu_cell{:});
    max_B = max([B_cell{:}]);
    
    RHS_vec = zeros(numEdges * numHarmA, 1);
    for k = 1:numChunks
        Block = ElementBlocks{k};
        rhs_data = RHS_cell{k}; 
        for h = 1:numHarmA
            vals = squeeze(rhs_data(:, h, :)); 
            rows = Block.T2E; 
            offset = (h-1)*numEdges;
            S = sparse(rows(:), 1, vals(:), numEdges, 1);
            RHS_vec(offset+1 : offset+numEdges) = RHS_vec(offset+1 : offset+numEdges) + full(S);
        end
    end
end

function B_t = calc_element_B_time(P, A_t)
    [~, Nt] = size(A_t);
    v21=P(:,2)-P(:,1); v31=P(:,3)-P(:,1); v41=P(:,4)-P(:,1);
    J = [v21, v31, v41]; invJ_T = inv(J)';
    G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1]; G_phy = invJ_T * G_ref;
    edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    Curl_N = zeros(3, 6);
    for i=1:6
        na=edge_pairs(i,1); nb=edge_pairs(i,2);
        Curl_N(:,i) = 2 * cross(G_phy(:,na), G_phy(:,nb));
    end
    B_t = Curl_N * A_t;
end