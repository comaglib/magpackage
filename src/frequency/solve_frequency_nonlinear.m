function [Solution, Info] = solve_frequency_nonlinear(Model, Coil, FreqParams)
% SOLVE_FREQUENCY_NONLINEAR 非线性时谐磁场求解器 (最终版 - 自适应阻尼)
% 
% 方法: Effective Reluctivity (Fixed-Point Iteration)
% 特性:
% 1. 利用重构后的组装器直接注入 Nu_vec，代码高度复用。
% 2. 引入自适应阻尼调整 (Adaptive Relaxation)，确保收敛稳定性。

    fprintf('==============================================\n');
    fprintf('   Nonlinear Frequency Domain Solver          \n');
    fprintf('   (Effective Reluctivity + Adaptive Damp)    \n');
    fprintf('==============================================\n');

    freq = FreqParams.Frequency;
    omega = 2 * pi * freq;
    max_iter = 50; if isfield(FreqParams, 'MaxIter'), max_iter = FreqParams.MaxIter; end
    tol = 1e-4; if isfield(FreqParams, 'Tol'), tol = FreqParams.Tol; end
    
    % 1. 准备不变量 (拓扑与常数矩阵)
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    % 涡流与耦合矩阵 (假设几何不变)
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    M_reg = assemble_mass_matrix(Model, 'All');
    
    % 感应源项 (恒定部分)
    CoilUnit = Coil; CoilUnit.I(:)=1;
    As_unit = project_source_A_on_edges(Model, CoilUnit);
    % b_ind = -jw * sigma * As. 注意 I 是相量
    I_phasor = Coil.I(1);
    b_ind_base = -1i * omega * (M_sigma * (As_unit * I_phasor));
    
    % 初始磁阻率 (Linear)
    Nu_vec = []; % 空表示使用 Model 默认值
    
    % 迭代变量
    x_curr = zeros(numEdges + numActiveV, 1);
    
    % 自适应阻尼参数
    relax = 0.5;        % 初始松弛因子
    prev_diff = inf;    % 上一步误差
    min_relax = 0.01;   % 最小阻尼
    max_relax = 1.0;    % 最大阻尼
    
    % 2. 迭代主循环
    for iter = 1:max_iter
        
        % --- A. 组装刚度矩阵 K(nu) ---
        % 利用重构特性: 直接传入 Nu_vec
        K = assemble_magnetic_stiffness(Model, Nu_vec);
        
        % --- B. 组装磁化源 b_mag(nu) ---
        % RHS 也会随 nu 变化 (nu0 - nu)
        b_mag = assemble_rhs_reduced(Model, Coil, Nu_vec);
        
        % 总 RHS
        SysRHS = [b_mag + b_ind_base; sparse(numActiveV, 1)];
        
        % --- C. 构建系统矩阵 ---
        scale_K = mean(abs(nonzeros(K)));
        eps_A = 1e-6 * scale_K;
        
        K_eff = K + 1i * omega * M_sigma + eps_A * M_reg;
        C_eff = C_sigma; 
        Ct_eff = 1i * omega * C_sigma';
        G_eff = 1i * omega * G_sigma;
        
        SysK = [K_eff, C_eff; Ct_eff, G_eff];
        
        % --- D. 边界条件与求解 ---
        if isfield(Model.Runtime,'FixedEdges'), fixed_edges=Model.Runtime.FixedEdges; else, fixed_edges=[]; end
        [K_bc, R_bc] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, zeros(size(fixed_edges)));
        
        Model.Solver.Linear.Interface='MUMPS'; 
        Model.Solver.Linear.Symmetric=true;
        if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
        Model.Solver.Linear.MumpsICNTL.i8=77;
        
        x_new = linear_solve(K_bc, R_bc, Model);
        
        % --- E. 更新材料属性 ---
        A_new = x_new(1:numEdges);
        
        % 计算 B_eff 并更新 Nu
        [Nu_new, max_B] = update_reluctivity_parallel(Model, A_new, Coil);
        
        % --- F. 检查收敛与自适应更新 ---
        if isempty(Nu_vec)
            % 第一步 (Linear Solution)
            Nu_vec = Nu_new; % 初始化
            diff = 1.0;
            fprintf('    Step 0 (Linear): Max|B|=%.4f T\n', max_B);
        else
            % 计算相对误差
            diff = norm(Nu_new - Nu_vec) / norm(Nu_vec);
            
            % 自适应调整策略
            if diff < prev_diff
                % 误差下降，尝试增大步长 (加速)
                relax = min(relax * 1.1, max_relax);
                status = 'Acc';
            else
                % 误差上升或震荡，减小步长 (稳定)
                relax = max(relax * 0.33, min_relax);
                status = 'Damp';
            end
            
            fprintf('    Iter %2d: Max|B|=%.4f T, d(nu)=%.4e, relax=%.3f [%s]\n', ...
                iter, max_B, diff, relax, status);
            
            if diff < tol
                fprintf('  [Converged]\n');
                x_curr = x_new;
                break;
            end
            
            % 更新 Nu_vec (Relaxed Update)
            Nu_vec = relax * Nu_new + (1-relax) * Nu_vec;
            prev_diff = diff;
        end
        
        x_curr = x_new;
    end
    
    if iter == max_iter
        warning('Nonlinear iteration reached max steps.');
    end
    
    Solution.A = x_curr(1:numEdges);
    Solution.V_active = x_curr(numEdges+1:end);
    Info.Iterations = iter;
    Info.FinalNu = Nu_vec; 
end

function [Nu_new, max_B] = update_reluctivity_parallel(Model, A_phasor, Coil)
    % 并行计算每个单元的有效磁阻率
    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    
    % Pre-chunking
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw=Mesh.T; T2E_raw=Mesh.T2E; Signs_raw=double(Mesh.T2E_Sign);
    MatMap_raw = Model.Materials.ActiveMap;
    
    for k=1:numChunks
        idx = (k-1)*chunkSize+1 : min(k*chunkSize, numElems);
        ElementBlocks{k}.T = T_raw(:, idx);
        ElementBlocks{k}.T2E = T2E_raw(:, idx);
        ElementBlocks{k}.Signs = Signs_raw(:, idx);
        ElementBlocks{k}.MatIDs = MatMap_raw(idx);
        ElementBlocks{k}.Count = length(idx);
    end
    
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(A_phasor);
    C_Coil = parallel.pool.Constant(Coil);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    
    Nu_cell = cell(numChunks, 1);
    B_cell = cell(numChunks, 1);
    
    parfor k=1:numChunks
        Block = ElementBlocks{k};
        cnt = Block.Count;
        loc_P = C_P.Value; loc_A = C_A.Value; loc_Coil = C_Coil.Value; loc_MatLib = C_MatLib.Value;
        
        nu_local = zeros(cnt, 1);
        b_max_local = 0;
        
        % Centers
        Centers = zeros(3, cnt);
        for e=1:cnt
            nodes = Block.T(:,e);
            Centers(:,e) = mean(loc_P(:,nodes), 2);
        end
        Bs_all = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e=1:cnt
            nodes = Block.T(:,e);
            pts = loc_P(:, nodes);
            edges = Block.T2E(:,e);
            s = Block.Signs(:,e);
            
            a_vec = loc_A(edges) .* s;
            Br = calc_element_B(pts, a_vec);
            
            B_tot = Br + Bs_all(:,e);
            
            % 计算有效值 B_eff (假设 BH 曲线定义为 Peak B)
            B_mag = norm([abs(B_tot(1)), abs(B_tot(2)), abs(B_tot(3))]);
            
            if B_mag > b_max_local, b_max_local = B_mag; end
            
            mat_id = Block.MatIDs(e);
            mat_info = loc_MatLib(mat_id);
            
            [nu, ~] = eval_material_nu(B_mag^2, mat_info);
            nu_local(e) = nu;
        end
        Nu_cell{k} = nu_local;
        B_cell{k} = b_max_local;
    end
    
    Nu_new = vertcat(Nu_cell{:});
    max_B = max([B_cell{:}]);
end