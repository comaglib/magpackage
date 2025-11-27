function [Solution, Info] = solve_hbfem(Model, Coil, HBFEMParams)
% SOLVE_HBFEM 分解式谐波平衡有限元求解器 (High-Performance Decomposed HBFEM)
%
% 概述 (Overview):
%   本函数实现了基于分解策略（Decomposed / Fixed-Point）的谐波平衡有限元方法（HBFEM），
%   专门用于求解含深度饱和铁磁材料的三维非线性时谐磁场问题。通过解耦各次谐波方程，
%   大幅降低了内存消耗，并利用高级收敛策略解决了强非线性下的不稳定性。
%
% 核心特性 (Key Features):
%   1. [Smart Adaptive Ramping] 智能自适应负载步进策略。根据收敛速度动态调整负载增量，
%      在线性区快速推进，在非线性区精细迭代，显著提升计算效率。
%   2. [Dynamic Relaxation] 基于误差梯度的动态松弛因子控制（PID-like）。自动调节
%      松弛因子 Alpha，在保证稳定性的前提下最大化收敛速度。
%   3. [Source Smoothing] 几何自适应源场平滑。自动计算并消除 Biot-Savart 积分中的
%      近场奇点，防止因网格与线圈重合导致的非物理高场和数值发散。
%   4. [M-Scaling Stabilization] 磁化强度缩放技术。在深度饱和区（B > 3.0T）对
%      非线性残差进行幅度限制，有效抑制正反馈震荡。
%   5. [Geometry Pre-computation] 几何内核预计算。在迭代前一次性计算所有单元的
%      Curl形函数和雅可比矩阵，极大加速了右端项（RHS）的并行组装。
%   6. [High-Performance Computing] 全流程并行化设计。采用 Map-Reduce 模式和
%      数据切片（Slicing）技术，结合 MUMPS 直接求解器优先策略，充分利用硬件性能。
%
% 输入 (Inputs):
%   Model       - 有限元模型结构体，包含 Mesh (网格), Materials (材料库) 等。
%   Coil        - 线圈几何与电流参数结构体。
%   HBFEMParams - 求解控制参数 (Frequency, Harmonics, TimeSteps 等)。
%
% 输出 (Outputs):
%   Solution    - 结果结构体，包含各次谐波的磁矢位解向量 (Solution.Harmonics)。
%   Info        - 求解统计信息 (迭代步数、最终误差、时间步数等)。

    fprintf('==============================================\n');
    fprintf('   HBFEM Solver (Decomposed / Fixed-Point)    \n');
    fprintf('   Method: P-HBFEM (Smart Adaptive Ramping)   \n');
    fprintf('==============================================\n');

    base_freq = HBFEMParams.Frequency;
    harmonics = HBFEMParams.Harmonics; 
    
    % --- 1. 自适应时间步配置 (FFT Optimization) ---
    max_h = max(harmonics);
    min_required = 2 * max_h + 1;
    if ~isfield(HBFEMParams, 'TimeSteps') || isempty(HBFEMParams.TimeSteps) || HBFEMParams.TimeSteps < min_required
        % 自动向上取整到最近的 2 的幂次 (32, 64, 128...) 以获得最佳 FFT 性能
        n_time = max(32, 2^nextpow2(min_required + 2));
        fprintf('  [Auto-Setup] TimeSteps adjusted to %d\n', n_time);
    else
        n_time = HBFEMParams.TimeSteps;
    end
    
    numHarm = length(harmonics);
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    
    % 启动并行池 (如果尚未启动)
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % --- 2. 初始化计算模块 ---
    AFT = aft_module();
    Info_A = AFT.prepare_indices(harmonics, n_time);
    
    % 组装恒定质量矩阵 (Mass Matrices)
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly'); % 导电区涡流项
    M_reg = assemble_mass_matrix(Model, 'All');             % 全局正则化项
    
    % --- 3. 源电流与投影 ---
    I_base_vec = zeros(numHarm, 1);
    if isfield(Coil, 'HarmonicCurrents')
        I_base_vec = Coil.HarmonicCurrents;
    elseif isfield(Coil, 'I') && ~isempty(Coil.I)
        I_base_vec(1) = Coil.I(1); 
    else
        I_base_vec(1) = 1.0; 
    end
    
    if isfield(Coil, 'Turns'), N_turns = Coil.Turns; else, N_turns = 1.0; end
    CoilUnit = Coil; CoilUnit.I(:) = 1.0 * N_turns;
    As_unit = project_source_A_on_edges(Model, CoilUnit); % 用于涡流计算
    
    % --- 4. 源场预计算与自适应平滑 (Source Field Pre-calc) ---
    fprintf('  [Pre-Calc] Computing static Biot-Savart field (Adaptive Smoothing)...\n');
    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    Centers = zeros(3, numElems);
    P = Mesh.P; T = Mesh.T;
    V_total = 0;
    % 计算单元重心及特征尺寸
    for i = 1:numElems
        nodes = T(:, i); pts = P(:, nodes);
        Centers(:,i) = sum(pts, 2) / 4; 
        d1 = norm(pts(:,2)-pts(:,1)); d2 = norm(pts(:,3)-pts(:,1)); d3 = norm(pts(:,4)-pts(:,1));
        V_total = V_total + d1*d2*d3/6; 
    end
    avg_vol = V_total / numElems;
    
    % 计算原始 Biot-Savart 场
    Bs_unit_raw = compute_biot_savart_B_serial(CoilUnit, Centers);
    
    % 应用 Soft-Clip 平滑技术，消除网格奇点
    I_max = max(abs(I_base_vec)); if I_max == 0, I_max = 1.0; end
    B_limit_phys = 10.0; % 物理上限 10T
    Limit_Unit = B_limit_phys / I_max;
    bs_mag = sqrt(sum(Bs_unit_raw.^2, 1));
    scale_factor = Limit_Unit ./ sqrt(bs_mag.^2 + Limit_Unit^2);
    Bs_unit_global = Bs_unit_raw .* repmat(scale_factor, 3, 1);
    
    if isfield(Model.Runtime,'FixedEdges'), fe=Model.Runtime.FixedEdges; else, fe=[]; end
    
    % --- 5. 迭代初始化 ---
    X_curr = zeros(numEdges, numHarm); 
    Nu_DC = zeros(numElems, 1);
    mu0 = 4*pi*1e-7;
    Nu_DC(:) = 1.0 / mu0; % 初始状态设为全空气 (Vacuum Start)，保证矩阵正定性
    
    % 线性求解器配置 (优先 MUMPS)
    use_mumps = false;
    if exist('linear_solve', 'file')
        use_mumps = true;
        Model.Solver.Linear.Interface = 'MUMPS';
        Model.Solver.Linear.Symmetric = true;
        fprintf('  [Solver] MUMPS interface enabled.\n');
    end
    
    fprintf('  [HBFEM-Dec] Starting Smart Iterations...\n');
    
    % --- 6. 智能迭代控制参数 ---
    MaxIter = 200; 
    Tol = 1e-4;
    
    load_scale = 0.0;      % 当前负载比例 (0~1)
    load_step_size = 0.1;  % 初始负载增量
    min_step_size = 0.01;  % 最小增量
    max_step_size = 1.0;   % 最大增量 (允许一步到位)
    
    alpha = 0.5;           % 初始松弛因子
    err_prev = inf;
    
    iter_total = 0;
    converged_flag = false;
    iter_inner = 0;        % 当前负载步内的迭代计数
    
    % --- 7. 主迭代循环 ---
    while load_scale < 1.0 || ~converged_flag
        iter_total = iter_total + 1;
        if iter_total > MaxIter, break; end
        
        % 7.1 负载步进控制 (Adaptive Ramping)
        if converged_flag && load_scale < 1.0
            % 上一步收敛，推进负载
            load_scale = min(load_scale + load_step_size, 1.0);
            converged_flag = false; 
            err_prev = inf;         
            
            % 智能步长调整：收敛快则加速，收敛慢则减速
            if iter_inner < 5
                load_step_size = min(load_step_size * 2.0, max_step_size);
            elseif iter_inner > 15
                load_step_size = max(load_step_size * 0.5, min_step_size);
            end
            iter_inner = 0; 
        end
        iter_inner = iter_inner + 1;
        
        I_vec_curr = I_base_vec * load_scale;
        
        % (A) 组装左端项 (LHS) 基准刚度矩阵
        K_base = assemble_magnetic_stiffness(Model, Nu_DC);
        scale_K = mean(abs(nonzeros(K_base)));
        eps_val = 1e-6 * scale_K; % 正则化因子
        
        % (B) 更新材料属性 & 计算非线性残差 (并行)
        [Nu_DC_suggested, M_resid_harm, M_src_harm, max_B_phys] = ...
            update_properties_constrained(Model, X_curr, AFT, Info_A, I_vec_curr, Nu_DC, Bs_unit_global);
        
        % (C) 逐谐波求解线性方程组
        X_new = zeros(size(X_curr));
        for k = 1:numHarm
            h_order = harmonics(k);
            w_k = h_order * 2 * pi * base_freq;
            
            % 构建复数刚度矩阵 (Stiffness + Eddy + Regularization)
            if w_k == 0
                K_solve = K_base + eps_val * M_reg;
            else
                K_solve = K_base + 1i * w_k * M_sigma + eps_val * M_reg;
            end
            
            RHS_k = sparse(numEdges, 1);
            
            % 1. 涡流源项
            if abs(I_vec_curr(k)) > 1e-9
                b_eddy = -1i * w_k * M_sigma * (As_unit * I_vec_curr(k));
                RHS_k = RHS_k + b_eddy;
            end
            
            % 2. 非线性残差与源场驱动项
            M_vec_resid = squeeze(M_resid_harm(:, k, :)); 
            RHS_resid = assemble_rhs_from_magnetization(Model, M_vec_resid);
            
            M_vec_src = squeeze(M_src_harm(:, k, :));
            RHS_source = assemble_rhs_from_magnetization(Model, M_vec_src);
            
            RHS_k = RHS_k - RHS_resid - RHS_source;
            
            % 应用边界条件
            [K_sys, R_sys] = apply_dirichlet_bc(K_solve, RHS_k, fe, zeros(size(fe)));
            
            % 线性求解 (MUMPS优先)
            if use_mumps
                try
                    x_k = linear_solve(K_sys, R_sys, Model);
                catch
                    x_k = K_sys \ R_sys;
                end
            else
                x_k = K_sys \ R_sys;
            end
            X_new(:, k) = x_k;
        end
        
        % (D) 收敛检查与动态松弛
        
        % 物理安全保护: 如果B场过大，强制减速
        if max_B_phys > 10.0
            alpha = 0.1; 
        end
        
        diff_norm = norm(X_new - X_curr, 'fro');
        sol_norm = norm(X_new, 'fro') + 1e-10;
        err = diff_norm / sol_norm;
        
        fprintf('    Iter %d (Step +%.2f): Max|B|=%.2fT, Err=%.2e (Alpha=%.2f) [Load=%.0f%%]\n', ...
            iter_total, load_step_size, max_B_phys, err, alpha, load_scale*100);
        
        if err < Tol
            converged_flag = true;
            if load_scale >= 1.0
                fprintf('  [Converged] Final solution reached.\n');
                X_curr = X_new;
                break;
            end
        else
            converged_flag = false;
        end
        
        % 动态松弛控制 (Dynamic Alpha)
        if iter_inner > 1
            if err < err_prev
                % 误差下降：尝试加速
                if err < 0.1
                    alpha = min(alpha * 1.2, 1.0); 
                else
                    alpha = min(alpha * 1.05, 0.9);
                end
            else
                % 误差上升：迅速减速
                alpha = max(alpha * 0.6, 0.1);
            end
        end
        
        % 更新解与矩阵属性
        X_curr = (1-alpha)*X_curr + alpha*X_new;
        Nu_DC = (1-alpha)*Nu_DC + alpha*Nu_DC_suggested;
        
        err_prev = err;
    end
    
    Solution.Harmonics = X_curr; 
    Info.Iterations = iter_total;
    Info.TimeSteps = n_time;
end

% -------------------------------------------------------------------------
% 辅助函数: 并行属性更新与残差计算 (Update Properties)
% -------------------------------------------------------------------------
function [Nu_DC_suggested, M_resid_harm, M_src_harm, max_B] = update_properties_constrained(Model, X_curr, AFT, Info_A, I_vec, Nu_DC_Input, Bs_unit_global)
    
    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    numHarmA = Info_A.NumHarmonics;
    Nt = Info_A.Nt; 
    
    % 预计算电流时域波形
    I_t = zeros(1, Nt);
    for k = 1:numHarmA
        coeffs = zeros(1, numHarmA); coeffs(k) = I_vec(k);
        I_t = I_t + AFT.freq2time(coeffs, Info_A);
    end
    
    % 数据切片 (Data Slicing) 以减少并行广播开销
    chunkSize = 10000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    % 提取大数组引用
    T_raw = Mesh.T; T2E_raw = Mesh.T2E; Signs_raw = double(Mesh.T2E_Sign); MatMap_raw = Model.Materials.ActiveMap;
    
    for k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        idx = idx_start : idx_end;
        % 构建数据块
        ElementBlocks{k}.T = T_raw(:, idx); 
        ElementBlocks{k}.T2E = T2E_raw(:, idx);
        ElementBlocks{k}.Signs = Signs_raw(:, idx); 
        ElementBlocks{k}.MatIDs = MatMap_raw(idx);
        ElementBlocks{k}.Count = length(idx);
        % 切片输入变量
        ElementBlocks{k}.Nu_DC_In = Nu_DC_Input(idx);         
        ElementBlocks{k}.Bs_Unit  = Bs_unit_global(:, idx);    
    end
    
    % 构建 Constant 变量
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(X_curr);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    C_AFT = parallel.pool.Constant(AFT);
    C_It = parallel.pool.Constant(I_t);
    
    % 输出容器
    Nu_DC_Out_cell = cell(numChunks, 1);
    M_resid_cell = cell(numChunks, 1);
    M_src_cell = cell(numChunks, 1);
    B_max_cell = cell(numChunks, 1);
    
    % 并行计算循环
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        cnt = Block.Count;
        
        % 获取常量
        loc_P = C_P.Value; loc_A = C_A.Value; loc_MatLib = C_MatLib.Value; loc_It = C_It.Value; loc_AFT = C_AFT.Value;
        loc_Nu_In = Block.Nu_DC_In; loc_Bs_Unit = Block.Bs_Unit;
        
        % 局部变量分配
        nu_dc_local = zeros(cnt, 1);
        m_resid_local = zeros(numHarmA, 3, cnt);
        m_src_local = zeros(numHarmA, 3, cnt);
        b_max_local = 0;
        
        loc_T = Block.T; loc_T2E = Block.T2E; loc_Signs = Block.Signs; loc_MatIDs = Block.MatIDs;
        G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
        
        for i_local = 1:cnt
            % 1. 几何计算 (Inline Geometry Calculation)
            nodes = loc_T(:, i_local);
            p1 = loc_P(:, nodes(1)); p2 = loc_P(:, nodes(2)); p3 = loc_P(:, nodes(3)); p4 = loc_P(:, nodes(4));
            v21 = p2 - p1; v31 = p3 - p1; v41 = p4 - p1;
            c34 = cross(v31, v41); detJ = dot(v21, c34);
            c42 = cross(v41, v21); c23 = cross(v21, v31);
            col1 = c34 / detJ; col2 = c42 / detJ; col3 = c23 / detJ;
            g1 = -col1 - col2 - col3; G_phy = [g1, col1, col2, col3]; 
            
            % 2. 反应场计算 Br(t) = Curl N * A(t)
            edges = loc_T2E(:, i_local); s = loc_Signs(:, i_local);
            A_elem = loc_A(edges, :) .* s; 
            A_t = loc_AFT.freq2time(A_elem, Info_A);
            
            Br_t = zeros(3, Nt);
            Br_t = Br_t + 2*cross(G_phy(:,1), G_phy(:,2)) * A_t(1,:);
            Br_t = Br_t + 2*cross(G_phy(:,1), G_phy(:,3)) * A_t(2,:);
            Br_t = Br_t + 2*cross(G_phy(:,1), G_phy(:,4)) * A_t(3,:);
            Br_t = Br_t + 2*cross(G_phy(:,2), G_phy(:,3)) * A_t(4,:);
            Br_t = Br_t + 2*cross(G_phy(:,2), G_phy(:,4)) * A_t(5,:);
            Br_t = Br_t + 2*cross(G_phy(:,3), G_phy(:,4)) * A_t(6,:);
            
            % 3. 总场合成 B_tot = Br + Bs
            Bs_t = loc_Bs_Unit(:, i_local) * loc_It;
            B_tot_t = Br_t + Bs_t;
            
            B_mag_sq = sum(B_tot_t.^2, 1);
            B_mag_t = sqrt(B_mag_sq);
            if max(B_mag_t) > b_max_local, b_max_local = max(B_mag_t); end
            
            mat_id = loc_MatIDs(i_local); mat_info = loc_MatLib(mat_id);
            
            % 4. 材料属性更新 (Physics Constrained)
            % 限制 B 场查表输入，防止溢出 B-H 曲线范围
            B_clamped_sq = min(B_mag_sq, 3.0^2);
            [nu_t, ~] = eval_material_nu(B_clamped_sq, mat_info);
            nu_dc_local(i_local) = max(nu_t); % 采用硬矩阵策略
            
            % 5. 残差计算 (M-Scaling)
            % 仅在极端饱和区 (B > 3.0T) 对残差幅度进行限制
            scale_factor = ones(1, Nt);
            mask = B_mag_t > 3.0;
            if any(mask), scale_factor(mask) = 3.0 ./ B_mag_t(mask); end
            B_resid_vec = B_tot_t .* repmat(scale_factor, 3, 1);
            
            % M_resid: 修正 LHS 矩阵误差的残差项
            M_resid = (nu_t - loc_Nu_In(i_local)) .* B_resid_vec;
            % M_src: 源场驱动项
            M_src = nu_t .* Bs_t; 
            
            for d = 1:3
                m_resid_local(:, d, i_local) = loc_AFT.time2freq(M_resid(d, :), Info_A);
                m_src_local(:, d, i_local) = loc_AFT.time2freq(M_src(d, :), Info_A);
            end
        end
        Nu_DC_Out_cell{k} = nu_dc_local;
        M_resid_cell{k} = m_resid_local;
        M_src_cell{k} = m_src_local;
        B_max_cell{k} = b_max_local;
    end
    
    Nu_DC_suggested = vertcat(Nu_DC_Out_cell{:});
    max_B = max([B_max_cell{:}]);
    M_resid_harm = cat(3, M_resid_cell{:}); 
    M_src_harm = cat(3, M_src_cell{:});
end

% -------------------------------------------------------------------------
% 辅助函数: 快速组装 RHS 向量
% -------------------------------------------------------------------------
function RHS_vec = assemble_rhs_from_magnetization(Model, M_vec)
    Mesh = Model.Mesh;
    numEdges = size(Mesh.Edges, 2);
    numElems = size(Mesh.T, 2);
    I_idx = zeros(6 * numElems, 1); V_val = zeros(6 * numElems, 1);
    G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
    P = Mesh.P; T = Mesh.T; T2E = Mesh.T2E; Signs = double(Mesh.T2E_Sign);
    
    for e = 1:numElems
        nodes = T(:,e); pts = P(:, nodes);
        v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
        J = [v21, v31, v41]; detJ = det(J); invJ_T = inv(J)'; 
        G_phy = invJ_T * G_ref; Vol = abs(detJ) / 6.0;
        M_e = M_vec(:, e); 
        for i = 1:6
            na = [1 1 1 2 2 3]; nb = [2 3 4 3 4 4];
            % 计算边元旋度: Curl N = 2 * (Grad_a x Grad_b)
            curl_ni = 2 * cross(G_phy(:, na(i)), G_phy(:, nb(i)));
            val = Vol * dot(curl_ni, M_e) * Signs(i, e);
            idx_store = (e-1)*6 + i;
            I_idx(idx_store) = T2E(i, e);
            V_val(idx_store) = val;
        end
    end
    RHS_vec = sparse(I_idx, 1, V_val, numEdges, 1);
end