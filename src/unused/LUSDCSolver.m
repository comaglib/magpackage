classdef LUSDCSolver < handle
    % LUSDCSolver 基于非线性序列扫描的 SDC 求解器 (v12.5 - Force Recovery Optimized)
    % 
    % 核心机制: Semi-Implicit SDC (SISDC) + GLL Nodes
    %   1. 节点策略: Gauss-Lobatto-Legendre (GLL) 节点，包含区间两端，插值更自然。
    %   2. 扫描策略: 序列扫描 (Sequential Sweeping)，利用最新的前序节点信息。
    %   3. [优化] Force Recovery: 利用牛顿迭代后的残差直接反推物理力，避免额外组装。
    
    properties
        Assembler       % Finite Element Assembler (有限元组装器)
        LinearSolver    % Linear Solver Wrapper (线性求解器)
        
        PolyOrder = 3        % 多项式阶数 P (建议 3 或 4)。
        MaxSDCIters = 5      % SDC 扫描次数 K (通常 5 次足够)。
        SDCTolerance = 1e-4  % SDC 收敛容差。
    end
    
    methods
        function obj = LUSDCSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; % 非对称
            obj.LinearSolver.MumpsICNTL.i14 = 300; % 增加工作内存
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            % SOLVE 主求解函数
            
            fprintf('==================================================\n');
            fprintf('   SDC Solver (Nonlinear SISDC, P=%d, Force Recovery)\n', obj.PolyOrder);
            fprintf('==================================================\n');
            
            % --- 1. 初始化 (Initialization) ---
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            ctx.space_A = space_A; ctx.space_P = space_P;
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            postProc = []; 
            if ~isempty(probePoint), postProc = PostProcessor(obj.Assembler); end
            
            full_time_list = []; full_current_list = []; full_B_list = [];
            probeB_History = zeros(length(timeSteps), 1);
            currentHistory = zeros(length(timeSteps), 1);
            
            % --- 2. GLL 节点与积分矩阵 ---
            fprintf('   [Init] Computing GLL nodes and integration matrix...\n');
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            NumPlotPoints = 50;
            tau_plot = linspace(-1, 1, NumPlotPoints)'; 
            
            % --- 3. 组装不变矩阵 ---
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', ctx.CircuitRowScale);
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
            globalTime = 0;
            
            % ==================================================
            % 主时间循环
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n=== Slab %d/%d (dt=%.1e, T=%.4f -> %.4f) ===\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % Initialization: Constant Prediction
                X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                
                % SDC Sweep
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    % 调用优化后的扫描函数
                    [X_slab, max_diff_norm, residual_norm] = obj.performSISDCSweep(X_old, x_start, ctx);
                    
                    if max_diff_norm < 1e-4
                        fprintf('   Iter %d: Update=%.2e, Res=%.2e', sdc_iter, max_diff_norm, residual_norm);
                    else
                        fprintf('   Iter %d: Update=%.4f, Res=%.4f', sdc_iter, max_diff_norm, residual_norm);
                    end
                    
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' -> Converged.\n');
                        break;
                    else
                        fprintf('\n');
                    end
                end
                
                % --- Data Collection ---
                t_plot_phys = t_start + (tau_plot + 1) / 2 * dt_slab;
                I_nodes = X_slab(end, :)';
                I_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, I_nodes, tau_plot);
                
                B_nodes = zeros(ctx.Spectral.NumNodes, 1);
                if ~isempty(postProc)
                    for k = 1:ctx.Spectral.NumNodes
                        B_val = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_nodes(k) = norm(B_val);
                    end
                end
                B_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, B_nodes, tau_plot);
                
                if slab_idx == 1, idx_start = 1; else, idx_start = 2; end
                
                full_time_list = [full_time_list; t_plot_phys(idx_start:end)];
                full_current_list = [full_current_list; I_plot(idx_start:end)];
                full_B_list = [full_B_list; B_plot(idx_start:end)];
                
                x_start = X_slab(:, end); 
                globalTime = t_end;
                currentHistory(slab_idx) = x_start(end);
                if ~isempty(postProc), probeB_History(slab_idx) = B_nodes(end); end
                
                if ~isempty(postProc)
                    fprintf('      [End] |B|=%.4f T, I=%.4f A\n', B_nodes(end), full_current_list(end));
                end
                
                if ~isempty(monitorHandle)
                    try
                        monitorHandle(globalTime, x_start(end), full_time_list, full_current_list);
                        drawnow limitrate;
                    catch; end
                end
            end
            
            solutionResults = x_start;
            info.FinalTime = globalTime;
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc), info.ProbeB_History = probeB_History; end
            info.Time_Full = full_time_list;
            info.Current_Full = full_current_list;
            info.ProbeB_Full = full_B_list;
        end
    end
    
    methods (Access = private)
        
        function [X_new, max_diff_norm, final_res_norm] = performSISDCSweep(obj, X_k, x0, ctx)
            % PERFORMSISDCSWEEP (Force Recovery Optimized)
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            % 1. 计算上一轮解的物理力 F(X_k) (Batch Mode)
            F_k = obj.evaluateF(X_k, ctx); 
            
            % 2. 计算积分项 I_accum
            dt_half = dt_slab / 2.0;
            I_accum = (F_k * Q.') * dt_half;
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            % --- [Optimization] 预组装常数质量矩阵 ---
            % 用于 Force Recovery: F = (M*u - RHS)/dt
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma; 
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_full_const = M_sys + M_cir_diag + M_row_coup;
            
            % 3. 序列扫描 (m = 2 to N)
            for m = 2:num_nodes
                tau_curr = nodes(m);
                tau_prev = nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half;
                
                u_prev_new = X_new(:, m-1);   
                f_curr_old = F_k(:, m);       
                
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                % RHS Construction
                M_u_prev = M_sys * u_prev_new;
                I_prev = u_prev_new(end);
                term_L = ctx.circuit.L * I_prev;
                term_W = ctx.W_vec' * u_prev_new;
                M_u_prev(end) = M_u_prev(end) + (term_L + term_W) * Scale;
                
                RHS_Vector = M_u_prev + int_delta - dt_sub * f_curr_old;
                
                % Smart Guess
                dist_history = norm(X_k(:, m) - X_k(:, m-1));
                if dist_history < 1e-9 
                    u_guess = u_prev_new; % Rolling
                else
                    u_guess = X_k(:, m); % Refining
                end
                
                % Implicit Solve
                [u_next, res_norm, iter_count] = obj.solveImplicitStep(u_guess, dt_sub, RHS_Vector, nodes(m), ctx);
                
                X_new(:, m) = u_next;
                
                % =========================================================
                % [Optimization] Force Recovery (力恢复)
                % 替代昂贵的 F_new_col = obj.evaluateF(u_next, ctx, m);
                % 原理: M*u - dt*F = RHS => F = (M*u - RHS) / dt
                % =========================================================
                if abs(dt_sub) > 1e-12
                    % 使用预组装的 M_full_const 进行矩阵向量乘法 (极快)
                    F_recovered = (M_full_const * u_next - RHS_Vector) / dt_sub;
                    
                    % 注意: 如果 evaluateF 内部包含时变项(如时变 L)，这里需额外处理。
                    % 但当前实现中 M 是常数，非线性只在 F(K*u) 中，所以该公式严格成立。
                    
                    % 更新 F_k 的当前列，供下一个节点使用 (SISDC 需要 F_m)
                    % 注意: SISDC 依赖的是"最新算出的力"，所以我们直接覆盖 F_k(:, m)
                    F_k(:, m) = F_recovered;
                else
                    % 保护: 极小步长回退到标准计算
                    F_k(:, m) = obj.evaluateF(u_next, ctx, m);
                end
                % =========================================================
                
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
                
                fprintf('        -> Node %d Solved: %d Iters, Res=%.2e\n', m, iter_count, res_norm);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx, single_col_idx)
            % EVALUATEF 计算物理方程 F(u)
            % 支持 Batch 模式 (X_mat 为矩阵) 和 Single Column 模式 (需提供 single_col_idx)
            
            if nargin < 4, single_col_idx = []; end
            
            num_cols = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_cols);
            
            C_tot = ctx.C_AP + ctx.C_AP';
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            nodes = ctx.Spectral.Nodes;
            
            for k = 1:num_cols
                x = X_mat(:, k);
                
                % 确定当前列对应的 Radau 节点索引
                if isempty(single_col_idx)
                    node_idx = k; 
                else
                    node_idx = single_col_idx; 
                end
                
                % 计算精确物理时间
                t_phys = ctx.t_start + (nodes(node_idx) + 1) * dt_half;
                
                A_vec = x(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(R_fem)) = -R_fem; 
                f_vec = f_vec - C_tot * x;
                f_vec = f_vec + ctx.W_vec * x(end); 
                
                V_src = ctx.circuit.V_source_func(t_phys);
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                
                F_val(:, k) = f_vec;
            end
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_node, ctx)
            % SOLVEIMPLICITSTEP 局部牛顿求解器 (集成修正牛顿法 Modified Newton)
            % 优化策略: 仅在收敛变慢时更新雅可比矩阵 (J)，否则只更新残差 (R)。
            
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            % 预组装常数矩阵 (避免在循环内重复创建 sparse)
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_tot = M_sys + M_cir_diag + M_row_coup;
            
            % 常数刚度部分 (电路电阻 & 耦合项的负值)
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            K_const_part = ctx.C_AP + ctx.C_AP' + K_col_coup + K_cir_diag;
            
            final_res = 1e10;
            res_norm_init = 1.0; 
            res_norm_prev = 1e10;
            
            % --- 控制参数 ---
            REL_TOL = 1e-6;      
            STAG_TOL = 1e-3;     
            MAX_ITERS = 15;
            
            % [优化] 雅可比矩阵冻结标志
            update_jacobian = true; % 第一步必须更新
            J_cache = [];           % 缓存的雅可比矩阵
            
            % 物理时间
            t_phys = ctx.t_start + (t_node + 1) * ctx.dt_slab / 2.0; 
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                
                % 1. 组装阶段
                % 根据 update_jacobian 决定是否组装刚度矩阵 K
                A_vec = u_curr(1:ctx.num_A);
                
                if update_jacobian
                    % [Full Mode] 计算 K 和 R
                    [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                else
                    % [Fast Mode] 仅计算 R (速度快很多)
                    [~, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                    K_fem_small = []; % 占位
                end
                
                % 2. 构建残差 (Residual Construction)
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                I_val = u_curr(end);
                neg_Force = neg_Force - ctx.W_vec * I_val; 
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
                
                Res = M_tot * u_curr + dt * neg_Force - RHS_Vector;
                res_norm = norm(Res);
                
                % --- 收敛判断 ---
                if iter == 1
                    res_norm_init = res_norm;
                    if res_norm_init < 1e-20, res_norm_init = 1.0; end 
                end
                
                if res_norm / res_norm_init < REL_TOL || res_norm < 1e-10
                    final_res = res_norm; break;
                end
                
                % 停滞检测
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / res_norm_prev;
                    if abs(improvement) < STAG_TOL
                        final_res = res_norm; break;
                    end
                end
                
                % [优化关键] 决定下一次是否需要更新 Jacobian
                % 如果收敛速度良好 (例如误差减少了 5 倍以上)，则下一次继续冻结 J
                convergence_rate = res_norm / res_norm_prev;
                res_norm_prev = res_norm;
                
                if convergence_rate > 0.2 % 收敛变慢 (Rate > 0.2)
                    force_update_next = true;
                    % 如果当前处于冻结状态且收敛很差，可能需要立即更新并重算(这里简化为下步更新)
                else
                    force_update_next = false;
                end

                % 3. 线性求解阶段
                if update_jacobian
                    % 只有在 flag 为 true 时才重新构建 J
                    [ki, kj, kv] = find(K_fem_small);
                    K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                    K_tan_tot = K_fem + K_const_part;
                    J_cache = M_tot + dt * K_tan_tot; % 更新缓存
                    
                    % 标记下一次尝试冻结 (除非本轮收敛太慢)
                    update_jacobian = false; 
                else
                    % fprintf('          [Opt] Skipping Jacobian Assembly (Iter %d)\n', iter);
                end
                
                % 如果上一步收敛太慢，强制下一次更新
                if force_update_next, update_jacobian = true; end
                
                % 求解 (使用缓存的 J_cache)
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_cache, Res, ctx.fixedDofs);
                
                % 注意: 如果 LinearSolver 支持重用因子 (Reuse Factors)，这里会有巨大提升。
                % 即使不支持，能省去上面的 sparse 组装和 J_cache 构建也是很大的优化。
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                
                u_curr = u_curr + du;
            end
            u_new = u_curr;
            
            if final_res / res_norm_init > 1e-1 && final_res > 1.0
                fprintf('          [Warn] Local Newton stalled: Rel=%.2e\n', final_res / res_norm_init);
            end
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            k_target = sqrt(k_hard_median * (k_hard_median / 1000));
            max_cir = max(abs(ctx.circuit.R + ctx.circuit.L/dt), max(abs(nonzeros(ctx.W_vec)))/dt);
            scaleFactor = k_target / max(max_cir, 1.0);
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            meshTags = obj.Assembler.Mesh.RegionTags;
            nu_vec = zeros(length(meshTags), 1);
            matLib = ctx.matLib; nu0 = 1 / (4 * pi * 1e-7);
            uniqueTags = unique(meshTags);
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k); if ~matLib.isKey(tag), continue; end
                mat = matLib(tag);
                if strcmp(mat.Type, 'Linear'), val = mat.Nu_Linear;
                else, if strcmp(mode, 'hard'), val = nu0; else, val = nu0 / 1000.0; end; end
                nu_vec(meshTags == tag) = val;
            end
        end
        
        function vals_interp = interpolatePolynomial(~, nodes, values, target_nodes)
            N = length(nodes); M = length(target_nodes); vals_interp = zeros(M, 1);
            for k = 1:M
                t = target_nodes(k); y = 0;
                for j = 1:N
                    l_j = 1;
                    for i = 1:N
                        if i ~= j, l_j = l_j * (t - nodes(i)) / (nodes(j) - nodes(i)); end
                    end
                    y = y + values(j) * l_j;
                end
                vals_interp(k) = y;
            end
        end
    end
end