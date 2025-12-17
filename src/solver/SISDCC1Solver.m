classdef SISDCC1Solver < handle
    % SISDCC1SOLVER 支持 C1 连续预测的 SISDC 求解器 (v13.0 - C1 Predictor)
    % 
    % 继承自 SISDCSolver 架构，改进了时间片连接处的连续性。
    % 核心改进:
    %   1. C1 Predictor: 利用上一时间步的导数信息进行线性外推作为初值猜测。
    %      消除时间片连接处的 "Kinks" (折角)，显著减少启动时的迭代误差。
    %   2. Spectral Differentiation: 利用谱微分矩阵高精度计算末端导数。
    
    properties
        Assembler       % Finite Element Assembler
        LinearSolver    % Linear Solver Wrapper
        
        PolyOrder = 3        % 多项式阶数 P
        MaxSDCIters = 5      % SDC 扫描次数
        SDCTolerance = 1e-4  % SDC 收敛容差
    end
    
    methods
        function obj = SISDCC1Solver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 300; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            % SOLVE 主求解函数 (带 C1 预测逻辑)
            
            fprintf('==================================================\n');
            fprintf('   SDC Solver (C1 Continuous Predictor, P=%d)\n', obj.PolyOrder);
            fprintf('==================================================\n');
            
            % --- 1. 初始化与预计算 ---
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
            
            % --- 2. 谱离散参数计算 ---
            fprintf('   [Init] Computing GLL nodes, integration & differentiation matrices...\n');
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            % [新增] 计算谱微分矩阵 D (用于计算导数)
            D_mat = obj.computeDifferentiationMatrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.D = D_mat; % 存入 Context
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            NumPlotPoints = 20;
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
            
            % [新增] 上一时间步终点的物理导数 (用于 C1 预测)
            u_dot_prev = zeros(ctx.numTotalDofs, 1);
            
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
                
                % --- [改进] 1. C1 Continuous Initialization ---
                if slab_idx == 1
                    % 第一步无法利用历史导数，使用常数预测
                    X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                else
                    % 后续步：利用 u_dot_prev 进行线性外推
                    % t_local 在 [0, dt_slab] 之间
                    t_local_nodes = (ctx.Spectral.Nodes + 1) / 2 * dt_slab; 
                    
                    % Prediction: x(t) = x_0 + v_0 * t
                    % x_start: [Dofs x 1], u_dot_prev: [Dofs x 1], t_local: [Nodes x 1]
                    % Result: [Dofs x Nodes]
                    X_slab = x_start + u_dot_prev * t_local_nodes';
                    
                    % fprintf('   [Predictor] Applied C1 extrapolation.\n');
                end
                
                % --- 2. SDC Sweep ---
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
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
                
                % --- 3. 计算末端导数 (为下一步做准备) ---
                % --- 修改为：鲁棒的有限差分导数 (Robust Finite Difference) ---
                % 利用上一 Slab 最后两个点 (u_N 和 u_{N-1}) 计算斜率
                % 这种方法对高频数值噪声不敏感，预测更保守、更安全。
                
                % 1. 获取倒数两个节点的解
                u_final = X_slab(:, end);
                u_penultimate = X_slab(:, end-1);
                
                % 2. 获取倒数两个节点对应的物理时间
                t_nodes_phys = (ctx.Spectral.Nodes + 1)/2 * dt_slab;
                dt_last_segment = t_nodes_phys(end) - t_nodes_phys(end-1);
                
                % 3. 计算有限差分斜率
                u_dot_prev = (u_final - u_penultimate) / dt_last_segment;
                
                % [可选] 4. 增加阻尼限制 (Damping/Limiter)
                % 如果预测导致 B 场变化太剧烈，可能会发散。
                % 简单的工程技巧：稍微衰减导数 (例如 0.8 倍)，让预测值 "Under-shoot" 而不是 "Over-shoot"
                u_dot_prev = u_dot_prev * 0.8;
                
                % --- [Data Collection] ---
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
                
                idx_start = 1; if slab_idx > 1, idx_start = 2; end
                full_time_list = [full_time_list; t_plot_phys(idx_start:end)];
                full_current_list = [full_current_list; I_plot(idx_start:end)];
                full_B_list = [full_B_list; B_plot(idx_start:end)];
                
                % Update State
                x_start = X_slab(:, end); 
                globalTime = t_end;
                currentHistory(slab_idx) = x_start(end);
                if ~isempty(postProc), probeB_History(slab_idx) = B_nodes(end); end
                
                if ~isempty(postProc)
                    fprintf('      [End] |B|=%.4f T, I=%.4f A\n', B_nodes(end), full_current_list(end));
                end
                
                if ~isempty(monitorHandle)
                    try, monitorHandle(globalTime, x_start(end), full_time_list, full_current_list); drawnow limitrate; catch; end
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
        
        function D = computeDifferentiationMatrix(~, nodes)
            % COMPUTEDIFFERENTIATIONMATRIX 计算 GLL 节点的谱微分矩阵 D
            % D_ij = dL_j(x_i) / dx
            
            N = length(nodes) - 1;
            D = zeros(N+1, N+1);
            
            % 获取每个节点处的 Legendre 多项式值 L_N(x)
            L_vals = zeros(N+1, 1);
            for i = 1:N+1
                [L_vals(i), ~] = SpectralTimeUtils.legendre_poly(nodes(i), N);
            end
            
            for i = 1:N+1
                for j = 1:N+1
                    if i ~= j
                        % 非对角元公式
                        D(i, j) = (L_vals(i) / L_vals(j)) / (nodes(i) - nodes(j));
                    else
                        % 对角元公式
                        if i == 1 % x = -1
                            D(i, j) = -0.25 * N * (N + 1);
                        elseif i == N+1 % x = 1
                            D(i, j) = 0.25 * N * (N + 1);
                        else
                            D(i, j) = 0.0;
                        end
                    end
                end
            end
        end
        
        function [X_new, max_diff_norm, final_res_norm] = performSISDCSweep(obj, X_k, x0, ctx)
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            F_k = obj.evaluateF(X_k, ctx); 
            dt_half = dt_slab / 2.0;
            I_accum = (F_k * Q.') * dt_half;
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            for m = 2:num_nodes
                tau_curr = nodes(m);
                tau_prev = nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half; 
                
                u_prev_new = X_new(:, m-1);   
                f_curr_old = F_k(:, m);       
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                Scale = ctx.CircuitRowScale;
                M_sys = ctx.M_sigma; 
                [w_idx, ~, w_val] = find(ctx.W_vec);
                
                M_u_prev = M_sys * u_prev_new;
                I_prev = u_prev_new(end);
                term_L = ctx.circuit.L * I_prev;
                term_W = ctx.W_vec' * u_prev_new;
                M_u_prev(end) = M_u_prev(end) + (term_L + term_W) * Scale;
                
                RHS_Vector = M_u_prev + int_delta - dt_sub * f_curr_old;
                
                dist_history = norm(X_k(:, m) - X_k(:, m-1));
                if dist_history < 1e-9 
                    u_guess = u_prev_new;
                else
                    u_guess = X_k(:, m);
                end
                
                [u_next, res_norm, iter_count] = obj.solveImplicitStep(u_guess, dt_sub, RHS_Vector, nodes(m), ctx);
                
                X_new(:, m) = u_next;
                
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
                
                fprintf('        -> Node %d Solved: %d Iters, Res=%.2e\n', m, iter_count, res_norm);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            
            C_tot = ctx.C_AP + ctx.C_AP';
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            
            for m = 1:num_nodes
                x = X_mat(:, m);
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half; 
                
                A_vec = x(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(R_fem)) = -R_fem; 
                f_vec = f_vec - C_tot * x;
                f_vec = f_vec + ctx.W_vec * x(end); 
                
                V_src = ctx.circuit.V_source_func(t);
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                
                F_val(:, m) = f_vec;
            end
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_node, ctx)
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_tot = M_sys + M_cir_diag + M_row_coup;
            
            final_res = 1e10;
            res_norm_init = 1.0; 
            res_norm_prev = 1e10;
            iter_count = 0;
            
            REL_TOL = 1e-6;      
            STAG_TOL = 1e-3;     
            MAX_ITERS = 15;
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                
                A_vec = u_curr(1:ctx.num_A);
                [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                
                I_val = u_curr(end);
                neg_Force = neg_Force - ctx.W_vec * I_val; 
                
                t_phys = ctx.t_start + (t_node + 1) * ctx.dt_slab / 2.0; 
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
                
                Res = M_tot * u_curr + dt * neg_Force - RHS_Vector;
                res_norm = norm(Res);
                
                if iter == 1
                    res_norm_init = res_norm;
                    if res_norm_init < 1e-20, res_norm_init = 1.0; end 
                end
                
                rel_err = res_norm / res_norm_init;
                
                if rel_err < REL_TOL
                    final_res = res_norm;
                    fprintf('          Newton Converged: RelErr=%.2e\n', rel_err);
                    break;
                end
                
                if res_norm < 1e-10
                    final_res = res_norm;
                    break;
                end
                
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / res_norm_prev;
                    if abs(improvement) < STAG_TOL
                        final_res = res_norm;
                        break;
                    end
                end
                
                res_norm_prev = res_norm;
                
                [ki, kj, kv] = find(K_fem_small);
                K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                
                K_base = K_fem + ctx.C_AP + ctx.C_AP';
                K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
                K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
                
                K_tan_tot = K_base + K_col_coup + K_cir_diag;
                Jac = M_tot + dt * K_tan_tot;
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(Jac, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                
                u_curr = u_curr + du;
            end
            
            u_new = u_curr;
            
            if final_res / res_norm_init > 1e-1 && final_res > 1.0
                fprintf('          [Warn] Local Newton stalled at high residual: Rel=%.2e\n', final_res / res_norm_init);
            end
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            
            k_soft_rep = k_hard_median / 1000; 
            
            k_target = sqrt(k_hard_median * k_soft_rep);
            Z_cir = abs(ctx.circuit.R + ctx.circuit.L/dt); 
            max_W = max(abs(nonzeros(ctx.W_vec)));         
            max_coupling_val = max_W / dt;               
            max_circuit_val = max(Z_cir, max_coupling_val);
            if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
            scaleFactor = k_target / max_circuit_val;
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            meshTags = obj.Assembler.Mesh.RegionTags;
            numElems = length(meshTags);
            nu_vec = zeros(numElems, 1);
            matLib = ctx.matLib;
            nu0 = 1 / (4 * pi * 1e-7);
            uniqueTags = unique(meshTags);
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k);
                if ~matLib.isKey(tag), continue; end
                mat = matLib(tag);
                val = 0;
                if strcmp(mat.Type, 'Linear')
                    val = mat.Nu_Linear;
                else
                    if strcmp(mode, 'hard')
                        val = nu0; 
                    else
                        val = nu0 / 1000.0; 
                    end
                end
                nu_vec(meshTags == tag) = val;
            end
        end
        
        function vals_interp = interpolatePolynomial(~, nodes, values, target_nodes)
            N = length(nodes);
            M = length(target_nodes);
            vals_interp = zeros(M, 1);
            for k = 1:M
                t = target_nodes(k);
                y = 0;
                for j = 1:N
                    l_j = 1;
                    for i = 1:N
                        if i ~= j
                            l_j = l_j * (t - nodes(i)) / (nodes(j) - nodes(i));
                        end
                    end
                    y = y + values(j) * l_j;
                end
                vals_interp(k) = y;
            end
        end
    end
end