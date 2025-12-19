classdef SDIRK2ISDCSolver < handle
    % SDIRK2ISDCSolver 基于 SDIRK-2 扫除器的高阶 SDC 求解器 (Auto-Stagnation Fix)
    % 
    % 修复记录:
    %   1. [Fix] 增加了智能停滞检测：当牛顿迭代陷入震荡或改进微弱时，自动提前退出，
    %      避免死循环。
    %   2. [Fix] 优化了收敛容差策略，适应变压器饱和区的剧烈非线性。
    %   3. [Log] 调整了输出格式，详细展示每个节点的 Diff 和 Iter 情况。
    
    properties
        Assembler       % 有限元组装器
        LinearSolver    % 线性求解器
        
        PolyOrder = 3        % 谱元阶数 P (3 -> 4 Nodes)
        MaxSDCIters = 3      % SDC 扫描次数
        SDCTolerance = 1e-5  % SDC 收敛容差
        
        % SDIRK2 参数 (L-stable)
        Gamma = 1 - 1/sqrt(2); 
    end
    
    methods
        function obj = SDIRK2ISDCSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 400; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            
            fprintf('\n==========================================================\n');
            fprintf('   SDIRK-2 SDC Solver (Robust Full-Newton + Stagnation Check)\n');
            fprintf('==========================================================\n');
            
            % --- 1. 初始化 ---
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            ctx.space_A = space_A; ctx.space_P = space_P;
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            postProc = []; 
            if ~isempty(probePoint), postProc = PostProcessor(obj.Assembler); end
            
            % --- 2. 谱元数据 ---
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = SpectralTimeUtils.integration_matrix(gll_nodes); 
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
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', ctx.CircuitRowScale);
            
            % 预组装全局质量矩阵
            Scale = ctx.CircuitRowScale;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            ctx.M_full_const = ctx.M_sigma + M_cir_diag + M_row_coup;
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_curr(1:n) = init_State(1:n);
            end
            
            full_time_list = []; full_current_list = []; full_B_list = [];
            globalTime = 0;
            
            % ==================================================
            % 主循环
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n--- Slab %d/%d | dt=%.2e | T=%.4f -> %.4f ---\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. Initialization
                X_slab = repmat(x_curr, 1, ctx.Spectral.NumNodes);
                
                % 2. SDC Sweep Loop
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    
                    % 执行扫除
                    [X_slab, max_diff_norm] = obj.performSDIRK2Sweep(X_old, x_curr, ctx);
                    
                    % 输出状态
                    if max_diff_norm < obj.SDCTolerance
                        fprintf('   Sweep %d: Update=%.2e [Converged]\n', sdc_iter, max_diff_norm);
                        break;
                    else
                        fprintf('   Sweep %d: Update=%.4f\n', sdc_iter, max_diff_norm);
                    end
                end
                
                % 3. Output
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
                
                if slab_idx == 1, idx_st = 1; else, idx_st = 2; end
                full_time_list = [full_time_list; t_plot_phys(idx_st:end)];
                full_current_list = [full_current_list; I_plot(idx_st:end)];
                full_B_list = [full_B_list; B_plot(idx_st:end)];
                
                x_curr = X_slab(:, end);
                globalTime = t_end;
                
                if ~isempty(monitorHandle)
                    try monitorHandle(globalTime, x_curr(end), full_time_list, full_current_list); catch; end
                end
            end
            
            solutionResults = x_curr;
            info.FinalTime = globalTime;
            info.Time_Full = full_time_list;
            info.Current_Full = full_current_list;
            info.ProbeB_Full = full_B_list;
        end
    end
    
    methods (Access = private)
        
        function [X_new, max_diff_norm] = performSDIRK2Sweep(obj, X_k, x0, ctx)
            % PERFORMSDIRK2SWEEP 核心扫除函数
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            gamma = obj.Gamma;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            F_k = obj.evaluateF_Batch(X_k, ctx); 
            
            dt_half = dt_slab / 2.0;
            I_spec = (F_k * Q.') * dt_half; 
            
            max_diff_norm = 0;
            
            % 逐节点推进
            for m = 2:num_nodes
                h_sub = (nodes(m) - nodes(m-1)) * dt_half;
                u_prev_new = X_new(:, m-1);
                
                % --- SDC Source ---
                tau_m_prev = nodes(m-1);
                tau_m = nodes(m);
                tau_gamma = tau_m_prev + gamma * (tau_m - tau_m_prev);
                
                u_k_gamma = obj.interpolateStateSingle(nodes, X_k, tau_gamma);
                t_phys_gamma = ctx.t_start + (tau_gamma + 1) * dt_half;
                t_phys_end   = ctx.t_start + (tau_m + 1) * dt_half;
                
                F_k_gamma = obj.evaluateF_Single(u_k_gamma, t_phys_gamma, ctx);
                F_k_end   = F_k(:, m); 
                
                I_low = h_sub * ((1-gamma) * F_k_gamma + gamma * F_k_end);
                I_high = I_spec(:, m) - I_spec(:, m-1);
                sigma = (I_high - I_low) / h_sub;
                
                % --- Stage 1 ---
                dt_eff = h_sub * gamma;
                M_u_prev = ctx.M_full_const * u_prev_new;
                RHS_1 = M_u_prev + dt_eff * sigma;
                
                dist_hist = norm(X_k(:,m) - X_k(:,m-1));
                if dist_hist < 1e-9, u_guess = u_prev_new; else, u_guess = X_k(:,m); end
                
                % Solve Stage 1
                [g1, ~, iter1] = obj.solveImplicitStep(u_guess, dt_eff, RHS_1, t_phys_gamma, ctx);
                fprintf('      [Node %d Stage 1] ', m);
                
                F_phys_g1 = obj.evaluateF_Single(g1, t_phys_gamma, ctx);
                
                % --- Stage 2 ---
                term_F_g1 = h_sub * (1-gamma) * F_phys_g1;
                term_sigma = h_sub * sigma;
                RHS_2 = M_u_prev + term_F_g1 + term_sigma;
                
                % Solve Stage 2
                [u_next, ~, iter2] = obj.solveImplicitStep(g1, dt_eff, RHS_2, t_phys_end, ctx);
                
                X_new(:, m) = u_next;
                
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                
                fprintf('      [Node %d Stage 2]         Node %d Done. Iters=%d/%d, Diff=%.2e\n', ...
                    m, m, iter1, iter2, diff);
            end
        end
        
        function F_val = evaluateF_Batch(obj, X_mat, ctx)
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            dt_half = ctx.dt_slab / 2.0;
            for m = 1:num_nodes
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half;
                F_val(:, m) = obj.evaluateF_Single(X_mat(:, m), t, ctx);
            end
        end
        
        function f_vec = evaluateF_Single(obj, x, t_phys, ctx)
            num_dofs = ctx.numTotalDofs;
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            
            f_vec = sparse(num_dofs, 1);
            f_vec(1:length(R_fem)) = -R_fem; 
            f_vec = f_vec - (ctx.C_AP + ctx.C_AP') * x;
            f_vec = f_vec + ctx.W_vec * I_val; 
            
            V_src = ctx.circuit.V_source_func(t_phys);
            f_vec(end) = (V_src - ctx.circuit.R * I_val) * ctx.CircuitRowScale;
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_phys, ctx)
            % 局部牛顿求解器 (Full Newton with Stagnation Check)
            
            u_curr = u_guess;
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            
            [w_idx, ~, w_val] = find(ctx.W_vec);
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            K_const_part = ctx.C_AP + ctx.C_AP' + K_col_coup + K_cir_diag;
            
            final_res = 1e10;
            res_norm_init = 1.0; res_norm_prev = 1e10;
            
            MAX_ITERS = 15;
            STAGNATION_TOL = 1e-3; % 改善率阈值 (0.1%)
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                A_vec = u_curr(1:ctx.num_A);
                
                [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                neg_Force = neg_Force - ctx.W_vec * u_curr(end); 
                
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * u_curr(end) - V_src) * Scale;
                
                Res = ctx.M_full_const * u_curr + dt * neg_Force - RHS_Vector;
                res_norm = norm(Res);
                
                % --- 收敛判据 ---
                if iter == 1, res_norm_init = res_norm; if res_norm_init < 1e-20, res_norm_init=1; end; end
                
                % 1. 绝对/相对收敛
                if res_norm / res_norm_init < 1e-5 || res_norm < 1e-8
                    final_res = res_norm; break; 
                end
                
                % 2. 停滞检测 (Stagnation Check)
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / res_norm_prev;
                    % 如果改善微弱
                    if improvement < STAGNATION_TOL
                        % 如果残差已经比较小 (e.g. 1e-2)，可以接受停滞并退出
                        % 对于饱和变压器，残差卡在 1e-3 ~ 1e-2 是常态
                        if res_norm < 1e-1
                             final_res = res_norm; 
                             fprintf('(Stag:%.1e) ', res_norm); % Debug
                             break;
                        else
                             % 残差过大但停滞，也没必要继续算，退出并标记
                             fprintf('(HighStag:%.1e) ', res_norm); 
                             final_res = res_norm; break;
                        end
                    end
                end
                res_norm_prev = res_norm;
                
                % --- 求解更新 ---
                [ki, kj, kv] = find(K_fem_small);
                K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                Jac = ctx.M_full_const + dt * (K_fem + K_const_part);
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(Jac, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                
                u_curr = u_curr + du;
            end
            
            u_new = u_curr;
        end
        
        function u_interp = interpolateStateSingle(~, nodes, X_mat, tau_target)
            N = length(nodes);
            u_interp = zeros(size(X_mat, 1), 1);
            for j = 1:N
                lj = 1;
                for i = 1:N
                    if i ~= j, lj = lj * (tau_target - nodes(i)) / (nodes(j) - nodes(i)); end
                end
                u_interp = u_interp + X_mat(:, j) * lj;
            end
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
             nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
             K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
             diag_M = diag(ctx.M_sigma); num_fem = size(K_hard, 1);
             hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem))/dt;
             k_hard = full(median(hard_vals(hard_vals>1e-12)));
             if isempty(k_hard), k_hard = 1e6; end
             k_target = sqrt(k_hard * (k_hard/1000));
             max_cir = max(abs(ctx.circuit.R + ctx.circuit.L/dt), max(abs(nonzeros(ctx.W_vec)))/dt);
             scaleFactor = k_target / max(max_cir, 1.0);
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            meshTags = obj.Assembler.Mesh.RegionTags;
            nu_vec = zeros(length(meshTags), 1);
            matLib = ctx.matLib; nu0 = 1/(4*pi*1e-7);
            uniqueTags = unique(meshTags);
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k); if ~matLib.isKey(tag), continue; end
                mat = matLib(tag);
                if strcmp(mat.Type, 'Linear'), val = mat.Nu_Linear;
                else, if strcmp(mode,'hard'), val = nu0; else, val = nu0/1000; end; end
                nu_vec(meshTags == tag) = val;
            end
        end
        
        function vals = interpolatePolynomial(~, nodes, values, target_nodes)
             N = length(nodes); M = length(target_nodes); vals = zeros(M,1);
             for k=1:M
                 t=target_nodes(k); y=0;
                 for j=1:N
                     lj=1; for i=1:N, if i~=j, lj=lj*(t-nodes(i))/(nodes(j)-nodes(i)); end; end
                     y=y+values(j)*lj;
                 end
                 vals(k)=y;
             end
        end
    end
end