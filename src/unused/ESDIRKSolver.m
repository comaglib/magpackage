classdef ESDIRKSolver < handle
    % ESDIRKSolver 基于 ESDIRK 的刚性问题求解器 (v1.2 - Robust Divergence Guard)
    % 
    % 核心机制: ESDIRK (Explicit Singly Diagonally Implicit Runge-Kutta)
    % 算法: TR-BDF2 (3级, 2阶, L-stable)
    % 
    % [修复更新 v1.2]
    %   1. 修复了 Jacobian 永久冻结的 Bug。
    %   2. 增加了 [Divergence Guard]: 当检测到残差反弹时，自动回滚(Rollback)并强制更新 Jacobian。
    %   3. 增加了 [Slow Convergence Check]: 收敛过慢时自动更新 Jacobian。
    
    properties
        Assembler       % 有限元组装器
        LinearSolver    % 线性求解器
        
        % 牛顿迭代相关参数
        NewtonTolerance = 1e-6 
        MaxNewtonIters = 20 % 稍微增加最大迭代次数以允许 Retry
    end
    
    methods
        function obj = ESDIRKSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 400; % 进一步增加内存预留
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
            
            fprintf('==========================================================\n');
            fprintf('   ESDIRK Solver (TR-BDF2) | Robust Mode\n');
            fprintf('==========================================================\n');
            
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
            
            % 数据容器
            full_time_list = []; full_current_list = []; full_B_list = [];
            probeB_History = zeros(length(timeSteps), 1);
            currentHistory = zeros(length(timeSteps), 1);
            
            % --- 2. 定义 TR-BDF2 Butcher 表 ---
            gamma = 2 - sqrt(2);
            RK.c = [0; gamma; 1];
            w = 1 / (2 * (2 - gamma));
            d = (1 - gamma) / (2 - gamma); 
            
            RK.A = [0,       0,       0; 
                    gamma/2, gamma/2, 0; 
                    w,       w,       d]; 
            RK.b = RK.A(3, :);
            RK.s = 3;
            
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
            
            % 预组装总质量矩阵
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
            
            globalTime = 0;
            
            % ==================================================
            % 主时间循环
            % ==================================================
            for step_idx = 1:length(timeSteps)
                dt = timeSteps(step_idx);
                t_n = globalTime;
                t_next = t_n + dt;
                
                fprintf('\n----------------------------------------------------------\n');
                fprintf(' Step %d/%d | dt=%.2e | T=%.5f -> %.5f\n', ...
                    step_idx, length(timeSteps), dt, t_n, t_next);
                fprintf('----------------------------------------------------------\n');
                
                u_stages = zeros(ctx.numTotalDofs, RK.s);
                F_stages = zeros(ctx.numTotalDofs, RK.s);
                
                % Jacobian Context for SDIRK (Modified Newton)
                jaco_ctx.is_cached = false;
                jaco_ctx.J_factor = [];
                jaco_ctx.prev_dt_gamma = 0;
                
                for i = 1:RK.s
                    t_stage = t_n + RK.c(i) * dt;
                    
                    % 1. Construct RHS
                    if i == 1, M_un = ctx.M_full_const * x_curr; end
                    
                    accum_F = sparse(ctx.numTotalDofs, 1);
                    for j = 1 : i-1
                        if RK.A(i, j) ~= 0
                            accum_F = accum_F + RK.A(i, j) * F_stages(:, j);
                        end
                    end
                    RHS_explicit = M_un + dt * accum_F;
                    
                    % 2. Solve Stage
                    a_ii = RK.A(i, i);
                    
                    if a_ii == 0
                        % Explicit Stage
                        fprintf('   [Stage %d] Explicit Prediction\n', i);
                        u_i = x_curr; 
                        F_i = obj.evaluateF(u_i, t_stage, ctx);
                        
                    else
                        % Implicit Stage
                        dt_eff = dt * a_ii;
                        fprintf('   [Stage %d] Implicit Solve (dt_eff=%.2e)\n', i, dt_eff);
                        
                        if i > 1, u_guess = u_stages(:, i-1); else, u_guess = x_curr; end
                        
                        % Newton Solve with Verbose Logging
                        [u_i, ~, iter_count, jaco_ctx] = obj.solveImplicitStep(u_guess, dt_eff, RHS_explicit, t_stage, ctx, jaco_ctx);
                        
                        % Force Recovery
                        F_i = (ctx.M_full_const * u_i - RHS_explicit) / dt_eff;
                    end
                    
                    u_stages(:, i) = u_i;
                    F_stages(:, i) = F_i;
                end
                
                % Update
                x_curr = u_stages(:, end);
                globalTime = t_next;
                
                % Post-processing (Logging)
                I_val_next = x_curr(end);
                B_val_end = 0;
                if ~isempty(postProc)
                    B_vec = postProc.probeB(x_curr(1:ctx.num_A), probePoint);
                    B_val_end = norm(B_vec);
                end
                
                % Store History
                t_plot_pts = [t_n; t_next];
                I_plot = [u_stages(end, 1); I_val_next];
                B_plot = [B_val_end; B_val_end];
                
                if step_idx == 1, idx_st = 1; else, idx_st = 2; end
                full_time_list = [full_time_list; t_plot_pts(idx_st:end)];
                full_current_list = [full_current_list; I_plot(idx_st:end)];
                full_B_list = [full_B_list; B_plot(idx_st:end)];
                
                currentHistory(step_idx) = I_val_next;
                probeB_History(step_idx) = B_val_end;
                
                if ~isempty(postProc)
                    fprintf('   [Result] |B|=%.4f T, I=%.4f A\n', B_val_end, I_val_next);
                end
                
                if ~isempty(monitorHandle)
                    try
                        monitorHandle(globalTime, I_val_next, full_time_list, full_current_list);
                        drawnow limitrate;
                    catch; end
                end
            end
            
            solutionResults = x_curr;
            info.FinalTime = globalTime;
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc), info.ProbeB_History = probeB_History; end
            info.Time_Full = full_time_list;
            info.Current_Full = full_current_list;
            info.ProbeB_Full = full_B_list;
        end
    end
    
    methods (Access = private)
        
        function F_val = evaluateF(obj, x, t_phys, ctx)
            % Explicit Force Evaluation
            num_dofs = ctx.numTotalDofs;
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            
            f_vec = sparse(num_dofs, 1);
            f_vec(1:length(R_fem)) = -R_fem; 
            f_vec = f_vec - (ctx.C_AP + ctx.C_AP') * x;
            f_vec = f_vec + ctx.W_vec * I_val;
            
            V_src = ctx.circuit.V_source_func(t_phys);
            f_cir = (V_src - ctx.circuit.R * I_val) * ctx.CircuitRowScale;
            f_vec(end) = f_cir;
            F_val = f_vec;
        end
        
        function [u_new, final_res, iter_count, jaco_ctx] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_phys, ctx, jaco_ctx)
            % 局部牛顿求解器 (Robust Divergence Guard)
            
            u_curr = u_guess;
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            
            % Constant stiffness parts
            [w_idx, ~, w_val] = find(ctx.W_vec);
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            K_const_part = ctx.C_AP + ctx.C_AP' + K_col_coup + K_cir_diag;
            
            final_res = 1e10;
            res_norm_init = 1.0; 
            res_norm_prev = 1e10;
            
            % --- Rollback Buffer ---
            u_prev_iter_val = u_curr; % 用于回滚
            
            % Decide Jacobian Update
            update_jacobian = true; 
            if isfield(jaco_ctx, 'prev_dt_gamma') && abs(dt - jaco_ctx.prev_dt_gamma) < 1e-9 && ~isempty(jaco_ctx.J_factor)
                % 如果 dt 没变且有缓存，可以尝试先用冻结的
                 % update_jacobian = false; 
            end
            if isfield(jaco_ctx, 'prev_dt_gamma') && abs(dt - jaco_ctx.prev_dt_gamma) > 1e-9
                update_jacobian = true; 
            end
            
            force_update_next = false;
            
            for iter = 1:obj.MaxNewtonIters
                iter_count = iter;
                
                % 0. Check Force Update from prev iter
                if force_update_next
                    update_jacobian = true;
                    force_update_next = false;
                end
                
                A_vec = u_curr(1:ctx.num_A);
                
                % 1. Assemble
                if update_jacobian
                    [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                    jac_status = '[Update]';
                else
                    [~, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                    K_fem_small = []; 
                    jac_status = '[Frozen]';
                end
                
                % 2. Residual
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                I_val = u_curr(end);
                neg_Force = neg_Force - ctx.W_vec * I_val; 
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
                
                Res = ctx.M_full_const * u_curr + dt * neg_Force - RHS_Vector;
                res_norm = norm(Res);
                
                % --- Divergence Guard (回滚机制) ---
                % 如果不是第一步，且残差反弹严重(>1.5倍)，且当前是 Frozen 状态
                if iter > 1 && res_norm > res_norm_prev * 1.5 && ~update_jacobian
                    fprintf('      Iter %2d | Res=%.4e | [Divergence Detected] Rolling back & Forcing Update.\n', iter, res_norm);
                    
                    % 1. 回滚解向量
                    u_curr = u_prev_iter_val;
                    % 2. 恢复残差基准
                    res_norm = res_norm_prev; 
                    % 3. 强制下一次立即更新
                    update_jacobian = true;
                    force_update_next = false; 
                    
                    % 4. 立即进入下一次循环 (Retry with Update)
                    % 注意：此时 iter 会增加，这是对的，消耗一次迭代次数来修正
                    continue; 
                end
                
                % 保存当前解用于下一次可能的回滚
                u_prev_iter_val = u_curr;
                
                % 3. Check Convergence
                if iter == 1
                    res_norm_init = res_norm;
                    if res_norm_init < 1e-20, res_norm_init = 1.0; end 
                end
                rel_err = res_norm / res_norm_init;
                
                fprintf('      Iter %2d | Res=%.4e | Rel=%.4e %s', iter, res_norm, rel_err, jac_status);
                
                if rel_err < obj.NewtonTolerance || res_norm < 1e-10
                    final_res = res_norm;
                    fprintf(' -> Converged.\n');
                    break;
                end
                fprintf('\n'); 
                
                % Slow Convergence Check
                % 如果收敛速度慢 (Rate > 0.1)，下一次必须更新
                rate = res_norm / res_norm_prev;
                if iter > 1 && rate > 0.1
                    force_update_next = true;
                    % fprintf('      [Info] Slow convergence (Rate=%.2f). Force Update next.\n', rate);
                end
                
                res_norm_prev = res_norm;
                
                % 4. Solve Linear System
                if update_jacobian
                    [ki, kj, kv] = find(K_fem_small);
                    K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                    Jac = ctx.M_full_const + dt * (K_fem + K_const_part);
                    
                    jaco_ctx.J_factor = Jac; 
                    jaco_ctx.prev_dt_gamma = dt;
                    
                    update_jacobian = false; % Next one try frozen
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(jaco_ctx.J_factor, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                u_curr = u_curr + du;
            end
            
            u_new = u_curr;
            
            if final_res / res_norm_init > 1e-1 && final_res > 1.0
                fprintf('      [CRITICAL] Newton Failed! Final Res=%.2e\n', final_res);
            end
        end
        
        % ... (Helper functions keep same) ...
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
             nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
             K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
             diag_M = diag(ctx.M_sigma); num_fem = size(K_hard, 1);
             hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem))/dt;
             k_hard = full(median(hard_vals(hard_vals>1e-12)));
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
    end
end