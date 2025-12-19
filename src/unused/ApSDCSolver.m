classdef ApSDCSolver < handle
    % ApSDCSolver 变压器低频电磁场 SDC 求解器 (v15.2 - Hybrid Safe)
    %
    % 核心机制: Hybrid Newton Strategy (混合牛顿策略)
    %
    % [v15.2 修复与优化]:
    %   1. [Hybrid Strategy]: 
    %      - Sweep 1: 强制 Full Newton (步步更新 J)。确保在强非线性区不发散。
    %      - Sweep 2+: 强制 Modified Newton (冻结 J)。利用已建立的刚度场快速收敛。
    %   2. [Safety]: 修复了 v15.1 的崩溃问题。当线搜索失败(残差发散)时，
    %      拒绝本次步长更新，回滚状态，并强制下一式更新 Jacobian。
    %   3. [Speed]: 通过 Sweep 2+ 的冻结策略，在后续扫描中获得显著提速。
    
    properties
        Assembler       
        LinearSolver    
        
        PolyOrder = 3        % P=3
        MaxSDCIters = 5      % SDC 扫描次数
        SDCTolerance = 1e-4  % SDC 整体收敛容差
        
        RelaxFactor = 1.0    % 基础松弛因子
        StagnationTol = 0.01 % 停滞阈值
    end
    
    methods
        function obj = ApSDCSolver(assembler)
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
            
            fprintf('==========================================================\n');
            fprintf('   SISDC Solver (v15.2 - Hybrid Safe Strategy)\n');
            fprintf('   Strategy: Sweep1(Full) + Sweep2+(Frozen) + Diverge Protect\n');
            fprintf('==========================================================\n');
            
            % --- 1. 上下文构建 ---
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~isempty(space_P) && ~dofHandler.DofMaps.isKey(space_P.toString())
                dofHandler.distributeDofs(space_P); 
            end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            ctx.space_A = space_A; ctx.space_P = space_P; 
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            postProc = []; 
            if ~isempty(probePoint), postProc = PostProcessor(obj.Assembler); end
            full_time_list = []; full_current_list = []; full_B_list = [];
            
            % --- 2. SDC 积分矩阵 ---
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            tau_plot = linspace(-1, 1, 20)';
            
            % --- 3. 组装恒定矩阵 ---
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_field = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            if ~isempty(space_P)
                C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
                [ci, cj, cv] = find(C_AP_local);
                ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            else
                ctx.C_AP = sparse(ctx.numTotalDofs, ctx.numTotalDofs);
            end
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            Scale = ctx.CircuitRowScale;
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', Scale);
            
            M_cir_diag = sparse(ctx.numTotalDofs, ctx.numTotalDofs, ctx.circuit.L * Scale, ctx.numTotalDofs, ctx.numTotalDofs);
            M_coup_row = sparse(repmat(ctx.numTotalDofs, length(w_idx), 1), w_idx, w_val * Scale, ctx.numTotalDofs, ctx.numTotalDofs);
            ctx.M_sys_total = ctx.M_field + M_cir_diag + M_coup_row;
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
            x_prev_slab_start = x_start;
            dt_prev = 0;
            globalTime = 0;
            
            % ==================================================
            % Time Slab Loop
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                fprintf('\n=== Slab %d/%d (dt=%.2e, T=%.4f -> %.4f) ===\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                % --- 1. Predictor (Safe Linear Extrapolation) ---
                X_slab = zeros(ctx.numTotalDofs, ctx.Spectral.NumNodes);
                if slab_idx > 1 && dt_prev > 0
                    slope = (x_start - x_prev_slab_start) / dt_prev * 0.95; 
                    node_times = (ctx.Spectral.Nodes + 1) / 2 * dt_slab; 
                    for k = 1:ctx.Spectral.NumNodes
                        X_slab(:, k) = x_start + slope * node_times(k);
                    end
                else
                    X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                end
                
                % 2. SISDC Sweep
                diff_history = [];
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    
                    [X_slab, max_diff_norm] = obj.performSISDCSweep(X_old, x_start, ctx, sdc_iter);
                    
                    if sdc_iter == 1
                        fprintf('  Swp 1 Diff: %.2e', max_diff_norm);
                    else
                        fprintf(' -> %.2e', max_diff_norm);
                    end
                    
                    % SDC Level Stagnation Check
                    diff_history = [diff_history, max_diff_norm]; %#ok<AGROW>
                    if sdc_iter >= 3
                        improvement = diff_history(end-1) / (diff_history(end) + 1e-20);
                        if improvement < 1.05 % <5% improvement
                            fprintf(' [Stagnated]'); break;
                        end
                    end
                    
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' [Converged]'); break;
                    end
                    if isnan(max_diff_norm) || max_diff_norm > 1e10
                        warning('SDC Diverged!'); break;
                    end
                end
                fprintf('\n');
                
                x_prev_slab_start = x_start;
                dt_prev = dt_slab;
                
                t_plot_phys = t_start + (tau_plot + 1) / 2 * dt_slab;
                I_nodes = X_slab(end, :)';
                I_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, I_nodes, tau_plot);
                
                B_nodes = zeros(ctx.Spectral.NumNodes, 1);
                if ~isempty(postProc)
                    for k = 1:ctx.Spectral.NumNodes
                        B_val = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_nodes(k) = norm(B_val);
                    end
                    fprintf('      [End] |B|=%.4f T, I=%.4f A\n', B_nodes(end), I_nodes(end));
                end
                B_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, B_nodes, tau_plot);
                
                idx_s = 2; if slab_idx == 1, idx_s = 1; end
                full_time_list = [full_time_list; t_plot_phys(idx_s:end)];
                full_current_list = [full_current_list; I_plot(idx_s:end)];
                full_B_list = [full_B_list; B_plot(idx_s:end)];
                
                x_start = X_slab(:, end);
                globalTime = t_end;
                if ~isempty(monitorHandle)
                    try monitorHandle(globalTime, x_start(end), full_time_list, full_current_list); drawnow limitrate; catch; end;
                end
            end
            
            solutionResults = x_start;
            info.Time_Full = full_time_list;
            info.Current_Full = full_current_list;
            info.ProbeB_Full = full_B_list;
        end
    end
    
    methods (Access = private)
        
        function [X_new, max_diff_norm] = performSISDCSweep(obj, X_k, x0, ctx, sdc_iter)
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            X_new = X_k; X_new(:, 1) = x0; 
            
            F_k = obj.evaluateF(X_k, ctx); 
            dt_half = dt_slab / 2.0;
            I_accum = (F_k * Q.') * dt_half;
            max_diff_norm = 0;
            
            for m = 2:num_nodes
                tau_curr = ctx.Spectral.Nodes(m);
                tau_prev = ctx.Spectral.Nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half;
                
                u_prev_solved = X_new(:, m-1); 
                f_curr_old_k  = F_k(:, m);     
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                M_u_prev = ctx.M_sys_total * u_prev_solved;
                RHS_Const = M_u_prev + int_delta - dt_sub * f_curr_old_k;
                
                % [Hybrid Newton Strategy]
                if sdc_iter == 1
                    % Sweep 1: 必须稳健，全量更新 J，迭代次数给足
                    u_guess = u_prev_solved; 
                    tol_newton = 1e-1; % 初始粗容差
                    max_newton = 15;   
                else
                    % Sweep 2+: 此时解已在附近，可以冻结 J 提速
                    u_guess = X_k(:, m);     
                    tol_newton = 1e-5;
                    max_newton = 1; % 单步修正
                end
                
                u_next = obj.solveNewtonStep(u_guess, dt_sub, RHS_Const, tau_curr, ctx, tol_newton, max_newton, m, sdc_iter);
                
                X_new(:, m) = u_next;
                
                nd = norm(u_next - X_k(:, m));
                max_diff_norm = max(max_diff_norm, nd / (norm(u_next) + 1e-6));
            end
        end
        
        function u_new = solveNewtonStep(obj, u_guess, dt, RHS_Const, t_node, ctx, tol, max_iters, node_idx, sdc_iter)
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            t_phys = ctx.t_start + (t_node + 1) * ctx.dt_slab / 2.0; 
            
            res_norm_prev = 1e20;
            J_cache = []; % 雅可比缓存
            
            for iter = 1:max_iters
                % 1. Calc Residual
                [Res, neg_Force, ~] = obj.calcResidual(u_curr, dt, RHS_Const, t_phys, ctx);
                res_norm = norm(Res);
                
                if sdc_iter == 1 && max_iters > 1
                    % fprintf('      [Node %d] Newton %d: |Res|=%.4e\n', node_idx, iter, res_norm);
                end
                
                % Check Convergence
                if iter == 1 && res_norm < tol
                    % if sdc_iter == 1 && max_iters > 1, fprintf(' (Converged)\n'); end
                    u_new = u_curr; return; 
                end
                
                % Stagnation Check
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / (res_norm_prev + 1e-20);
                    if improvement < obj.StagnationTol && improvement > -1e-5
                        if sdc_iter == 1, fprintf(' (Stagnated: %.1f%%) -> Accept.\n', improvement*100); end
                        u_new = u_curr; return; 
                    elseif improvement < -0.1
                        % [Warning] Diverging detected
                        if sdc_iter == 1, fprintf(' (Diverging! Force Update J)\n'); end
                        J_cache = []; % 强行清除缓存，强制下一次更新
                    end
                end
                
                res_norm_prev = res_norm;
                
                % 2. Jacobian Update Strategy
                update_J = false;
                
                % Sweep 1: 始终更新 J (Full Newton)
                if sdc_iter == 1
                    update_J = true;
                else
                    % Sweep 2+: 默认冻结，仅第一步更新
                    if iter == 1, update_J = true; end
                end
                
                % 如果缓存为空，必须更新
                if isempty(J_cache), update_J = true; end
                
                if update_J
                    A_vec = u_curr(1:ctx.num_A);
                    [K_fem_small, ~] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                    
                    num_dofs = ctx.numTotalDofs;
                    [ki, kj, kv] = find(K_fem_small);
                    K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                    [w_idx, ~, w_val] = find(ctx.W_vec);
                    K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
                    K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
                    
                    K_tan = K_fem + (ctx.C_AP + ctx.C_AP') + K_col_coup + K_cir_diag;
                    Jac = ctx.M_sys_total + dt * K_tan;
                    
                    [J_bc, ~] = BoundaryCondition.applyDirichlet(Jac, Res, ctx.fixedDofs);
                    J_cache = J_bc;
                end
                
                % 3. Solve
                [~, Res_bc] = BoundaryCondition.applyDirichlet(J_cache, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_cache, -Res_bc);
                
                % 4. Backtracking Line Search (Safety Net)
                alpha = obj.RelaxFactor;
                u_try = u_curr;
                line_search_ok = false;
                
                for k = 1:5 
                    u_try = u_curr + alpha * du;
                    [Res_try, ~, ~] = obj.calcResidual(u_try, dt, RHS_Const, t_phys, ctx);
                    
                    if norm(Res_try) < res_norm || norm(Res_try) < tol
                        line_search_ok = true; break;
                    else
                        alpha = alpha * 0.5;
                    end
                end
                
                % [Critical Fix] 如果线搜索失败，拒绝本次更新！
                if ~line_search_ok
                    % if sdc_iter == 1, fprintf(' [LS Fail - Reject Step] \n'); end
                    % 保持 u_curr 不变，清除 J_cache，强制下一次重组
                    J_cache = [];
                    continue; 
                end
                
                u_curr = u_try;
            end
            if sdc_iter == 1 && max_iters > 1, fprintf('\n'); end
            u_new = u_curr;
        end
        
        function [Res, neg_Force, R_fem_force] = calcResidual(obj, u, dt, RHS_Const, t_phys, ctx)
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            A_vec = u(1:ctx.num_A);
            [~, R_fem_force] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            neg_Force = sparse(num_dofs, 1);
            neg_Force(1:length(R_fem_force)) = R_fem_force; 
            I_val = u(end);
            neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u;
            neg_Force = neg_Force - ctx.W_vec * I_val;
            V_src = ctx.circuit.V_source_func(t_phys);
            neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
            Res = ctx.M_sys_total * u + dt * neg_Force - RHS_Const;
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            dt_half = ctx.dt_slab / 2.0;
            for m = 1:num_nodes
                x = X_mat(:, m);
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half;
                [~, neg_Force, ~] = obj.calcResidual(x, 0, 0, t, ctx);
                F_val(:, m) = -neg_Force; 
            end
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            meshTags = obj.Assembler.Mesh.RegionTags;
            nu_vec = zeros(length(meshTags), 1);
            matLib = ctx.matLib;
            nu0 = 1/(4*pi*1e-7);
            if isa(matLib, 'containers.Map')
                keys = matLib.keys;
                for i=1:length(keys)
                    k=keys{i}; if ischar(k), kv=str2double(k); else, kv=k; end
                    mat=matLib(k);
                    val = nu0; if ~strcmp(mat.Type,'Linear'), val=nu0/1000; end
                    nu_vec(meshTags==kv) = val;
                end
            end
            
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec);
            diag_M = diag(ctx.M_field); 
            num_fem = size(K_hard, 1); if num_fem>length(diag_M), num_fem=length(diag_M); end
            
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem))/dt;
            vals_nz = hard_vals(hard_vals>1e-12);
            if isempty(vals_nz), med=1.0; else, med=full(median(vals_nz)); end
            
            target = sqrt(med * (med/1000));
            max_cir = max([abs(ctx.circuit.R+ctx.circuit.L/dt), max(abs(nonzeros(ctx.W_vec)))/dt]);
            if max_cir<1e-20, max_cir=1.0; end
            scaleFactor = target / max_cir;
        end
        
        function vals_interp = interpolatePolynomial(~, nodes, values, target_nodes)
            N = length(nodes); M = length(target_nodes);
            vals_interp = zeros(M, 1);
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