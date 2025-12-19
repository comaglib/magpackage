classdef LISDCSolver < handle
    % LISDCSolver 变压器低频电磁场 LI-SDC 求解器 (v13.3 - Damped & Stable)
    %
    % 核心算法: Linearly Implicit Spectral Deferred Correction (LI-SDC)
    % 
    % 修复日志 (v13.3):
    %   1. [Stability] 引入 RelaxFactor (松弛因子)。对于变压器饱和等刚性问题，
    %      全量修正(Newton Step)会导致震荡，必须引入阻尼 (Damping)。
    %   2. [Robustness] 增加了对线性求解器结果的 NaN/Inf 检查。
    %   3. [Logic] 修正了 RHS 计算中潜在的维度不匹配风险。
    
    properties
        Assembler       % Finite Element Assembler
        LinearSolver    % Linear Solver Wrapper (MUMPS/Pardiso)
        
        PolyOrder = 3        % 多项式阶数 P
        MaxSDCIters = 15     % SDC 最大扫描次数
        SDCTolerance = 1e-4  % SDC 收敛容差
        RelaxFactor = 0.1    % [重要] 松弛因子 (0.5 ~ 0.8)，防止刚度剧烈变化时发散
    end
    
    methods
        function obj = LISDCSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; % 非对称
            obj.LinearSolver.MumpsICNTL.i14 = 400; % 增加内存
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
            fprintf('   LI-SDC Solver (v13.3 - Damped LI-SDC)\n');
            fprintf('   PolyOrder=%d, Relax=%.2f\n', obj.PolyOrder, obj.RelaxFactor);
            fprintf('==========================================================\n');
            
            % --- 1. 上下文与预计算 ---
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~isempty(space_P) && ~dofHandler.DofMaps.isKey(space_P.toString())
                dofHandler.distributeDofs(space_P); 
            end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            ctx.space_A = space_A; 
            ctx.space_P = space_P; 
            ctx.matLib = matLib; 
            ctx.circuit = circuitProps;
            
            % 后处理
            postProc = []; 
            if ~isempty(probePoint), postProc = PostProcessor(obj.Assembler); end
            
            % 数据记录
            full_time_list = []; full_current_list = []; full_B_list = [];
            
            % --- 2. SDC 积分矩阵 ---
            fprintf('   [Init] Computing GLL nodes and integration matrix...\n');
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            tau_plot = linspace(-1, 1, 20)';
            
            % --- 3. 组装恒定矩阵 ---
            fprintf('   [Init] Assembling invariant system matrices...\n');
            
            % M_field (通常为0，除非有导电体)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_field = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % C_AP (Lagrange Multiplier)
            if ~isempty(space_P)
                C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
                [ci, cj, cv] = find(C_AP_local);
                ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            else
                ctx.C_AP = sparse(ctx.numTotalDofs, ctx.numTotalDofs);
            end
            
            % W_vec (Winding Coupling)
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            
            % Circuit Scale
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            Scale = ctx.CircuitRowScale;
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', Scale);
            
            % M_sys_total (Mass Matrix)
            % Row 2: L * dI/dt + W^T * dA/dt = ...
            M_cir_diag = sparse(ctx.numTotalDofs, ctx.numTotalDofs, ctx.circuit.L * Scale, ctx.numTotalDofs, ctx.numTotalDofs);
            M_coup_row = sparse(repmat(ctx.numTotalDofs, length(w_idx), 1), w_idx, w_val * Scale, ctx.numTotalDofs, ctx.numTotalDofs);
            ctx.M_sys_total = ctx.M_field + M_cir_diag + M_coup_row;
            
            % 边界条件
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            % 初始状态
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
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
                
                % 1. Predictor: Copy previous state (Constant)
                X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                
                % 2. SDC Sweep
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    
                    % *** 执行带阻尼的 LI-SDC ***
                    [X_slab, max_diff_norm] = obj.performLinearlyImplicitSweep(X_old, x_start, ctx);
                    
                    if sdc_iter == 1
                         fprintf('\n   Slab %d (dt=%.1e): Swp 1 (Diff=%.2e)\n', slab_idx, dt_slab, max_diff_norm);
                    else
                         fprintf(' -> %d (%.2e)\n', sdc_iter, max_diff_norm);
                    end
                    
                    % 稳定性熔断
                    if max_diff_norm > 1e10 || isnan(max_diff_norm)
                        warning('SDC Diverged. Breaking loop.');
                        break;
                    end
                    
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' [Converged]');
                        break;
                    end
                end
                fprintf('\n');
                
                % --- Post-Processing ---
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
        
        function [X_new, max_diff_norm] = performLinearlyImplicitSweep(obj, X_old, x0, ctx)
            % 线性化 SDC 扫描 (带阻尼)
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            
            X_new = X_old; 
            X_new(:, 1) = x0; 
            
            % 1. 全量评估 F(u_old)
            F_old = obj.evaluateF(X_old, ctx);
            
            % 2. 积分修正项 (Integral of F using old solution)
            dt_half = dt_slab / 2.0;
            I_accum = (F_old * Q.') * dt_half; 
            
            max_diff_norm = 0;
            
            % 3. 逐子步求解
            for m = 2:num_nodes
                tau_curr = ctx.Spectral.Nodes(m);
                tau_prev = ctx.Spectral.Nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half;
                
                u_old = X_old(:, m);          % Linearization point
                u_prev_new = X_new(:, m-1);   % Already updated previous node
                
                % A. LHS Jacobian: M - dt * J
                [K_total, ~] = obj.assembleTangentStiffness(u_old, ctx);
                LHS_Mat = ctx.M_sys_total + dt_sub * K_total;
                
                % B. RHS Defect: (M*u_prev + Integral_Delta) - M*u_old
                M_u_prev = ctx.M_sys_total * u_prev_new;
                M_u_old  = ctx.M_sys_total * u_old;
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                RHS_Vector = (M_u_prev + int_delta) - M_u_old;
                
                % C. Solve Linear System
                [LHS_bc, RHS_bc] = BoundaryCondition.applyDirichlet(LHS_Mat, RHS_Vector, ctx.fixedDofs);
                
                delta_u = obj.LinearSolver.solve(LHS_bc, RHS_bc);
                
                % Check for Failure
                if isempty(delta_u) || any(isnan(delta_u)) || norm(delta_u) > 1e15
                    warning('Linear solve failed or exploded at node %d', m);
                    delta_u = zeros(size(u_old));
                end
                
                % D. Apply Update with Relaxation (Damping)
                % 防止饱和区刚度突变导致的震荡
                u_new = u_old + obj.RelaxFactor * delta_u;
                X_new(:, m) = u_new;
                
                nd = norm(delta_u);
                rel_upd = nd / (norm(u_new) + 1e-6);
                max_diff_norm = max(max_diff_norm, rel_upd);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            % F(u) = -K*u - C*u + ...
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            
            C_sym = ctx.C_AP + ctx.C_AP'; 
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            
            for m = 1:num_nodes
                x = X_mat(:, m);
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half;
                
                % Field Force
                A_vec = x(1:ctx.num_A);
                [~, K_times_A] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(K_times_A)) = -K_times_A; 
                f_vec = f_vec - C_sym * x;           
                f_vec = f_vec + ctx.W_vec * x(end);  
                
                % Circuit Force
                V_src = ctx.circuit.V_source_func(t);
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                
                F_val(:, m) = f_vec;
            end
        end
        
        function [K_total, Res_check] = assembleTangentStiffness(obj, u, ctx)
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            
            A_vec = u(1:ctx.num_A);
            [K_fem_small, ~] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            [ki, kj, kv] = find(K_fem_small);
            K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
            
            K_base = K_fem + (ctx.C_AP + ctx.C_AP');
            
            [w_idx, ~, w_val] = find(ctx.W_vec);
            % F_field = ... + W*I.  dF/dI = W.  K_tan = -J = -W.
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            
            K_total = K_base + K_col_coup + K_cir_diag;
            Res_check = [];
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            
            diag_M = diag(ctx.M_field); 
            num_fem = size(K_hard, 1);
            if num_fem > length(diag_M), num_fem = length(diag_M); end
            
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            vals_nz = hard_vals(hard_vals > 1e-12);
            if isempty(vals_nz), k_hard_median = 1.0; else, k_hard_median = full(median(vals_nz)); end
            
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
            
            if isa(matLib, 'containers.Map')
                keys = matLib.keys;
                for i = 1:length(keys)
                    tag = keys{i};
                    if ischar(tag), tagVal = str2double(tag); else, tagVal = tag; end
                    mat = matLib(tag);
                    
                    if strcmp(mat.Type, 'Linear')
                        val = mat.Nu_Linear;
                    else
                        if strcmp(mode, 'hard'), val = nu0; else, val = nu0 / 1000.0; end
                    end
                    nu_vec(meshTags == tagVal) = val;
                end
            end
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