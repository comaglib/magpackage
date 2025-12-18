classdef SDIRK2ISDCSolver < handle
    % SDIRK2ISDCSolver 基于 SDIRK-2 扫除器的高阶 SDC 求解器
    % 
    % 核心机制: SDIRK-2 Sweeper SDC
    %   1. 基础积分器: 2级2阶 SDIRK (L-stable), gamma = 1 - 1/sqrt(2).
    %   2. 精度提升: 理论上每次扫描提升 2 阶精度 (Sweep 1 -> Order 2, Sweep 2 -> Order 4).
    %   3. 适用场景: 极度刚性且要求高精度的变压器瞬态 (如涌流、短路).
    
    properties
        Assembler       % 有限元组装器
        LinearSolver    % 线性求解器
        
        PolyOrder = 3        % 谱元阶数 P (建议 3)
        MaxSDCIters = 2      % SDC 扫描次数 (由于是二阶Sweeper，通常 2 次即可达到 P=3 的精度)
        SDCTolerance = 1e-5  % SDC 收敛容差
        
        % SDIRK2 参数
        Gamma = 1 - 1/sqrt(2); 
    end
    
    methods
        function obj = SDIRK2ISDCSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 400; % 增加内存预留
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
            fprintf('   SDIRK-2 SDC Solver (Order+2/Sweep, L-Stable)\n');
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
            
            % --- 2. 谱元基础数据 (GLL Nodes) ---
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = SpectralTimeUtils.integration_matrix(gll_nodes); 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % 预计算 SDIRK 扫除所需的插值矩阵 (从 GLL 节点插值到 t_gamma)
            % 我们需要在每个子区间 [t_{m-1}, t_m] 内的 t_gamma 处评估 u^k
            % t_gamma_local = -1 + (gamma * 2) ? 不，这是针对标准区间 [-1,1] 的映射
            % 在 performSweep 中动态插值更灵活。
            
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
            
            % 预组装质量矩阵
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
            
            % 数据记录
            full_time_list = []; full_current_list = []; full_B_list = [];
            probeB_History = zeros(length(timeSteps), 1);
            currentHistory = zeros(length(timeSteps), 1);
            
            globalTime = 0;
            
            % ==================================================
            % 主时间循环 (Slab Loop)
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n----------------------------------------------------------\n');
                fprintf(' Slab %d/%d | dt=%.2e | T=%.5f -> %.5f\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                fprintf('----------------------------------------------------------\n');
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. Initialization (Flat Prediction)
                X_slab = repmat(x_curr, 1, ctx.Spectral.NumNodes);
                
                % 2. SDC Sweep Loop (使用 SDIRK-2 扫除器)
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    fprintf('   [Sweep %d] SDIRK-2 Sweeping...\n', sdc_iter);
                    
                    [X_slab, max_diff_norm, ~] = obj.performSDIRK2Sweep(X_old, x_curr, ctx);
                    
                    fprintf('      -> Sweep Update Norm: %.2e\n', max_diff_norm);
                    
                    if max_diff_norm < obj.SDCTolerance
                        fprintf('      -> SDC Converged.\n');
                        break;
                    end
                end
                
                % 3. Output Data
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
                currentHistory(slab_idx) = x_curr(end);
                if ~isempty(postProc), probeB_History(slab_idx) = B_nodes(end); end
                
                if ~isempty(postProc)
                    fprintf('   [End Slab] |B|=%.4f T, I=%.4f A\n', B_nodes(end), full_current_list(end));
                end
                
                if ~isempty(monitorHandle)
                    try monitorHandle(globalTime, x_curr(end), full_time_list, full_current_list); catch; end
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
        
        function [X_new, max_diff_norm, final_res_norm] = performSDIRK2Sweep(obj, X_k, x0, ctx)
            % PERFORMSDIRK2SWEEP 核心函数: 使用 SDIRK-2 方法对误差进行扫除
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            gamma = obj.Gamma;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            % 1. 计算上一轮解的物理力 F(u^k) (在所有 GLL 节点)
            F_k = obj.evaluateF_Batch(X_k, ctx); 
            
            % 2. 计算高阶谱积分 I_spec = Q * F_k
            dt_half = dt_slab / 2.0;
            I_spec = (F_k * Q.') * dt_half; % [DoFs x Nodes]
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            % SDIRK 级间 Jacobian 复用缓存
            jaco_ctx.J_factor = [];
            jaco_ctx.prev_dt_gamma = 0;
            
            % 3. 逐区间扫描 (Sequential Sweep m=2:N)
            for m = 2:num_nodes
                % 当前子区间 [t_{m-1}, t_m]
                h_sub = (nodes(m) - nodes(m-1)) * dt_half;
                
                u_prev_new = X_new(:, m-1); % 最新解 u_{m-1}^{k+1}
                
                % --- 准备源项修正 (Source Correction) ---
                % 我们需要计算 integral(F(u^k)) 的 SDIRK 近似值
                % I_SDIRK = h * [ (1-gamma)*F(u^k(t_gamma)) + gamma*F(u^k(t_m)) ]
                
                % (a) 插值 u^k 到 t_gamma
                % t_gamma 在 [-1, 1] 空间中的坐标
                tau_m_prev = nodes(m-1);
                tau_m = nodes(m);
                tau_gamma = tau_m_prev + gamma * (tau_m - tau_m_prev);
                
                u_k_gamma = obj.interpolateStateSingle(nodes, X_k, tau_gamma);
                u_k_end   = X_k(:, m); % u^k(t_m)
                
                % (b) 计算 F(u^k) 在这两个点的值
                t_phys_gamma = ctx.t_start + (tau_gamma + 1) * dt_half;
                t_phys_end   = ctx.t_start + (tau_m + 1) * dt_half;
                
                F_k_gamma = obj.evaluateF_Single(u_k_gamma, t_phys_gamma, ctx);
                F_k_end   = F_k(:, m); % 已经算过了
                
                % (c) 计算 SDIRK 低阶积分近似
                I_low = h_sub * ((1-gamma) * F_k_gamma + gamma * F_k_end);
                
                % (d) 计算高阶谱积分增量 (Spectral Integral Increment)
                I_high = I_spec(:, m) - I_spec(:, m-1);
                
                % (e) 源项 Sigma (平均分布到每一步)
                % 我们将修正项 (I_high - I_low)/h 添加到 SDIRK 的 F 中
                sigma = (I_high - I_low) / h_sub;
                
                % --- 执行 SDIRK-2 步进 ---
                % 求解方程: u' = F(u) + sigma
                % 初值: u_prev_new
                
                % Stage 1: Solve for g1 at t_gamma
                % g1 = u_prev + h*gamma*( F(g1) + sigma )
                % Eq: M*g1 - h*gamma*F(g1) = M*u_prev + h*gamma*M*sigma (注意 sigma 是力，这里简化处理)
                % 实际上 RHS 应该是: M*u_prev + h*gamma*(sigma_force)
                
                dt_eff = h_sub * gamma;
                
                % 构建 RHS for Stage 1
                % M * u_prev + dt_eff * sigma
                % 注意: sigma 已经是力向量了，但它包含了 mass 矩阵的逆效应吗？不。
                % evaluateF 返回的是 M^{-1} * Force? 不，evaluateF 返回的是 Force 向量。
                % 我们的方程是 M u' = F_vec.
                % 所以 M * (u' - sigma) = F_vec. => M u' = F_vec + M*sigma.
                % 这里的 sigma 是从 (I_high - I_low)/h 算出来的，I 也是 Force*dt。
                % 所以 sigma 是 Force 维度的。
                % RHS_1 = M * u_prev + dt_eff * sigma
                
                M_u_prev = ctx.M_full_const * u_prev_new;
                RHS_1 = M_u_prev + dt_eff * sigma;
                
                % 猜测初值 (Rolling Guess)
                u_guess = u_prev_new; 
                
                % Solve Stage 1
                fprintf('      [Node %d Stage 1] ', m);
                [g1, ~, ~, jaco_ctx] = obj.solveImplicitStep(u_guess, dt_eff, RHS_1, t_phys_gamma, ctx, jaco_ctx);
                
                % Force Recovery for g1 (avoid re-evaluation)
                % F(g1) = (M*g1 - RHS_1)/dt_eff
                if abs(dt_eff) > 1e-12
                    F_g1 = (ctx.M_full_const * g1 - RHS_1) / dt_eff;
                else
                    F_g1 = obj.evaluateF_Single(g1, t_phys_gamma, ctx);
                end
                % 注意: 这里恢复出的 F_g1 实际上是 (F(g1) + sigma)。
                % 但我们在 Stage 2 公式中需要的是 (F(g1) + sigma)。所以直接用！
                F_stage1_effective = F_g1; 
                
                % Stage 2: Solve for u_new at t_m (End of step)
                % u_new = u_prev + h*(1-gamma)*(F(g1)+sigma) + h*gamma*(F(u_new)+sigma)
                % M*u_new - h*gamma*F(u_new) = M*u_prev + h*(1-gamma)*M*(F(g1)+sigma) + h*gamma*M*sigma
                % Wait, F_stage1_effective already includes sigma.
                % M*u_new - h*gamma*(F(u_new)+sigma) = M*u_prev + h*(1-gamma)*M*F_stage1_effective
                % RHS_2 = M*u_prev + h*(1-gamma)*F_stage1_effective*M? No, F is force vector.
                % RHS_2 = M_u_prev + h*(1-gamma)*F_stage1_effective + h*gamma*sigma.
                
                RHS_2 = M_u_prev + h_sub * (1-gamma) * F_stage1_effective + dt_eff * sigma;
                
                % Solve Stage 2
                fprintf('      [Node %d Stage 2] ', m);
                % Jacobina Reuse: dt_eff is same as Stage 1 (h*gamma)
                [u_next, res_norm, iter_count, jaco_ctx] = obj.solveImplicitStep(g1, dt_eff, RHS_2, t_phys_end, ctx, jaco_ctx);
                
                % Update
                X_new(:, m) = u_next;
                
                % Error check
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
                
                fprintf('        Node %d Done. Iters=%d\n', m, iter_count);
            end
        end
        
        function F_val = evaluateF_Batch(obj, X_mat, ctx)
            % Batch Evaluate F (Standard for SDC start)
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
            % Single State Evaluate F
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
        
        function [u_new, final_res, iter_count, jaco_ctx] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_phys, ctx, jaco_ctx)
            % 局部牛顿求解器 (SDIRK Optimized with Jacobian Reuse)
            
            u_curr = u_guess;
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            
            % Constant stiffness
            [w_idx, ~, w_val] = find(ctx.W_vec);
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            K_const_part = ctx.C_AP + ctx.C_AP' + K_col_coup + K_cir_diag;
            
            final_res = 1e10;
            res_norm_init = 1.0; res_norm_prev = 1e10;
            iter_count = 0;
            
            % Jacobian Update Logic
            update_jacobian = true;
            % 如果上一级传来的 dt 相同，且有缓存，尝试复用
            if isfield(jaco_ctx, 'prev_dt_gamma') && abs(dt - jaco_ctx.prev_dt_gamma) < 1e-9 && ~isempty(jaco_ctx.J_factor)
                 update_jacobian = false; % Reuse J from Stage 1
                 % fprintf('J');
            else
                 % fprintf('U');
            end
            
            MAX_ITERS = 15;
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                A_vec = u_curr(1:ctx.num_A);
                
                % 1. Assemble Residual (Fast)
                if update_jacobian
                    [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                else
                    [~, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                    K_fem_small = [];
                end
                
                % 2. Residual Vector
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                neg_Force = neg_Force - ctx.W_vec * u_curr(end); 
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * u_curr(end) - V_src) * Scale;
                
                Res = ctx.M_full_const * u_curr + dt * neg_Force - RHS_Vector;
                res_norm = norm(Res);
                
                % 3. Check Convergence
                if iter == 1, res_norm_init = res_norm; if res_norm_init < 1e-20, res_norm_init = 1; end; end
                
                if res_norm / res_norm_init < 1e-6 || res_norm < 1e-9
                    final_res = res_norm; break; 
                end
                
                if iter > 1 && abs(res_norm_prev - res_norm)/res_norm_prev < 1e-3
                    final_res = res_norm; break; % Stagnation accepted
                end
                res_norm_prev = res_norm;
                
                % 4. Solve Linear System
                if update_jacobian
                    [ki, kj, kv] = find(K_fem_small);
                    K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                    Jac = ctx.M_full_const + dt * (K_fem + K_const_part);
                    
                    jaco_ctx.J_factor = Jac;
                    jaco_ctx.prev_dt_gamma = dt;
                    update_jacobian = false; 
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(jaco_ctx.J_factor, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                u_curr = u_curr + du;
            end
            
            u_new = u_curr;
        end
        
        function u_interp = interpolateStateSingle(~, nodes, X_mat, tau_target)
            % Lagrange Interpolation for a single time point
            N = length(nodes);
            u_interp = zeros(size(X_mat, 1), 1);
            
            for j = 1:N
                lj = 1;
                for i = 1:N
                    if i ~= j
                        lj = lj * (tau_target - nodes(i)) / (nodes(j) - nodes(i));
                    end
                end
                u_interp = u_interp + X_mat(:, j) * lj;
            end
        end
        
        % ... Helper functions (calculateCircuitScale, getBoundNuVec, interpolatePolynomial) ...
        % (Use same implementation as SISDCSolver)
        
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