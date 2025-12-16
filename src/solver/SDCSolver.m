classdef SDCSolver < handle
    % SDCSOLVER 基于非线性序列扫描的 SDC 求解器 (v12.0 - Nonlinear SISDC + Full Logging)
    % 
    % 核心机制: Semi-Implicit SDC (SISDC)
    %   1. 不做全局线性化预测，而是沿着 GLL 节点逐个求解非线性方程。
    %   2. 每个子步本质上是一个 "带高阶源项修正的 Backward Euler"。
    %   3. 这种方法对变压器饱和等刚性问题极其稳健。
    
    properties
        Assembler       % 有限元组装器
        LinearSolver    % 线性求解器
        
        PolyOrder = 3        % 多项式阶数 (建议 3 或 4)
        MaxSDCIters = 5      % SDC 扫描次数 (SISDC 收敛很快，通常 5 次足够)
        SDCTolerance = 1e-4  % SDC 收敛容差
    end
    
    methods
        function obj = SDCSolver(assembler)
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
            
            fprintf('==================================================\n');
            fprintf('   SDC Solver (Nonlinear SISDC, P=%d)\n', obj.PolyOrder);
            fprintf('==================================================\n');
            
            % --- 1. 初始化与预计算 ---
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            ctx.space_A = space_A; ctx.space_P = space_P;
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            % 初始化探针
            postProc = []; probeB_History = []; 
            if ~isempty(probePoint)
                postProc = PostProcessor(obj.Assembler);
                probeB_History = zeros(length(timeSteps), 1);
            end
            
            % SDC 谱积分矩阵 (GLL)
            fprintf('   [Init] Computing GLL nodes and integration matrix...\n');
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; % 积分矩阵 [NumNodes x NumNodes]
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % 组装不变矩阵 (M, C, W)
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
            
            % 边界条件索引
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
            
            currentHistory = zeros(length(timeSteps), 1);
            globalTime = 0;
            
            % ==================================================
            % 主时间循环 (Time Slab Loop)
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n=== Slab %d/%d (dt=%.1e, T=%.4f -> %.4f) ===\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. Initialization (直接复制初值)
                X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                
                % 2. SDC Sweep (非线性序列扫描)
                for sdc_iter = 1:obj.MaxSDCIters
                    fprintf('   [SDC Sweep %d] Start...\n', sdc_iter);
                    
                    % 保存上一轮的解 (用于计算高阶积分项)
                    X_old = X_slab;
                    
                    % 执行一次完整的序列扫描
                    [X_slab, max_diff_norm, residual_norm] = obj.performSISDCSweep(X_old, x_start, ctx);
                    
                    fprintf('   [SDC Sweep %d] Result: UpdateNorm=%.2e, Residual=%.2e', ...
                        sdc_iter, max_diff_norm, residual_norm);
                    
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' -> Converged.\n');
                        break;
                    else
                        fprintf('\n');
                    end
                end
                
                % Advance
                x_start = X_slab(:, end); 
                globalTime = t_end;
                currentHistory(slab_idx) = x_start(end);
                
                % Probe & Log
                B_mag = 0;
                if ~isempty(postProc)
                    B_val = postProc.probeB(x_start(1:ctx.num_A), probePoint);
                    B_mag = norm(B_val);
                    probeB_History(slab_idx) = B_mag;
                    fprintf('      [End of Slab] |B|=%.4f T, I=%.4f A\n', B_mag, full(x_start(end)));
                end
                
                % Plot
                if ~isempty(monitorHandle)
                    try
                        monitorHandle(globalTime, x_start(end), ...
                            cumsum(timeSteps(1:slab_idx)), currentHistory(1:slab_idx));
                        drawnow limitrate;
                    catch; end
                end
            end
            
            solutionResults = x_start;
            info.FinalTime = globalTime;
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc), info.ProbeB_History = probeB_History; end
        end
    end
    
    methods (Access = private)
        
        function [X_new, max_diff_norm, final_res_norm] = performSISDCSweep(obj, X_k, x0, ctx)
            % PERFORMSISDCSWEEP 执行一次半隐式谱延迟修正扫描
            
            num_nodes = ctx.Spectral.NumNodes;
            num_dofs = ctx.numTotalDofs;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; % 起点固定
            
            % 1. 计算上一轮解 X_k 在所有节点处的物理力 F(X_k)
            fprintf('      - Calculating High-Order Integral Source...\n');
            F_k = obj.evaluateF(X_k, ctx); 
            
            % 2. 计算积分项 I_mat = F_k * Q' * dt/2
            dt_half = dt_slab / 2.0;
            I_accum = (F_k * Q.') * dt_half; % 全局积分
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            % 3. 序列扫描 (Sequential Sweep)
            for m = 2:num_nodes
                % 当前子步长
                tau_curr = nodes(m);
                tau_prev = nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half;
                
                fprintf('      - Solving Node %d/%d (dt_sub=%.2e)...\n', m, num_nodes, dt_sub);
                
                % 准备隐式方程的右端项
                u_prev_new = X_new(:, m-1);
                f_curr_old = F_k(:, m);
                
                % 积分修正增量 (Correction Integral)
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                % RHS = M * u_{m-1} - dt * F_old(u_m) + Integral
                % 这里为了传递给 solveImplicitStep，我们构建 Constant_Term
                % Constant_Term = M * u_prev_new + int_delta - dt_sub * f_curr_old
                % (注意: evaluateF 返回的是 Force，不是 M^{-1}*Force)
                
                % 计算 M * u_prev_new
                Scale = ctx.CircuitRowScale;
                M_sys = ctx.M_sigma; 
                [w_idx, ~, w_val] = find(ctx.W_vec);
                
                % M_tot 乘法 (手动处理以节省内存)
                M_u_prev = M_sys * u_prev_new;
                % 加上电路部分 M (L*I 和 W'*A)
                I_prev = u_prev_new(end);
                % 电路行: (W' * A) * Scale + (L * I) * Scale
                % W_vec 是 FEM 到电路的耦合，也就是 W'
                term_L = ctx.circuit.L * I_prev;
                term_W = ctx.W_vec' * u_prev_new;
                M_u_prev(end) = M_u_prev(end) + (term_L + term_W) * Scale;
                % FEM 行: 电路对场的耦合 (+W*I_prev * Scale)? 
                % 答: M 矩阵中，场方程没有对 I 的微分项(只有 K 有)。
                % 场方程: M*dA/dt + K*A - W*I = 0.
                % 所以场方程的 M 只有 M_sigma。
                % 但是！在 M_cir_row 构建时，我们之前加上了 w_val。
                % 让我们仔细检查 assemble invariant matrices:
                % ctx.M_sigma 只有 mass weighted.
                % 我们需要把 W * I 项加进去吗？不，W*I 是 stiffness 项 (无微分)。
                % 所以场方程的 M 乘法只涉及 M_sigma。
                
                % 构造 Constant RHS Vector
                RHS_Vector = M_u_prev + int_delta - dt_sub * f_curr_old;
                
                % 初始猜测 (使用上一轮的值)
                u_guess = X_k(:, m);
                
                % 执行局部牛顿迭代
                [u_next, res_norm, iter_count] = obj.solveImplicitStep(u_guess, dt_sub, RHS_Vector, nodes(m), ctx);
                
                X_new(:, m) = u_next;
                
                % 记录误差
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
                
                fprintf('        -> Node %d Solved: %d Iters, Res=%.2e, Diff=%.2e\n', ...
                    m, iter_count, res_norm, diff);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            % 计算物理力 F(u) = -K*u - W*I + Source
            % 对应方程: M * du/dt = F(u)
            
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            
            C_tot = ctx.C_AP + ctx.C_AP';
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            
            for m = 1:num_nodes
                x = X_mat(:, m);
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half;
                
                % 计算 K_fem * A
                A_vec = x(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                % 组装力向量
                % F_field = - (K_fem*A + C*A - W*I)
                % 注意符号: 标准方程 M*u' + K*u = Source => M*u' = Source - K*u
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(R_fem)) = -R_fem; 
                f_vec = f_vec - C_tot * x;
                f_vec = f_vec + ctx.W_vec * x(end); % +W*I
                
                % 电路部分
                % Eq: V - R*I - d/dt(Psi + L*I) = 0
                % => d/dt(...) = V - R*I
                % F_cir = (V - R*I)
                V_src = ctx.circuit.V_source_func(t);
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                
                F_val(:, m) = f_vec;
            end
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_node, ctx)
            % 求解局部非线性方程: M*u - dt*F(u) - RHS_Vector = 0
            % [v12.3 Fix] 纯相对误差控制 (Pure Relative Error Control)
            
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            % M_tot (质量矩阵部分)
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_tot = M_sys + M_cir_diag + M_row_coup;
            
            final_res = 1e10;
            res_norm_init = 1.0; % 初始残差
            res_norm_prev = 1e10;
            iter_count = 0;
            
            % --- 收敛参数 ---
            % 只要残差相比初始值下降了 6 个数量级，即认为物理收敛
            REL_TOL = 1e-6;      
            
            % 停滞容差: 如果某一步改善幅度小于 0.1%，说明到了数值底噪，强制停止
            STAG_TOL = 1e-3;     
            
            MAX_ITERS = 15;
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                
                % 1. 组装 Jacobian 和 Residual
                A_vec = u_curr(1:ctx.num_A);
                [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                
                % Residual = M*u + dt*(-Force) - RHS
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
                
                % --- 收敛判断逻辑 ---
                
                % 记录初始残差
                if iter == 1
                    res_norm_init = res_norm;
                    % 防止完美初值导致除零 (极罕见)
                    if res_norm_init < 1e-20, res_norm_init = 1.0; end
                end
                
                % 计算相对误差 (相对于本步迭代的起点)
                rel_err = res_norm / res_norm_init;
                
                % 1. 相对误差达标 -> 收敛
                if rel_err < REL_TOL
                    final_res = res_norm;
                    fprintf('          Newton Converged: RelErr=%.2e\n', rel_err);
                    break;
                end
                
                % 2. 绝对残差极小 -> 收敛 (兜底保护)
                if res_norm < 1e-10
                    final_res = res_norm;
                    break;
                end
                
                % 3. 停滞检测 -> 强制收敛
                % 如果当前残差和上一步几乎没变 (下降停滞)，说明已经到了机器精度或模型精度的极限
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / res_norm_prev;
                    if abs(improvement) < STAG_TOL
                        % 认为已尽力，停止迭代
                        final_res = res_norm;
                        break;
                    end
                end
                
                res_norm_prev = res_norm;
                
                % --- 求解更新 ---
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
            
            % 只有当相对误差依然很大 (例如 > 0.1) 且发散时才警告
            if final_res / res_norm_init > 1e-1 && final_res > 1.0
                fprintf('          [Warn] Local Newton stalled at high residual: Rel=%.2e\n', final_res / res_norm_init);
            end
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
             % (保持不变)
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            
            k_soft_rep = k_hard_median / 1000; % 简化估算
            
            k_target = sqrt(k_hard_median * k_soft_rep);
            Z_cir = abs(ctx.circuit.R + ctx.circuit.L/dt); 
            max_W = max(abs(nonzeros(ctx.W_vec)));         
            max_coupling_val = max_W / dt;               
            max_circuit_val = max(Z_cir, max_coupling_val);
            if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
            scaleFactor = k_target / max_circuit_val;
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
             % (保持不变)
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
    end
end