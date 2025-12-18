classdef ISDCAlignSolver < handle
    % ISDCAlignSolver 基于波形拐点对齐的自适应 SDC 求解器 (v1.0)
    % 
    % 核心特性:
    %   1. Grid Alignment: 自动检测波形拐点(如饱和起始点)，并将时间步边界对齐到拐点。
    %      这消除了高阶多项式拟合非光滑拐点时的 Gibbs 振荡，解决了收敛困难。
    %   2. Low-Cost Probe: 利用单次扫描(Sweep 1)的结果快速判断波形特征。
    %   3. Curvature Detection: 利用二阶导数极值定位拐点。
    %   4. Robust Predictor: 采用 C1 预测结合有限差分导数，保证刚性问题的稳定性。
    
    properties
        Assembler       % Finite Element Assembler
        LinearSolver    % Linear Solver Wrapper
        
        PolyOrder = 3        % 多项式阶数 (建议 3 或 4)
        MaxSDCIters = 10     % 最大扫描次数
        SDCTolerance = 1e-4  % 收敛容差
        
        % --- 自适应控制参数 ---
        CurvatureThreshold = 100.0  % 曲率阈值 (判定拐点灵敏度，需根据电流幅值调整)
        MinStepSize = 1e-6          % 最小允许步长 (防止过度切分)
    end
    
    methods
        function obj = ISDCAlignSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 300; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 dt_preset, t_max, ... 
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            % SOLVE 主求解函数
            % 输入变化: 接收 dt_preset (预设步长) 和 t_max (总时间)，而非固定数组
            
            fprintf('==================================================\n');
            fprintf('   SDC Solver (Waveform Alignment, P=%d)\n', obj.PolyOrder);
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
            
            % 预计算谱矩阵
            fprintf('   [Init] Computing GLL nodes & matrices...\n');
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            D_mat = obj.computeDifferentiationMatrix(gll_nodes); % 用于导数计算
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.D = D_mat;
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % 组装不变矩阵
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
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            % 初始状态
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 12 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
            % 数据收集容器
            full_time_list = []; 
            full_current_list = []; 
            full_B_list = [];
            
            % 运行状态变量
            current_time = 0;
            slab_count = 0;
            u_dot_prev = zeros(ctx.numTotalDofs, 1); % 上一步的导数(用于预测)
            
            % ==================================================
            % 主时间循环 (Adaptive Loop)
            % ==================================================
            while current_time < t_max - 1e-9
                slab_count = slab_count + 1;
                
                % 1. 确定试探性步长 (Trial Step)
                % 默认尝试走完预设步长，或者直到仿真结束
                dt_try = min(dt_preset, t_max - current_time);
                
                % [核心逻辑] 探测阶段 (Probe Phase)
                fprintf('\n--- Slab %d Probe (Start t=%.4f, Try dt=%.1e) ---\n', slab_count, current_time, dt_try);
                
                % 执行快速探测 (仅 1 次 Sweep)
                [~, has_inflection, t_split] = obj.probeForInflection(x_start, u_dot_prev, dt_try, current_time, ctx);
                
                % 2. 决策阶段 (Decision Phase)
                if has_inflection && t_split > current_time + obj.MinStepSize && t_split < current_time + dt_try - obj.MinStepSize
                    % >> 发现拐点，执行分割 <<
                    dt_final = t_split - current_time;
                    fprintf('   [ALIGN] Inflection detected at t=%.5f. Splitting step (dt=%.2e).\n', t_split, dt_final);
                else
                    % >> 无明显拐点，接受试探步长 <<
                    dt_final = dt_try;
                    fprintf('   [ACCEPT] Smooth waveform. Proceeding with full step.\n');
                end
                
                % 3. 精细求解 (Refine Solve)
                % 使用确定的 dt_final 进行计算直到收敛
                ctx.dt_slab = dt_final;
                ctx.t_start = current_time;
                ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_final);
                
                % 初始猜测 (C1 Predictor)
                if current_time == 0
                     X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                else
                     % 线性外推: x = x0 + v*t
                     t_local = (ctx.Spectral.Nodes + 1)/2 * dt_final;
                     X_slab = x_start + u_dot_prev * t_local';
                end
                
                % SDC 迭代
                converged = false;
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    [X_slab, max_diff, ~] = obj.performISDCSweep(X_old, x_start, ctx);
                    
                    if max_diff < obj.SDCTolerance
                        fprintf('   -> Converged at Iter %d (Update=%.2e).\n', sdc_iter, max_diff);
                        converged = true;
                        break;
                    end
                end
                
                if ~converged
                    fprintf('   [WARN] Step did not fully converge (Update=%.2e).\n', max_diff);
                end
                
                % 4. 后处理与数据收集
                % 计算物理时间节点
                t_nodes_phys = current_time + (ctx.Spectral.Nodes + 1)/2 * dt_final;
                
                % 提取物理量
                I_nodes = X_slab(end, :)';
                
                B_nodes = zeros(ctx.Spectral.NumNodes, 1);
                if ~isempty(postProc)
                    for k = 1:ctx.Spectral.NumNodes
                        B_val = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_nodes(k) = norm(B_val);
                    end
                end
                
                % 存入全量列表 (用于 Visual Reconstruction)
                % 避免重复存储连接点
                if isempty(full_time_list)
                    full_time_list = t_nodes_phys;
                    full_current_list = I_nodes;
                    full_B_list = B_nodes;
                else
                    full_time_list = [full_time_list; t_nodes_phys(2:end)];
                    full_current_list = [full_current_list; I_nodes(2:end)];
                    full_B_list = [full_B_list; B_nodes(2:end)];
                end
                
                % 计算末端导数 (Robust Finite Difference)
                % 这种方法比谱导数更稳定，避免过冲
                u_dot_prev = (X_slab(:, end) - X_slab(:, end-1)) / (t_nodes_phys(end) - t_nodes_phys(end-1));
                
                % 更新状态
                x_start = X_slab(:, end);
                current_time = current_time + dt_final;
                
                % 5. 实时绘图 (使用 PCHIP 平滑重建)
                if ~isempty(monitorHandle)
                   try
                       % 生成用于画图的密集网格
                       t_plot = linspace(full_time_list(1), full_time_list(end), length(full_time_list)*3);
                       I_plot = pchip(full_time_list, full_current_list, t_plot);
                       monitorHandle(current_time, x_start(end), t_plot, I_plot); 
                       drawnow limitrate; 
                   catch; end
                end
            end
            
            % 返回结果
            solutionResults = x_start;
            info.FinalTime = current_time;
            
            % 进行最终的 C1 平滑插值用于输出
            t_out_grid = linspace(0, t_max, 1000)'; % 默认输出 1000 个点
            % 只有当数据点足够时才插值
            if length(full_time_list) > 2
                info.Time_Full = t_out_grid;
                info.Current_Full = pchip(full_time_list, full_current_list, t_out_grid);
                info.ProbeB_Full = pchip(full_time_list, full_B_list, t_out_grid);
            else
                info.Time_Full = full_time_list;
                info.Current_Full = full_current_list;
                info.ProbeB_Full = full_B_list;
            end
        end
    end
    
    methods (Access = private)
        
        function [X_probe, found, t_split] = probeForInflection(obj, x_start, u_dot_prev, dt_try, t_start, ctx)
            % PROBEFORINFLECTION 快速探测波形并查找拐点
            
            % 1. 设置上下文
            ctx.dt_slab = dt_try;
            ctx.t_start = t_start;
            % 重新计算 Scale，因为 dt 变了
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_try);
            
            % 2. 构造初始猜测 (Predictor)
            if t_start == 0
                X_guess = repmat(x_start, 1, ctx.Spectral.NumNodes);
            else
                t_local = (ctx.Spectral.Nodes + 1)/2 * dt_try;
                X_guess = x_start + u_dot_prev * t_local';
            end
            
            % 3. 执行单次扫描 (One Sweep Cost)
            % 仅迭代一次，用最小成本获取波形趋势
            [X_probe, ~, ~] = obj.performISDCSweep(X_guess, x_start, ctx);
            
            % 4. 分析波形 (Analyze)
            I_probe = X_probe(end, :); % 获取电流向量 [1 x Nodes]
            
            % 计算二阶导数 (Curvature)
            % dI/dtau = I * D'
            dI_dtau = I_probe * ctx.Spectral.D';
            d2I_dtau2 = dI_dtau * ctx.Spectral.D';
            
            % 转换回物理时间曲率: d2I/dt2 = d2I/dtau2 * (2/dt)^2
            curvature = abs(d2I_dtau2 * (2/dt_try)^2);
            
            % 寻找最大曲率点
            % 忽略首尾节点（边界处导数计算通常最不稳定）
            nodes_range = 2:length(curvature)-1;
            if isempty(nodes_range)
                 max_k = 0;
            else
                [max_k, idx_local] = max(curvature(nodes_range)); 
                real_idx = nodes_range(idx_local);
            end
            
            % 5. 判断逻辑
            found = false;
            t_split = -1;
            
            if max_k > obj.CurvatureThreshold
                found = true;
                % 计算拐点物理时间
                tau_split = ctx.Spectral.Nodes(real_idx);
                t_split = t_start + (tau_split + 1)/2 * dt_try;
            end
        end
        
        function D = computeDifferentiationMatrix(~, nodes)
            % COMPUTEDIFFERENTIATIONMATRIX 计算 GLL 节点的谱微分矩阵 D
            N = length(nodes) - 1;
            D = zeros(N+1, N+1);
            
            L_vals = zeros(N+1, 1);
            for i = 1:N+1
                [L_vals(i), ~] = SpectralTimeUtils.legendre_poly(nodes(i), N);
            end
            
            for i = 1:N+1
                for j = 1:N+1
                    if i ~= j
                        D(i, j) = (L_vals(i) / L_vals(j)) / (nodes(i) - nodes(j));
                    else
                        if i == 1 
                            D(i, j) = -0.25 * N * (N + 1);
                        elseif i == N+1 
                            D(i, j) = 0.25 * N * (N + 1);
                        else
                            D(i, j) = 0.0;
                        end
                    end
                end
            end
        end
        
        % --- 复用核心 ISDC 逻辑 (源自 ISDCSolver) ---
        
        function [X_new, max_diff_norm, final_res_norm] = performISDCSweep(obj, X_k, x0, ctx)
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
                
                [u_next, res_norm, ~] = obj.solveImplicitStep(u_guess, dt_sub, RHS_Vector, nodes(m), ctx);
                
                X_new(:, m) = u_next;
                
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
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