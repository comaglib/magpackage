classdef LUSDCSolver < handle
    % LUSDCSolver 基于 LU分解优化的 SDC 求解器 (LU-SDC + Radau IIA)
    % 
    % 核心机制: LU-Decomposed SDC (LU-SDC) / Q-ImEx SDC
    %   1. 节点策略: 使用 Radau IIA 节点 (L-stable)，特别适合强刚性场路耦合问题。
    %   2. 预条件子: 对稠密谱积分矩阵 Q 进行 LU 分解 (Q = L * U)。
    %   3. 扫描策略: 使用 L 矩阵的下三角部分进行"三角化扫描"。
    %      这使得每次扫描的收敛速度由线性提升至接近二次收敛。
    %
    % 包含特性:
    %   - 完整的日志输出 (fprintf)
    %   - 空气区域奇异性自动修复 (Adaptive Regularization)
    %   - 严格的动量守恒 RHS 构造 (Momentum Consistent)
    
    properties
        Assembler       % Finite Element Assembler (有限元组装器)
        LinearSolver    % Linear Solver Wrapper (线性求解器)
        
        PolyOrder = 3        % 多项式阶数 P (建议 3, 对应 3 个 Radau 节点, 5阶精度)
        MaxSDCIters = 3      % SDC 扫描次数 (LU-SDC 收敛极快，通常 3 次足够)
        SDCTolerance = 1e-5  % SDC 收敛容差
    end
    
    methods
        function obj = LUSDCSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 400; % 大幅增加内存分配百分比
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
            fprintf('   LU-SDC Solver (Radau IIA, P=%d, L-Stable)\n', obj.PolyOrder);
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
            
            % 数据容器
            full_time_list = []; full_current_list = []; full_B_list = [];
            probeB_History = zeros(length(timeSteps), 1);
            currentHistory = zeros(length(timeSteps), 1);
            
            % --- 2. LU-SDC 核心预计算 (Pre-computation) ---
            fprintf('   [Init] Computing Radau IIA nodes and LU decomposition...\n');
            
            % [关键] 获取 Radau IIA 节点 (区间 (-1, 1])
            [radau_nodes, radau_weights] = obj.getRadauNodes(obj.PolyOrder);
            
            % 计算积分矩阵 Q_mat (大小 P x P)
            Q_mat = obj.computeIntegrationMatrix(radau_nodes);
            
            % [关键] 对 Q 进行 LU 分解
            [L_mat, U_mat] = lu(Q_mat);
            
            ctx.Spectral.Nodes = radau_nodes; 
            ctx.Spectral.Weights = radau_weights;
            ctx.Spectral.Q_mat = Q_mat;
            ctx.Spectral.L_mat = L_mat; 
            ctx.Spectral.U_mat = U_mat; 
            ctx.Spectral.NumNodes = length(radau_nodes); 
            
            % 绘图插值点
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
                
                fprintf('\n=== Slab %d/%d (dt=%.1e, LU-SDC) ===\n', slab_idx, length(timeSteps), dt_slab);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. Initialization
                % Radau 节点不包含 t=-1 (x_start)，所以 X_slab 第一列存储 x_start
                X_slab = zeros(ctx.numTotalDofs, ctx.Spectral.NumNodes + 1);
                X_slab(:, 1) = x_start; 
                
                % 平坦预测 (Flat Prediction)
                for k = 1:ctx.Spectral.NumNodes
                    X_slab(:, k+1) = x_start;
                end
                
                % 2. LU-SDC Sweep
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    fprintf('   [Sweep %d] ', sdc_iter);
                    
                    % 执行扫描
                    [X_slab, max_diff_norm, residual_norm] = obj.performLUSDCSweep(X_old, x_start, ctx);
                    
                    fprintf('Update=%.2e, Res=%.2e', max_diff_norm, residual_norm);
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' -> Converged.\n');
                        break;
                    else
                        fprintf('\n');
                    end
                end
                
                % --- High-Res Output ---
                t_plot_phys = t_start + (tau_plot + 1) / 2 * dt_slab;
                all_nodes = [-1; ctx.Spectral.Nodes]; % 包含起点用于插值
                
                % 电流插值
                I_vals = X_slab(end, :)'; 
                I_plot = obj.interpolatePolynomial(all_nodes, I_vals, tau_plot);
                
                % 磁密插值
                B_vals = zeros(length(all_nodes), 1);
                if ~isempty(postProc)
                    for k = 1:length(all_nodes)
                        B_vec = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_vals(k) = norm(B_vec);
                    end
                end
                B_plot = obj.interpolatePolynomial(all_nodes, B_vals, tau_plot);
                
                % 数据拼接
                if slab_idx == 1, idx_start = 1; else, idx_start = 2; end
                full_time_list = [full_time_list; t_plot_phys(idx_start:end)];
                full_current_list = [full_current_list; I_plot(idx_start:end)];
                full_B_list = [full_B_list; B_plot(idx_start:end)];
                
                % Update State
                x_start = X_slab(:, end); % Radau 最后一个节点即为 t_end
                globalTime = t_end;
                currentHistory(slab_idx) = x_start(end);
                if ~isempty(postProc), probeB_History(slab_idx) = B_vals(end); end
                
                if ~isempty(postProc)
                    fprintf('      [End] |B|=%.4f T, I=%.4f A\n', B_vals(end), full_current_list(end));
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
        
        function [X_new, max_diff_norm, final_res_norm] = performLUSDCSweep(obj, X_k, x0, ctx)
            % PERFORMLUSDCSWEEP 执行 LU-SDC 扫描 (Dimensionally Corrected)
            %
            % 方程: M * u_m - dt_eff * F(u_m) = RHS_Vector
            % RHS_Vector = (M * u_0) + (Accumulated_Impulse_L) + (Correction_Impulse_Q)
            % 所有项必须统一为 [Mass * State] 或 [Force * Time] 量纲 (动量/磁链)
            
            num_radau = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            dt_half = dt_slab / 2.0; 
            
            L_mat = ctx.Spectral.L_mat;
            Q_mat = ctx.Spectral.Q_mat;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            % 1. 计算上一轮解的物理力 F(X_k)
            %    F_k_radau 对应 Radau 节点 (t1...tP)
            F_k_radau = obj.evaluateF(X_k(:, 2:end), ctx); 
            
            % 2. 计算积分修正项 (Defect / Correction)
            %    I_full = Integral(F) dt.  单位: [Force * Time]
            I_full = (F_k_radau * Q_mat.') * dt_half; 
            
            % 3. 预计算初始动量项 (Initial Momentum: M * u0)
            %    这个基准项对所有节点都一样，是 SDC 的锚点
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma; 
            [w_idx, ~, w_val] = find(ctx.W_vec);
            
            M_u0 = M_sys * x0;
            I_0 = x0(end);
            term_L = ctx.circuit.L * I_0;
            term_W = ctx.W_vec' * x0; 
            M_u0(end) = M_u0(end) + (term_L + term_W) * Scale;
            
            % 准备 F_new 存储
            F_new_radau = zeros(size(F_k_radau));
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            % 4. 序列扫描 (从 m=1 到 P)
            for m = 1:num_radau
                
                % A. 有效时间步长
                eff_dt = L_mat(m,m) * dt_half;
                
                % B. 显式累积项 (Accumulated Impulse from Lower Triangular L)
                %    Sum_{j<m} (L_mj * dt * F_new_j)
                accum_impulse = zeros(ctx.numTotalDofs, 1);
                for j = 1 : m-1
                    if L_mat(m, j) ~= 0
                        accum_impulse = accum_impulse + (L_mat(m, j) * dt_half) * F_new_radau(:, j);
                    end
                end
                
                % C. 修正项 (Defect)
                %    Defect = I_full_m - (L * F_k)_m
                %    这一项代表了"上一轮迭代得到的高阶积分"与"当前步骤期望的低阶近似"之间的差
                L_times_Fk = (L_mat(m, :) * F_k_radau')' * dt_half;
                correction_impulse = I_full(:, m) - L_times_Fk;
                
                % D. 组装总右端项 (RHS_Vector)
                %    RHS = M*u0 + Accum + Correction
                %    物理含义: t_m 时刻的系统总动量/磁链
                RHS_Vector = M_u0 + accum_impulse + correction_impulse;
                
                % E. 求解
                %    Solve: M*u - eff_dt*F(u) = RHS_Vector
                u_guess = X_k(:, m+1); 
                
                [u_next, res_norm, ~] = obj.solveImplicitStep(u_guess, eff_dt, RHS_Vector, ctx.Spectral.Nodes(m), ctx);
                
                X_new(:, m+1) = u_next;
                
                % F. 更新力 F_new (用于下一个节点的显式累积)
                F_new_col = obj.evaluateF(u_next, ctx);
                F_new_radau(:, m) = F_new_col;
                
                % 统计
                diff = norm(u_next - X_k(:, m+1)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
            end
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_node, ctx)
            % 局部 Newton 求解器 (Robust Version with Adaptive Regularization)
            % 求解方程: M*u - dt*F(u) - RHS_Vector = 0
            
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_tot = M_sys + M_cir_diag + M_row_coup;
            
            final_res = 1e10; res_norm_init = 1.0;
            iter_count = 0;
            
            t_phys = ctx.t_start + (t_node + 1) * ctx.dt_slab / 2.0;
            V_src = ctx.circuit.V_source_func(t_phys);
            
            % [dt 保护] 防止极小步长导致除零，但不能太大改变物理
            if abs(dt) < 1e-12, dt = 1e-12; end
            
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
                neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
                
                % Residual = M*u + dt*(-Force) - RHS
                Res = M_tot * u_curr + dt * neg_Force - RHS_Vector;
                
                res_norm = norm(Res);
                if iter == 1, res_norm_init = max(res_norm, 1e-6); end
                rel_err = res_norm / res_norm_init;
                
                % 收敛判断
                if rel_err < 1e-6 || res_norm < 1e-9
                    final_res = res_norm; break; 
                end
                
                % 组装切线刚度 K_tan
                [ki, kj, kv] = find(K_fem_small);
                K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                K_tan = K_fem + ctx.C_AP + ctx.C_AP' + ...
                        sparse(w_idx, repmat(num_dofs, length(w_idx),1), -w_val, num_dofs, num_dofs) + ...
                        sparse(num_dofs, num_dofs, ctx.circuit.R*Scale, num_dofs, num_dofs);
                
                % Jacobian = M + dt * K_tan
                Jac = M_tot + dt * K_tan;
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(Jac, Res, ctx.fixedDofs);
                
                % --- [自适应正则化求解] ---
                solve_success = false;
                try
                    du = obj.LinearSolver.solve(J_bc, -Res_bc);
                    solve_success = true;
                catch
                    % 第一次尝试失败，进入正则化流程
                end
                
                if ~solve_success
                    % 计算自适应 Shift: 1e-6 * Max_Diagonal
                    % 为什么用最大对角元？因为电路耦合项可能达到 1e10，
                    % 如果只用平均值(1e5)去修补空气项(0)，会被数值误差淹没。
                    diag_ref = max(abs(diag(J_bc)));
                    if diag_ref < 1e-12, diag_ref = 1.0; end
                    shift_val = 1e-6 * diag_ref;
                    
                    if iter == 1
                        fprintf('          [Info] Matrix singular. Adaptive Shift: %.1e\n', shift_val);
                    end
                    
                    J_reg = J_bc + shift_val * speye(size(J_bc));
                    
                    try
                        du = obj.LinearSolver.solve(J_reg, -Res_bc);
                        solve_success = true;
                    catch
                        fprintf('          [Error] Solve failed even with regularization.\n');
                        u_new = u_curr; final_res = res_norm; return;
                    end
                end
                % -----------------------
                
                u_curr = u_curr + du;
            end
            u_new = u_curr;
            
            % 严格收敛警告: 如果残差依然巨大，说明物理无解或参数错误
            if final_res > 1e-2 && final_res/res_norm_init > 1e-1
                 fprintf('          [CRITICAL] Newton failed. Res=%.2e. Results may be invalid.\n', final_res);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            % 计算物理方程 F(u) (Batch Mode)
            num_cols = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_cols);
            
            C_tot = ctx.C_AP + ctx.C_AP';
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            
            for k = 1:num_cols
                x = X_mat(:, k);
                % 简化处理: 使用 Slab 中点时间近似计算源项
                t_approx = ctx.t_start + dt_half; 
                V_src = ctx.circuit.V_source_func(t_approx); 
                
                A_vec = x(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(R_fem)) = -R_fem; 
                f_vec = f_vec - C_tot * x;
                f_vec = f_vec + ctx.W_vec * x(end);
                
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                F_val(:, k) = f_vec;
            end
        end
        
        function [nodes, weights] = getRadauNodes(~, P)
            % 获取 Radau IIA 节点 (映射到 [-1, 1])
            if P == 3
                c = [(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
                nodes = 2 * c - 1;
                w = [(16+sqrt(6))/36; (16-sqrt(6))/36; 1/9];
                weights = 2 * w;
            elseif P == 2
                 c = [1/3; 1];
                 nodes = 2 * c - 1;
                 w = [3/4; 1/4];
                 weights = 2 * w;
            else
                error('Radau nodes for P=%d not implemented locally.', P);
            end
        end
        
        function Q = computeIntegrationMatrix(~, nodes)
            % 计算积分矩阵 Q (使用数值积分替代符号计算)
            N = length(nodes);
            Q = zeros(N, N);
            
            for i = 1:N 
                t_end = nodes(i);
                for j = 1:N 
                    % 使用 MATLAB 内置 integral 进行高精度数值积分
                    func_j = @(t) lagrange_basis(t, nodes, j);
                    Q(i, j) = integral(func_j, -1, t_end, 'ArrayValued', true);
                end
            end
            
            function val = lagrange_basis(t, all_nodes, j)
                val = ones(size(t));
                xj = all_nodes(j);
                for k = 1:length(all_nodes)
                    if k ~= j
                        xk = all_nodes(k);
                        val = val .* (t - xk) / (xj - xk);
                    end
                end
            end
        end
        
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
        
        function vals = interpolatePolynomial(~, nodes, values, t_targets)
             N = length(nodes); M = length(t_targets); vals = zeros(M,1);
             for k=1:M
                 t=t_targets(k); y=0;
                 for j=1:N
                     lj=1; for i=1:N, if i~=j, lj=lj*(t-nodes(i))/(nodes(j)-nodes(i)); end; end
                     y=y+values(j)*lj;
                 end
                 vals(k)=y;
             end
        end
    end
end