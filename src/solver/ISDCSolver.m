classdef ISDCSolver < handle
    % ISDCSolver 基于非线性序列扫描的 SDC 求解器 (v12.0 - Nonlinear ISDC + Full Logging)
    % 
    % 核心机制: Implicit SDC (ISDC)
    %   1. 不做全局线性化预测，而是沿着 GLL 节点逐个求解非线性方程 (Sequential Sweeping)。
    %   2. 每个子步本质上是一个 "带高阶源项修正的 Backward Euler" (L-stable)。
    %   3. 这种方法对变压器饱和等刚性问题极其稳健 (Stiff & Nonlinear)。
    
    properties
        Assembler       % Finite Element Assembler (有限元组装器，负责空间离散)
        LinearSolver    % Linear Solver Wrapper (线性求解器，如 MUMPS)
        
        PolyOrder = 3        % 多项式阶数 P (建议 3 或 4)。决定了每个 Time Slab 内有 P+1 个 GLL 节点。
        MaxSDCIters = 5      % SDC 扫描次数 K (ISDC 收敛很快，通常 5 次足够将误差降至底噪)。
        SDCTolerance = 1e-4  % SDC 收敛容差 (用于提前退出扫描循环)。
    end
    
    methods
        function obj = ISDCSolver(assembler)
            % 构造函数：初始化组装器和线性求解器配置
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; % 非对称矩阵 (因为包含电路耦合项)
            obj.LinearSolver.MumpsICNTL.i14 = 300; % 增加工作内存百分比，防止溢出
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
            % 输入:
            %   space_A/P: 有限元空间对象
            %   matLib/sigmaMap: 材料属性库和电导率映射
            %   circuitProps/windingObj: 电路参数和绕组几何信息
            %   timeSteps: 时间步长数组 (每个元素代表一个 Time Slab 的长度)
            %   fixedDofs: 边界条件
            %   init_State: 初始状态向量
            %   monitorHandle: 用于实时绘图的回调函数
            %   probePoint: 磁密探测点坐标
            
            fprintf('==================================================\n');
            fprintf('   SDC Solver (Nonlinear ISDC, P=%d)\n', obj.PolyOrder);
            fprintf('==================================================\n');
            
            % --- 1. 初始化与预计算 (Initialization & Pre-computation) ---
            % 分配自由度 (DoFs)
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            % 构建上下文结构体 (Context)，用于在各函数间传递参数
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; % +1 是为了电路电流自由度 I
            ctx.space_A = space_A; ctx.space_P = space_P;
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            % 初始化后处理探针
            postProc = []; 
            if ~isempty(probePoint)
                postProc = PostProcessor(obj.Assembler);
            end
            
            % [Data Collection] 全量数据容器 (High-Res Interpolated)
            % 用于存储插值后的平滑波形数据
            full_time_list = [];
            full_current_list = [];
            full_B_list = [];
            
            % 简略历史 (仅记录每个 Slab 的终点值)
            probeB_History = zeros(length(timeSteps), 1);
            currentHistory = zeros(length(timeSteps), 1);
            
            % --- 2. SDC 谱积分矩阵计算 (Spectral Integration Matrix) ---
            fprintf('   [Init] Computing GLL nodes and integration matrix...\n');
            % 计算 Gauss-Lobatto-Legendre 节点和权重
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            % 计算积分矩阵 Q (Q_ij = integral of Lagrange basis j from -1 to node i)
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat; 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % [Plotting] 预计算插值点 (用于稠密输出)
            % 在每个 Time Slab 内部生成 50 个均匀分布点，用于绘制平滑曲线
            NumPlotPoints = 50;
            tau_plot = linspace(-1, 1, NumPlotPoints)'; % [-1, 1] 均匀分布
            
            % --- 3. 组装不变矩阵 (Invariant Matrices Assembly) ---
            fprintf('   [Init] Assembling invariant matrices...\n');
            % 质量矩阵 M_sigma (加权电导率)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 耦合矩阵 C (如果有 Lagrange 乘子)
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 绕组耦合向量 W (FEM -> Circuit)
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            % 电路方程缩放因子 (Geo-balancing)
            % 用于平衡场方程(1e-6量级)和电路方程(1e2量级)的数值差异，改善条件数
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', ctx.CircuitRowScale);
            
            % 处理 Dirichlet 边界条件
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            % 设置初始状态向量 x_0
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
            globalTime = 0;
            
            % ==================================================
            % 主时间循环 (Time Slab Loop)
            % 每个循环处理一个大的时间片 [t_n, t_n + dt_slab]
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n=== Slab %d/%d (dt=%.1e, T=%.4f -> %.4f) ===\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. Initialization (初值填充)
                % 将整个 Slab 内所有节点的初值设为上一个 Slab 的终值 (Constant Prediction)
                X_slab = repmat(x_start, 1, ctx.Spectral.NumNodes);
                
                % 2. SDC Sweep (非线性序列扫描)
                % 核心迭代过程：反复"扫描"时间片，提升精度
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    % 调用核心扫描函数
                    [X_slab, max_diff_norm, residual_norm] = obj.performISDCSweep(X_old, x_start, ctx);
                    
                    % 打印收敛信息
                    if max_diff_norm < 1e-4
                        fprintf('   Iter %d: Update=%.2e, Res=%.2e', sdc_iter, max_diff_norm, residual_norm);
                    else
                        fprintf('   Iter %d: Update=%.4f, Res=%.4f', sdc_iter, max_diff_norm, residual_norm);
                    end
                    
                    % 检查收敛条件
                    if max_diff_norm < obj.SDCTolerance
                        fprintf(' -> Converged.\n');
                        break;
                    else
                        fprintf('\n');
                    end
                end
                
                % --- [High-Res Data Collection] ---
                % 使用谱系数 (Lagrange Basis) 插值得到平滑曲线 (Dense Output)
                
                % 1. 计算插值点的物理时间
                t_plot_phys = t_start + (tau_plot + 1) / 2 * dt_slab;
                
                % 2. 插值电流 (X_slab 的最后一行是电流 I)
                I_nodes = X_slab(end, :)';
                I_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, I_nodes, tau_plot);
                
                % 3. 插值磁密 (先在 GLL 节点计算 B，再插值 B 的模值)
                B_nodes = zeros(ctx.Spectral.NumNodes, 1);
                if ~isempty(postProc)
                    for k = 1:ctx.Spectral.NumNodes
                        B_val = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_nodes(k) = norm(B_val);
                    end
                end
                B_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, B_nodes, tau_plot);
                
                % 4. 拼接数据 (跳过第一个点以避免与上一 Slab 重复，除非是第一个 Slab)
                if slab_idx == 1
                    idx_start = 1;
                else
                    idx_start = 2;
                end
                
                full_time_list = [full_time_list; t_plot_phys(idx_start:end)];
                full_current_list = [full_current_list; I_plot(idx_start:end)];
                full_B_list = [full_B_list; B_plot(idx_start:end)];
                
                % --- Update State (更新状态，准备下一 Slab) ---
                x_start = X_slab(:, end); % 取 Slab 最后一个节点作为下一起点
                globalTime = t_end;
                currentHistory(slab_idx) = x_start(end);
                if ~isempty(postProc), probeB_History(slab_idx) = B_nodes(end); end
                
                if ~isempty(postProc)
                    fprintf('      [End] |B|=%.4f T, I=%.4f A\n', B_nodes(end), full_current_list(end));
                end
                
                % 实时绘图
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
            
            % 返回全波形数据
            info.Time_Full = full_time_list;
            info.Current_Full = full_current_list;
            info.ProbeB_Full = full_B_list;
        end
    end
    
    methods (Access = private)
        
        function [X_new, max_diff_norm, final_res_norm] = performISDCSweep(obj, X_k, x0, ctx)
            % PERFORMISDCSWEEP 执行一次隐式谱延迟修正扫描 (One ISDC Sweep)
            %
            % 原理: 
            %   求解方程: u_{m} - u_{m-1} = \int_{t_{m-1}}^{t_m} F(u) dt
            %   ISDC迭代格式: 
            %   M*(u_{m}^{k+1} - u_{m-1}^{k+1}) - dt*F(u_{m}^{k+1}) = -dt*F(u_{m}^k) + Integral_Correction
            %   即: Implicit_Euler(u^{k+1}) = Explicit_Terms(u^k) + Spectral_Integral_Error
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; % 起点固定 (IC)
            
            % 1. 计算上一轮解 X_k 在所有节点处的物理力 F(X_k)
            %    这是计算高阶积分修正项的基础
            % fprintf('      - Calculating High-Order Integral Source...\n');
            F_k = obj.evaluateF(X_k, ctx); 
            
            % 2. 计算全段积分项 I_accum = Integral_{0}^{t} F(u^k)
            dt_half = dt_slab / 2.0;
            I_accum = (F_k * Q.') * dt_half;
            
            max_diff_norm = 0;
            final_res_norm = 0;
            
            % 3. 序列扫描 (Sequential Sweep from m=2 to N)
            for m = 2:num_nodes
                tau_curr = nodes(m);
                tau_prev = nodes(m-1);
                dt_sub = (tau_curr - tau_prev) * dt_half; % 当前子步长
                
                % 准备隐式方程的右端项 (RHS)
                u_prev_new = X_new(:, m-1);   % 刚算出来的最新前一点值
                f_curr_old = F_k(:, m);       % 上一轮迭代在本点的力
                
                % 积分修正增量: (Integral_{0}^{m} - Integral_{0}^{m-1}) = Integral_{m-1}^{m}
                int_delta = I_accum(:, m) - I_accum(:, m-1);
                
                % 构造常数右端向量 RHS_Vector
                % 公式: M*u_{new} - dt*F(u_{new}) = M*u_{prev} + int_delta - dt*F(u_{old})
                % 这里的 RHS_Vector 对应等式右边
                
                Scale = ctx.CircuitRowScale;
                M_sys = ctx.M_sigma; 
                [w_idx, ~, w_val] = find(ctx.W_vec);
                
                % 计算 M * u_prev
                M_u_prev = M_sys * u_prev_new;
                I_prev = u_prev_new(end);
                % 加上电路部分的 Mass (L*I 和 W'*A)
                term_L = ctx.circuit.L * I_prev;
                term_W = ctx.W_vec' * u_prev_new;
                M_u_prev(end) = M_u_prev(end) + (term_L + term_W) * Scale;
                
                RHS_Vector = M_u_prev + int_delta - dt_sub * f_curr_old;
                
                % =========================================================
                % [关键改进] 智能初值猜测 (Smart Initial Guess)
                % =========================================================
                % 策略:
                %   1. 如果是 Sweep 1 (解还是平的): 利用时间连续性，猜测 u_m = u_{m-1} (Rolling Guess)。
                %      这能极大加速涌流上升沿的求解 (避免从0猜起)。
                %   2. 如果是 Sweep > 1: 使用上一轮迭代的值 u_m^k (Refining Guess)。
                
                dist_history = norm(X_k(:, m) - X_k(:, m-1));
                
                if dist_history < 1e-9 
                    % Sweep 1: Rolling Guess
                    u_guess = u_prev_new;
                    % fprintf('        [Guess] Using Prev Node (Rolling)\n');
                else
                    % Sweep > 1: Refining Guess
                    u_guess = X_k(:, m);
                    % fprintf('        [Guess] Using Prev Sweep (Refining)\n');
                end
                
                % 执行局部牛顿迭代 (Local Newton Solve)
                [u_next, res_norm, iter_count] = obj.solveImplicitStep(u_guess, dt_sub, RHS_Vector, nodes(m), ctx);
                
                X_new(:, m) = u_next;
                
                % 记录相对更新量
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                final_res_norm = max(final_res_norm, res_norm);
                
                fprintf('        -> Node %d Solved: %d Iters, Res=%.2e\n', m, iter_count, res_norm);
            end
        end
        
        function F_val = evaluateF(obj, X_mat, ctx)
            % EVALUATEF 计算物理方程的右端力项 F(u)
            % 方程形式: M * du/dt = F(u)
            % F(u) = -K(u)*u - C*u + W*I + V_source (Field)
            %        V_source - R*I               (Circuit)
            
            num_nodes = size(X_mat, 2);
            num_dofs = ctx.numTotalDofs;
            F_val = zeros(num_dofs, num_nodes);
            
            C_tot = ctx.C_AP + ctx.C_AP';
            Scale = ctx.CircuitRowScale;
            dt_half = ctx.dt_slab / 2.0;
            
            for m = 1:num_nodes
                x = X_mat(:, m);
                t = ctx.t_start + (ctx.Spectral.Nodes(m) + 1) * dt_half; % 物理时间
                
                % 计算非线性刚度项 K_fem * A
                A_vec = x(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                % 组装场方程力向量
                % Force = -K*x - C*x + W*I
                f_vec = sparse(num_dofs, 1);
                f_vec(1:length(R_fem)) = -R_fem; 
                f_vec = f_vec - C_tot * x;
                f_vec = f_vec + ctx.W_vec * x(end); % 电路对场的耦合 (+W*I)
                
                % 组装电路方程力向量
                % Circuit Eq: V - R*I - d/dt(...) = 0 => M_cir * dI/dt = V - R*I
                V_src = ctx.circuit.V_source_func(t);
                I_val = x(end);
                f_cir = (V_src - ctx.circuit.R * I_val) * Scale;
                f_vec(end) = f_cir;
                
                F_val(:, m) = f_vec;
            end
        end
        
        function [u_new, final_res, iter_count] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_node, ctx)
            % SOLVEIMPLICITSTEP 求解局部非线性方程 (Local Nonlinear Solver)
            % 方程: M*u - dt*F(u) - RHS_Vector = 0
            % 方法: Newton-Raphson 迭代
            % [v12.3 Fix] 纯相对误差控制，防止在数值底噪处死循环
            
            u_curr = u_guess;
            Scale = ctx.CircuitRowScale;
            M_sys = ctx.M_sigma;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            
            % 组装总质量矩阵 M_tot (Field Mass + Circuit Inductance + Coupling)
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            M_tot = M_sys + M_cir_diag + M_row_coup;
            
            final_res = 1e10;
            res_norm_init = 1.0; % 初始残差
            res_norm_prev = 1e10;
            iter_count = 0;
            
            % --- 牛顿迭代收敛参数 ---
            % 只要残差相比初始值下降了 6 个数量级，即认为物理收敛
            REL_TOL = 1e-6;      
            
            % 停滞容差: 如果某一步改善幅度小于 0.1%，说明到了数值底噪，强制停止
            STAG_TOL = 1e-3;     
            
            MAX_ITERS = 15;
            
            for iter = 1:MAX_ITERS
                iter_count = iter;
                
                % 1. 组装当前状态的 Jacobian 和 Residual
                A_vec = u_curr(1:ctx.num_A);
                [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                
                % 计算 -Force 项 (用于 Residual)
                % -Force = K*u + C*u - W*I + R*I - V
                neg_Force = sparse(num_dofs, 1);
                neg_Force(1:length(R_fem_base)) = R_fem_base;
                neg_Force = neg_Force + (ctx.C_AP + ctx.C_AP') * u_curr;
                
                I_val = u_curr(end);
                neg_Force = neg_Force - ctx.W_vec * I_val; 
                
                t_phys = ctx.t_start + (t_node + 1) * ctx.dt_slab / 2.0; 
                V_src = ctx.circuit.V_source_func(t_phys);
                neg_Force(end) = (ctx.circuit.R * I_val - V_src) * Scale;
                
                % Residual = M*u + dt*(-Force) - RHS
                Res = M_tot * u_curr + dt * neg_Force - RHS_Vector;
                
                res_norm = norm(Res);
                
                % --- 收敛判断逻辑 ---
                
                % 记录初始残差
                if iter == 1
                    res_norm_init = res_norm;
                    if res_norm_init < 1e-20, res_norm_init = 1.0; end % 防止除零
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
                
                % --- 2. 组装 Jacobian Matrix ---
                % J = M + dt * d(-Force)/du
                [ki, kj, kv] = find(K_fem_small);
                K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                
                K_base = K_fem + ctx.C_AP + ctx.C_AP';
                K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
                K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
                
                K_tan_tot = K_base + K_col_coup + K_cir_diag;
                Jac = M_tot + dt * K_tan_tot;
                
                % --- 3. 求解线性系统 ---
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
            % 计算电路方程的缩放因子 (Geo-Balance Heuristic)
            % 目的: 使电路方程的对角元与 FEM 方程的对角元在数值上接近，减少条件数。
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
            % 获取极限磁导率分布 (用于预条件子或缩放计算)
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
                        val = nu0; % 空气/饱和
                    else
                        val = nu0 / 1000.0; % 高导磁
                    end
                end
                nu_vec(meshTags == tag) = val;
            end
        end
        
        function vals_interp = interpolatePolynomial(~, nodes, values, target_nodes)
            % INTERPOLATEPOLYNOMIAL 拉格朗日多项式插值 (Lagrange Interpolation)
            % 作用: 基于 GLL 节点上的离散解，重构出时间片内连续的多项式解。
            % 输入:
            %   nodes: 已知节点 (GLL points), [N x 1]
            %   values: 已知节点上的值, [N x 1]
            %   target_nodes: 需要插值的点 (Plot points), [M x 1]
            % 输出:
            %   vals_interp: 插值结果 [M x 1]
            
            N = length(nodes);
            M = length(target_nodes);
            vals_interp = zeros(M, 1);
            
            % 对每个目标点进行求值
            for k = 1:M
                t = target_nodes(k);
                y = 0;
                
                % Lagrange Basis Sum: L(t) = sum( y_j * l_j(t) )
                for j = 1:N
                    % Compute l_j(t) = product( (t - x_i) / (x_j - x_i) )
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