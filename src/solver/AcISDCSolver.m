classdef AcISDCSolver < handle
    % AcISDCSolver - Anderson Accelerated Backward Euler SDC Solver
    % 
    % 算法描述:
    %   本求解器专为解决强非线性、极度刚性的电磁瞬态问题（如变压器励磁涌流、短路）设计。
    %   它结合了以下两种核心技术：
    % 
    %   1. 基础积分器 (Base Integrator): Backward Euler (一阶后向欧拉)
    %      - 优势: 具有最大的数值阻尼 (L-stable)，能强力抑制高频振荡和 Gibbs 现象。
    %      - 作用: 在每个 SDC 子步中，作为内层求解器，保证非光滑过程的稳定性。
    % 
    %   2. 加速技术 (Acceleration): Anderson Acceleration (安德森加速)
    %      - 优势: 将不动点迭代的线性收敛转化为超线性收敛。
    %      - 作用: 利用过去 m 次 Sweep 的历史残差，构建最优线性组合，显著减少 SDC 
    %        外层扫描次数 (通常从 10-15 次减少到 3-5 次)，大幅提升计算效率。
    % 
    % 适用场景:
    %   - 变压器空载合闸（励磁涌流）
    %   - 铁磁材料深度饱和与磁滞特性计算
    %   - 包含尖峰波形的刚性电路耦合问题
    
    properties
        Assembler       % 有限元组装器 (负责刚度矩阵、质量矩阵计算)
        LinearSolver    % 线性求解器 (MUMPS 包装器)
        
        % --- SDC 基础参数 ---
        PolyOrder = 4        % 谱元阶数 P (决定了 Time Slab 内的 GLL 节点数 = P+1)
                             % 对于饱和问题，建议 P=4 (5个节点)，兼顾精度与稳定性。
                             
        MaxSDCIters = 10     % 最大扫描次数 (启用加速后，通常 5 次以内收敛)
        SDCTolerance = 1e-5  % SDC 收敛容差 (相对误差)
        
        % --- Anderson 加速参数 ---
        EnableAnderson = true;  % 开关: 是否启用安德森加速
        AndersonDepth = 3;      % 历史深度 m: 利用过去多少步的信息 (建议 2-5)
        AndersonMixing = 1.0;   % 混合因子 beta: 类似于松弛因子 (通常 1.0)
    end
    
    methods
        function obj = AcISDCSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            
            % 初始化线性求解器
            % 针对非对称矩阵 (电路耦合) 进行优化
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            
            % 增加内存预留比例 (400%)，防止在复杂网格下因 Fill-in 导致内存溢出
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
            % SOLVE 主求解流程
            %
            % 输入参数:
            %   space_A, space_P: 有限元空间对象 (Nedelec, Lagrange)
            %   matLib: 材料库对象 (包含 B-H 曲线)
            %   sigmaMap: 电导率映射
            %   circuitProps: 电路参数 (R, L, V_source)
            %   windingObj: 绕组几何信息
            %   timeSteps: 时间步长数组 (定义了 Time Slabs)
            %   fixedDofs_A/P: 边界条件索引
            %   init_State: 初始状态向量 (A, P, I)
            
            fprintf('\n==========================================================\n');
            if obj.EnableAnderson
                fprintf('   ISDC Solver (Backward Euler + Anderson Acceleration)\n');
            else
                fprintf('   ISDC Solver (Standard Backward Euler)\n');
            end
            fprintf('   PolyOrder=%d, Tol=%.1e, HistDepth=%d\n', ...
                obj.PolyOrder, obj.SDCTolerance, obj.AndersonDepth);
            fprintf('==========================================================\n');
            
            % --- 1. 初始化自由度与上下文 ---
            dofHandler = obj.Assembler.DofHandler;
            % 确保自由度已分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            % 构建计算上下文 (Context)，减少函数间参数传递
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; % +1 为电路电流自由度 I
            ctx.space_A = space_A; ctx.space_P = space_P;
            ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            % 初始化后处理探针 (如果有)
            postProc = []; 
            if ~isempty(probePoint), postProc = PostProcessor(obj.Assembler); end
            
            % --- 2. 准备谱元数据 (GLL Nodes) ---
            % 生成 Gauss-Lobatto-Legendre 节点和积分矩阵 Q
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            ctx.Spectral.Nodes = gll_nodes; 
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = SpectralTimeUtils.integration_matrix(gll_nodes); 
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % 绘图用的插值点 (用于生成平滑曲线)
            NumPlotPoints = 30; 
            tau_plot = linspace(-1, 1, NumPlotPoints)'; 
            
            % --- 3. 组装时不变矩阵 (Pre-assembly) ---
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % 导电率质量矩阵 M_sigma
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 规范势耦合矩阵 C_AP (用于 A-V 公式)
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 绕组耦合向量 W (FEM -> Circuit)
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            % 计算电路方程缩放因子 (Geo-balancing)
            % 用于改善矩阵条件数，平衡电磁场方程(1e-6)与电路方程(1e2)的数量级差异
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            fprintf('   [Init] Circuit Scale Factor = %.2e\n', ctx.CircuitRowScale);
            
            % 预组装全局常数质量矩阵 (包含 Sigma, L, W_coupling)
            Scale = ctx.CircuitRowScale;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_dofs = ctx.numTotalDofs;
            M_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.L * Scale, num_dofs, num_dofs);
            M_row_coup = sparse(repmat(num_dofs, length(w_idx), 1), w_idx, w_val * Scale, num_dofs, num_dofs);
            ctx.M_full_const = ctx.M_sigma + M_cir_diag + M_row_coup;
            
            % 处理边界条件索引
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            % 初始状态设置
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_curr(1:n) = init_State(1:n);
            end
            
            % 数据记录容器
            full_time_list = []; full_current_list = []; full_B_list = [];
            globalTime = 0;
            
            % ==================================================
            % 主时间循环 (Time Slab Loop)
            % ==================================================
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n--- Slab %d/%d | dt=%.2e | T=%.4f -> %.4f ---\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                % 1. 初始预测 (Initialization)
                % 简单的常数预测：假设整个 Slab 内的状态等于初始时刻状态
                X_slab = repmat(x_curr, 1, ctx.Spectral.NumNodes);
                
                % 清空 Anderson 加速的历史数据
                AA_Data.X_hist = [];
                AA_Data.G_hist = [];
                
                % 2. SDC 扫描循环 (Sweep Loop)
                for sdc_iter = 1:obj.MaxSDCIters
                    X_old = X_slab;
                    
                    % (A) 执行标准 ISDC 扫描 (Backward Euler)
                    % 计算一次物理上的修正，得到 X_sweep
                    [X_sweep, max_diff_norm] = obj.performISDCSweep(X_old, x_curr, ctx);
                    
                    % (B) 应用 Anderson 加速
                    if obj.EnableAnderson
                        % 将 X_old -> X_sweep 视为不动点迭代 G(x)
                        % 利用历史数据预测更优解 X_accel
                        [X_accel, AA_Data] = obj.applyAnderson(X_old, X_sweep, AA_Data);
                        X_slab = X_accel;
                        
                        % 重新计算加速后的收敛范数
                        diff_accel = norm(X_accel - X_old, 'fro') / (norm(X_accel, 'fro') + 1e-6);
                        
                        fprintf('   Sweep %d (AA): Update=%.2e (Raw=%.2e)\n', ...
                            sdc_iter, diff_accel, max_diff_norm);
                        
                        final_diff = diff_accel;
                    else
                        % 不加速
                        X_slab = X_sweep;
                        fprintf('   Sweep %d: Update=%.2e\n', sdc_iter, max_diff_norm);
                        final_diff = max_diff_norm;
                    end
                    
                    % (C) 收敛检查
                    if final_diff < obj.SDCTolerance
                        fprintf('      -> Converged.\n');
                        break;
                    end
                end
                
                % --- 3. 输出数据处理 ---
                % 计算绘图用的物理时间点
                t_plot_phys = t_start + (tau_plot + 1) / 2 * dt_slab;
                
                % 对解进行拉格朗日插值，获得平滑曲线
                I_nodes = X_slab(end, :)';
                I_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, I_nodes, tau_plot);
                
                % 如果有探针，计算并插值 B 场
                B_nodes = zeros(ctx.Spectral.NumNodes, 1);
                if ~isempty(postProc)
                    for k = 1:ctx.Spectral.NumNodes
                        B_val = postProc.probeB(X_slab(1:ctx.num_A, k), probePoint);
                        B_nodes(k) = norm(B_val);
                    end
                end
                B_plot = obj.interpolatePolynomial(ctx.Spectral.Nodes, B_nodes, tau_plot);
                
                % 拼接数据 (跳过第一个点以避免重复)
                if slab_idx == 1, idx_st = 1; else, idx_st = 2; end
                full_time_list = [full_time_list; t_plot_phys(idx_st:end)];
                full_current_list = [full_current_list; I_plot(idx_st:end)];
                full_B_list = [full_B_list; B_plot(idx_st:end)];
                
                % 更新下一 Slab 的初值
                x_curr = X_slab(:, end);
                globalTime = t_end;
                
                % 实时绘图回调
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
        
        function [X_accel, AA_Data] = applyAnderson(obj, X_old, X_new, AA_Data)
            % APPLYANDERSON 执行安德森加速逻辑
            %
            % 原理:
            %   定义残差 g(x) = G(x) - x = X_new - X_old
            %   目标是找到历史迭代的线性组合，使得残差范数最小化。
            %   x_{new} = x_k + beta * g_k - \sum [ gamma_i * (delta_x_i + delta_g_i) ]
            
            m = obj.AndersonDepth;
            beta = obj.AndersonMixing;
            
            % 将矩阵展平为向量，方便进行线性代数运算
            x_vec_old = X_old(:);
            x_vec_new = X_new(:);
            g_vec = x_vec_new - x_vec_old; % 当前残差
            
            % --- 1. 更新历史库 ---
            if isempty(AA_Data.X_hist)
                % 第一步，无历史，无法加速
                AA_Data.X_hist = x_vec_old;
                AA_Data.G_hist = g_vec;
                X_accel = X_new; 
                return;
            end
            
            % 追加当前步数据
            AA_Data.X_hist = [AA_Data.X_hist, x_vec_old];
            AA_Data.G_hist = [AA_Data.G_hist, g_vec];
            
            % 限制历史深度 (滑动窗口)
            if size(AA_Data.X_hist, 2) > m
                AA_Data.X_hist(:, 1) = [];
                AA_Data.G_hist(:, 1) = [];
            end
            
            mk = size(AA_Data.G_hist, 2);
            if mk < 2
                 X_accel = X_new; % 至少需要两步才能进行混合
                 return;
            end
            
            % --- 2. 求解最小二乘优化问题 ---
            % 目标: 最小化 || g_k - \sum gamma_i * (g_{k-i} - g_k) ||
            % 等价于在残差空间寻找最优组合系数
            
            current_g = g_vec;
            DG = zeros(length(g_vec), mk-1);
            
            % 构建差分矩阵 DG = [g_k - g_{k-1}, ..., g_k - g_{k-m}]
            for i = 1:mk-1
                DG(:, i) = current_g - AA_Data.G_hist(:, end-i);
            end
            
            % 求解系数 gamma (QR分解或直接除法)
            gamma = DG \ current_g; 
            
            % --- 3. 计算加速后的解 ---
            DX = zeros(length(x_vec_old), mk-1);
            for i = 1:mk-1
                DX(:, i) = x_vec_old - AA_Data.X_hist(:, end-i);
            end
            
            % 应用修正公式
            correction = (DX + DG) * gamma;
            x_acc_vec = x_vec_old + beta * current_g - correction;
            
            % 恢复形状
            X_accel = reshape(x_acc_vec, size(X_new));
        end
        
        function [X_new, max_diff_norm] = performISDCSweep(obj, X_k, x0, ctx)
            % PERFORMISDCSWEEP 执行一次标准的 ISDC 扫描 (Backward Euler)
            %
            % 流程:
            %   1. 计算全时间片的高阶谱积分 I_spec。
            %   2. 逐个子步求解: M(u_m - u_{m-1}) - dt*F(u_m) = I_spec - I_low
            
            num_nodes = ctx.Spectral.NumNodes;
            dt_slab = ctx.dt_slab;
            Q = ctx.Spectral.Q;
            nodes = ctx.Spectral.Nodes;
            
            X_new = zeros(size(X_k));
            X_new(:, 1) = x0; 
            
            % 1. 计算上一轮解的物理力 F(u^k) 并积分
            F_k = obj.evaluateF_Batch(X_k, ctx); 
            dt_half = dt_slab / 2.0;
            I_spec = (F_k * Q.') * dt_half; % 谱积分 [DoFs x Nodes]
            
            max_diff_norm = 0;
            
            % 2. 逐节点推进 (Sequential Sweep)
            for m = 2:num_nodes
                % 计算当前子区间的谱积分增量
                int_corr = I_spec(:, m) - I_spec(:, m-1);
                
                u_prev_new = X_new(:, m-1);
                F_k_curr   = F_k(:, m);
                h_curr = (nodes(m) - nodes(m-1)) * dt_half; % 子步长
                
                % --- 构造隐式方程右端项 (RHS) ---
                % SDC 修正公式 (Backward Euler):
                % M * u_m - h * F(u_m) = M * u_{m-1} + int_corr - h * F(u_m^k)
                % 
                % 解释:
                %   M * u_{m-1}: 来自 BE 的显式部分
                %   int_corr: 高阶积分带来的"正确"物理量变化
                %   - h * F(u_m^k): 减去 BE 的低阶积分近似，避免重复计算
                
                M_u_prev = ctx.M_full_const * u_prev_new;
                RHS = M_u_prev + int_corr - h_curr * F_k_curr;
                
                % --- 求解局部非线性方程 ---
                dt_eff = h_curr;
                
                % 初值猜测 (Rolling vs Refining)
                dist = norm(X_k(:,m) - X_k(:,m-1));
                if dist < 1e-9, u_guess = u_prev_new; else, u_guess = X_k(:,m); end
                
                t_phys = ctx.t_start + (nodes(m)+1)*dt_half;
                
                % 调用牛顿求解器
                [u_next, iter_count, ~] = obj.solveImplicitStep(u_guess, dt_eff, RHS, t_phys, ctx);
                
                X_new(:, m) = u_next;
                
                % 记录相对变化
                diff = norm(u_next - X_k(:, m)) / (norm(u_next) + 1e-6);
                max_diff_norm = max(max_diff_norm, diff);
                
                fprintf('        Node %d Done. Iters=%d, Diff=%.2e\n', m, iter_count, max_diff_norm);
            end
        end
        
        function F_val = evaluateF_Batch(obj, X_mat, ctx)
            % 批量评估物理力 F (用于积分计算)
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
            % 单点评估物理力 F(u, t)
            % 方程: M u' = F(u)
            % F 包含: -Stiffness*A, -Coupling*P, +Winding*I, +VoltageSrc
            
            num_dofs = ctx.numTotalDofs;
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            % 1. 非线性磁阻力 -K(nu)*A
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            
            f_vec = sparse(num_dofs, 1);
            f_vec(1:length(R_fem)) = -R_fem; 
            
            % 2. 耦合项与源项
            f_vec = f_vec - (ctx.C_AP + ctx.C_AP') * x;
            f_vec = f_vec + ctx.W_vec * I_val; 
            
            % 3. 电路方程力项 (V - R*I)
            V_src = ctx.circuit.V_source_func(t_phys);
            f_vec(end) = (V_src - ctx.circuit.R * I_val) * ctx.CircuitRowScale;
        end
        
        function [u_new, iter, final_res] = solveImplicitStep(obj, u_guess, dt, RHS_Vector, t_phys, ctx)
            % SOLVEIMPLICITSTEP 求解局部非线性方程 (Damped Newton + Stagnation Check)
            % 求解: M*u - dt*F(u) = RHS
            % 特性: 集成了回溯线搜索 (防止发散) 和 停滞检测 (防止死循环)
            
            u_curr = u_guess;
            num_dofs = ctx.numTotalDofs;
            Scale = ctx.CircuitRowScale;
            
            % 缓存常数矩阵部分 (加速计算)
            [w_idx, ~, w_val] = find(ctx.W_vec);
            K_col_coup = sparse(w_idx, repmat(num_dofs, length(w_idx), 1), -w_val, num_dofs, num_dofs);
            K_cir_diag = sparse(num_dofs, num_dofs, ctx.circuit.R * Scale, num_dofs, num_dofs);
            K_const_part = ctx.C_AP + ctx.C_AP' + K_col_coup + K_cir_diag;
            M_full = ctx.M_full_const;
            
            final_res = 1e10; 
            res_norm_init = 1.0; 
            res_norm_prev = 1e10;
            
            MAX_ITERS = 15; 
            STAGNATION_TOL = 1e-2; % 停滞容差
            MIN_ALPHA = 1e-2;      % 线搜索最小步长因子
            
            % 预先计算初始残差和刚度矩阵 (Iter 1 准备)
            A_vec = u_curr(1:ctx.num_A);
            [K_fem_small, R_fem_base] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            Res = obj.computeResidual(u_curr, R_fem_base, dt, RHS_Vector, t_phys, ctx, M_full, K_const_part);
            res_norm = norm(Res);
            
            for iter = 1:MAX_ITERS
                % --- 收敛判据 ---
                if iter == 1, res_norm_init = res_norm; if res_norm_init < 1e-20, res_norm_init=1; end; end
                
                % 1. 绝对/相对误差收敛
                if res_norm / res_norm_init < 1e-5 || res_norm < 1e-8
                    final_res = res_norm; 
                    % fprintf('          Newton Converged: Rel=%.2e\n', res_norm/res_norm_init);
                    break; 
                end
                
                % 2. 停滞检测 (用户要求保留，防止在数值底噪附近死循环)
                if iter > 1
                    improvement = (res_norm_prev - res_norm) / res_norm_prev;
                    if improvement < STAGNATION_TOL
                        % 认为已尽力，停止迭代
                        final_res = res_norm; 
                        % fprintf('          [Stagnation] Improvement %.2e < Tol\n', improvement);
                        break;
                    end
                end
                res_norm_prev = res_norm;
                
                % --- 牛顿步计算 ---
                
                % 1. 求解牛顿方向: Jac * du = -Res
                [ki, kj, kv] = find(K_fem_small);
                K_fem = sparse(ki, kj, kv, num_dofs, num_dofs);
                Jac = M_full + dt * (K_fem + K_const_part);
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(Jac, Res, ctx.fixedDofs);
                du = obj.LinearSolver.solve(J_bc, -Res_bc);
                
                % 2. 回溯线搜索 (Backtracking Line Search)
                alpha = 1.0;
                u_next = u_curr;
                Res_next = Res;
                valid_step = false;
                
                % 缓存旧的 K_fem，如果线搜索回退，可能不需要重新组装
                % 但为了简单起见，我们在接受步长后再组装下一次的 K
                
                while alpha > MIN_ALPHA
                    u_try = u_curr + alpha * du;
                    
                    % 仅计算残差 (False: 不组装刚度矩阵，只计算力向量 R_fem)
                    A_vec_try = u_try(1:ctx.num_A);
                    [~, R_fem_try] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec_try, ctx.matLib, false);
                    
                    Res_try = obj.computeResidual(u_try, R_fem_try, dt, RHS_Vector, t_phys, ctx, M_full, K_const_part);
                    res_norm_try = norm(Res_try);
                    
                    % 单调下降判据 (Monotonicity Check)
                    if res_norm_try < res_norm
                        % 接受步长
                        u_next = u_try;
                        Res_next = Res_try;
                        res_norm = res_norm_try; % 更新 res_norm 供下一次循环检查
                        
                        % **关键**: 线搜索成功后，立即为下一次迭代组装 K_fem (True)
                        % 这样回到循环开头时，K_fem_small 已经是基于 u_next 的最新值
                        [K_fem_small, ~] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec_try, ctx.matLib, true);
                        
                        valid_step = true;
                        break;
                    else
                        % 拒绝，减半步长
                        alpha = alpha * 0.3;
                    end
                end
                
                if ~valid_step
                    % 如果线搜索彻底失败 (极罕见)，通常意味着 Jacobian 病态或陷入局部极小
                    % 强制接受全步长并退出，或者保留上一步结果退出
                    % 这里选择强制更新一次并退出，避免死锁
                    u_curr = u_curr + du;
                    % fprintf('          [LineSearch Failed] Step forced.\n');
                    final_res = res_norm; % 这里的 res_norm 还是旧的，但也只能这样了
                    break;
                else
                    u_curr = u_next;
                    Res = Res_next;
                    % fprintf('(alpha=%.4f) ', alpha);
                end
            end
            
            u_new = u_curr;
            
            % if final_res / res_norm_init > 1e-1 && final_res > 1.0
            %     fprintf('          [Warn] Local Newton stalled at high residual: Rel=%.2e\n', final_res / res_norm_init);
            % end
        end
        
        function Res = computeResidual(obj, u, R_fem, dt, RHS_Vector, t_phys, ctx, M_mat, K_const)
            % 辅助函数: 计算残差 M*u - dt*F(u) - RHS
            % F(u) = -R_fem - K_const*u + V_src
            % R_fem: 非线性磁阻力向量 K(nu)*A (由 assembleJacobian 计算得到)
            
            % 1. 构造 -Force 向量
            % 场方程部分: K(u)*u -> R_fem (外部已传入)
            neg_Force = sparse(size(u,1), 1);
            neg_Force(1:length(R_fem)) = R_fem;
            
            % 加上线性刚度力 (K_const * u)
            % 注意符号: 物理方程是 M u' + K u = Src
            % SDC方程是 M u' = F(u) => F(u) = Src - K u
            % 这里的 neg_Force 对应 -F(u) = K u - Src
            % R_fem 已经是 K_nonlin * A
            neg_Force = neg_Force + K_const * u;
            
            % 修正电路源项 (V_source)
            % 电路方程: L I' + R I = V => L I' = V - R I
            % F_circuit = (V - R I) * Scale
            % neg_F_circuit = (R I - V) * Scale
            % K_const 中已经包含了 R*I 的项 (即 K_cir_diag * I)
            % 所以这里只需要减去 V * Scale
            
            V_src = ctx.circuit.V_source_func(t_phys);
            neg_Force(end) = neg_Force(end) - V_src * ctx.CircuitRowScale;
            
            % 2. 计算最终残差
            % Res = M*u + dt*neg_Force - RHS
            Res = M_mat * u + dt * neg_Force - RHS_Vector;
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
             % 电路方程缩放因子计算 (保持不变)
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
            % 获取极限磁导率分布
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
             % Lagrange 插值 (用于后处理平滑曲线)
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