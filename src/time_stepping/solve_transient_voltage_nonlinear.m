function [Solution, Info] = solve_transient_voltage_nonlinear(Model, Coil, CircuitProps, TimeParams)
% SOLVE_TRANSIENT_VOLTAGE_NONLINEAR 非线性瞬态场路耦合求解器 (Picard迭代版)
%
% 算法: 
%   采用定点迭代 (Fixed-Point / Picard Iteration) + 割线刚度矩阵 (Secant Stiffness)。
%   相比牛顿法，此方法在处理磁饱和与场路耦合时具有极高的数值稳定性。
%
% 特性:
%   1. 自适应阻尼松弛 (Adaptive Damping) 加速收敛并防止震荡。
%   2. MUMPS 线性求解器优先 (自动处理缩放)。
%
% 输入:
%   Model        - 有限元模型结构体
%   Coil         - 线圈对象 (含 Turns, P1, P2)
%   CircuitProps - 电路参数 (.R, .L_leak, .V_func)
%   TimeParams   - 时间参数 (.dt, .N_steps)
%
% 输出:
%   Solution     - 结果结构体 (.A, .I_history)
%   Info         - 运行信息 (.FinalTime)

    fprintf('==============================================\n');
    fprintf('   Nonlinear Transient Solver (Voltage Driven)\n');
    fprintf('   Method: Picard Iteration (Secant Stiffness)\n');
    fprintf('==============================================\n');

    dt = TimeParams.dt;
    N_steps = TimeParams.N_steps;
    V_func = CircuitProps.V_func;
    
    % 1. 准备数据
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    numFieldDoFs = numEdges + numActiveV;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 恒定矩阵
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    M_reg = assemble_mass_matrix(Model, 'All');
    
    % 绕组向量
    if isfield(Coil, 'Turns'), N_turns = Coil.Turns; else, N_turns = 1.0; end
    fprintf('  - Coil Turns: %d\n', N_turns);
    
    fprintf('正在组装绕组耦合向量 (Winding Vector)...\n');
    W_vec_A = assemble_winding_vector(Model, Coil);
    W_vec_A = W_vec_A * N_turns; 
    W_vec = [W_vec_A; sparse(numActiveV, 1)];
    
    % 感应源单位向量
    CoilUnit = Coil; CoilUnit.I(:) = 1.0 * N_turns;
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    
    % 2. 初始化
    x_curr = zeros(numFieldDoFs + 1, 1); % [A; V; I]
    x_prev = x_curr;
    I_prev = 0;
    A_prev = zeros(numEdges, 1);
    TimeHist.I = zeros(N_steps, 1);
    
    % 求解配置
    MaxIter = 30; Tol = 1e-4;
    if isfield(Model, 'Solver') && isfield(Model.Solver, 'Nonlinear')
        if isfield(Model.Solver.Nonlinear, 'MaxIter'), MaxIter=Model.Solver.Nonlinear.MaxIter; end
        if isfield(Model.Solver.Nonlinear, 'Tolerance'), Tol=Model.Solver.Nonlinear.Tolerance; end
    end
    
    % 初始磁阻率 (空表示默认线性)
    Nu_vec = [];
    
    fprintf('  [Transient] Starting Steps (dt=%.1e)...\n', dt);
    
    for n = 1:N_steps
        t_now = n * dt;
        V_source = V_func(t_now);
        
        % 松弛因子初始化
        alpha = 1.0; 
        err_prev = inf;
        
        for iter = 1:MaxIter
            % (A) 提取当前解
            A_k = x_curr(1:numEdges);
            I_k = x_curr(end);
            
            % (B) 更新物理参数 (Nu)
            CoilCurrent = Coil; CoilCurrent.I(:) = I_k * N_turns;
            
            % 使用上一各迭代步的 A 和 I 更新材料属性
            if isempty(Nu_vec) || iter > 1
                [Nu_new, max_B] = update_reluctivity_parallel(Model, A_k, CoilCurrent);
                % 松弛更新 Nu
                if isempty(Nu_vec), Nu_vec = Nu_new; end
                Nu_vec = 0.5 * Nu_new + 0.5 * Nu_vec; 
            end
            
            % (C) 组装线性化系统矩阵 K(nu)
            K_sec = assemble_magnetic_stiffness(Model, Nu_vec);
            
            % 磁化源 b_mag(nu)
            Q_vec = assemble_rhs_reduced(Model, CoilUnit, Nu_vec);
            
            % (D) 构建全局线性系统 LHS * x_new = RHS
            scale_K = mean(abs(nonzeros(K_sec)));
            eps_A = 1e-6 * scale_K;
            
            Block_AA = K_sec + M_sigma/dt + eps_A*M_reg;
            Block_AV = C_sigma;
            Block_VA = C_sigma' / dt;
            Block_VV = G_sigma;
            
            % 耦合项 (移到左边，符号变反)
            % F_col = Q*I + M*As/dt * I
            b_ind_unit = - (M_sigma * As_unit_vec) / dt;
            F_col = Q_vec + b_ind_unit;
            Block_AI = -F_col;
            
            Block_IA_A = W_vec_A' / dt;
            Block_IA_V = sparse(1, numActiveV);
            Block_II = CircuitProps.R + CircuitProps.L_leak/dt;
            
            GlobalK = [Block_AA,   Block_AV,   Block_AI;
                       Block_VA,   Block_VV,   sparse(numActiveV, 1);
                       Block_IA_A, Block_IA_V, sparse(Block_II)];
            
            % (E) 构建右端项 RHS
            RHS_A = (M_sigma * A_prev) / dt + (M_sigma * As_unit_vec * I_prev) / dt;
            RHS_V = (C_sigma' * A_prev) / dt;
            RHS_Circ = V_source + (W_vec_A' * A_prev)/dt + (CircuitProps.L_leak * I_prev)/dt;
            
            GlobalRHS = [RHS_A; RHS_V; RHS_Circ];
            
            % (F) 边界条件
            if isfield(Model.Runtime, 'FixedEdges'), fe=Model.Runtime.FixedEdges; else, fe=[]; end
            [K_solve, R_solve] = apply_dirichlet_bc(GlobalK, GlobalRHS, fe, zeros(size(fe)));
            
            % (G) 求解 x_new
            Model.Solver.Linear.Interface = 'MUMPS';
            Model.Solver.Linear.Symmetric = false;
            if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
            Model.Solver.Linear.MumpsICNTL.i8 = 77;
            
            try
                x_calc = linear_solve(K_solve, R_solve, Model);
                if any(isnan(x_calc)), error('MUMPS returned NaN'); end
            catch
                x_calc = K_solve \ R_solve; % Fallback
            end
            
            % (H) 自适应阻尼更新
            err = norm(x_calc - x_curr) / (norm(x_calc) + 1e-10);
            
            if iter > 1
                if err > err_prev
                    alpha = max(alpha * 0.33, 0.1); % 遇到困难，减速
                else
                    alpha = min(alpha * 1.1, 1.0); % 进展顺利，加速
                end
            end
            
            x_curr = (1 - alpha) * x_curr + alpha * x_calc;
            err_prev = err;
            
            if err < Tol
                break;
            end
        end
        
        if mod(n, 5) == 0 || n == 1
            fprintf('    Step %d: t=%.3fs, V=%.1fV, I=%.4fA (Max|B|=%.2fT, Iter=%d)\n', ...
                n, t_now, V_source, x_curr(end), max_B, iter);
        end
        
        x_prev = x_curr;
        I_prev = x_curr(end);
        A_prev = x_curr(1:numEdges);
        TimeHist.I(n) = I_prev;
    end
    
    Solution.A = x_curr(1:numEdges);
    Solution.I_history = TimeHist.I;
    Info.FinalTime = t_now;
end
