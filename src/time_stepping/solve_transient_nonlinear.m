function [Solution, Info] = solve_transient_nonlinear(Model, Coil, TimeParams)
% SOLVE_TRANSIENT_NONLINEAR 非线性瞬态磁场求解器 (电流驱动 - Picard版)
%
% 算法: 定点迭代 (Fixed-Point) + 割线刚度矩阵 (Secant Stiffness)
% 优势: 
%   1. 无需计算雅可比矩阵，避免了微分磁导率可能出现的负值/奇异问题。
%   2. 刚度矩阵 K(nu) 恒正定，MUMPS 求解极其稳定。

    fprintf('==============================================\n');
    fprintf('   Nonlinear Transient Solver (Current Driven)\n');
    fprintf('   Method: Picard Iteration (Secant Stiffness)\n');
    fprintf('==============================================\n');

    dt = TimeParams.dt;
    N_steps = TimeParams.N_steps;
    I_func = TimeParams.I_func;
    
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
    
    if isfield(Coil, 'Turns'), N_turns = Coil.Turns; else, N_turns = 1.0; end
    
    % 感应源单位向量
    CoilUnit = Coil; CoilUnit.I(:) = 1.0 * N_turns;
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    
    % 初始化
    x_curr = zeros(numFieldDoFs, 1); 
    x_prev = x_curr;
    I_prev = 0;
    if isfield(TimeParams, 'I_init'), I_prev = TimeParams.I_init; end
    
    TimeHist.I = zeros(N_steps, 1);
    
    MaxIter = 30; Tol = 1e-4;
    if isfield(Model, 'Solver') && isfield(Model.Solver, 'Nonlinear')
        if isfield(Model.Solver.Nonlinear, 'MaxIter'), MaxIter=Model.Solver.Nonlinear.MaxIter; end
        if isfield(Model.Solver.Nonlinear, 'Tolerance'), Tol=Model.Solver.Nonlinear.Tolerance; end
    end
    
    % 初始磁阻率
    Nu_vec = [];
    
    fprintf('  [Transient] Starting Steps (dt=%.1e)...\n', dt);
    
    for n = 1:N_steps
        t_now = n * dt;
        I_now = I_func(t_now);
        dI_dt = (I_now - I_prev) / dt;
        
        alpha = 1.0; err_prev = inf;
        
        for iter = 1:MaxIter
            A_k = x_curr(1:numEdges);
            
            % (A) 更新物理参数 (Nu)
            CoilCurrent = Coil; CoilCurrent.I(:) = I_now * N_turns;
            [Nu_new, max_B] = update_reluctivity_parallel(Model, A_k, CoilCurrent);
            
            if isempty(Nu_vec), Nu_vec = Nu_new; end
            Nu_vec = 0.5 * Nu_new + 0.5 * Nu_vec; 
            
            % (B) 组装割线刚度 K_sec
            K_sec = assemble_magnetic_stiffness(Model, Nu_vec);
            b_mag = assemble_rhs_reduced(Model, CoilCurrent, Nu_vec);
            
            % (C) 构建线性系统
            scale_K = mean(abs(nonzeros(K_sec)));
            eps_A = 1e-6 * scale_K;
            
            K_eff = K_sec + M_sigma/dt + eps_A*M_reg;
            Ct = C_sigma';
            G_eff = G_sigma * dt; % 对称形式
            
            SysK = [K_eff, C_sigma; Ct, G_eff];
            
            % (D) 构建 RHS
            b_ind = - M_sigma * (As_unit_vec * dI_dt);
            b_total = b_mag + b_ind;
            
            A_prev_vec = x_prev(1:numEdges);
            RHS_A = b_total + (M_sigma * A_prev_vec) / dt;
            % RHS_V = (C_sigma' * A_prev_vec) / dt; % 对称形式: C'/dt * A_prev ? 
            % 检查方程: C' * (A - Aold)/dt + G * V = 0
            % 乘 dt -> C' * A + G*dt * V = C' * Aold
            % 正确，这里 RHS_V = C' * Aold
            RHS_V = C_sigma' * A_prev_vec;
            
            SysRHS = [RHS_A; RHS_V];
            
            % (E) 边界条件与求解
            if isfield(Model.Runtime, 'FixedEdges'), fe=Model.Runtime.FixedEdges; else, fe=[]; end
            [K_solve, R_solve] = apply_dirichlet_bc(SysK, SysRHS, fe, zeros(size(fe)));
            
            Model.Solver.Linear.Interface = 'MUMPS';
            Model.Solver.Linear.Symmetric = true;
            if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
            Model.Solver.Linear.MumpsICNTL.i8 = 77;
            
            try
                x_calc = linear_solve(K_solve, R_solve, Model);
                if any(isnan(x_calc)), error('NaN'); end
            catch
                x_calc = K_solve \ R_solve;
            end
            
            % (F) 阻尼更新
            err = norm(x_calc - x_curr) / (norm(x_calc) + 1e-10);
            
            if iter > 1
                if err > err_prev
                    alpha = max(alpha * 0.33, 0.1);
                else
                    alpha = min(alpha * 1.1, 1.0);
                end
            end
            x_curr = (1 - alpha) * x_curr + alpha * x_calc;
            err_prev = err;
            
            if err < Tol, break; end
        end
        
        if mod(n, 5) == 0 || n == 1
            fprintf('    Step %2d: t=%.3fs, I=%.2fA (Max|B|=%.2fT, Iter=%d)\n', ...
                n, t_now, I_now, max_B, iter);
        end
        
        x_prev = x_curr;
        I_prev = I_now;
        TimeHist.I(n) = I_prev;
    end
    
    Solution.A = x_curr(1:numEdges);
    Solution.I_history = TimeHist.I;
    Info.FinalTime = t_now;
end
