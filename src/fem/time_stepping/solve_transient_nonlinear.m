function [Solution, Info] = solve_transient_nonlinear(Assembler, TimeParams, MatLibData, Space)
% SOLVE_TRANSIENT_NONLINEAR 非线性瞬态求解器 (v3.2 - Robust)
% 
% 求解方程: Curl(Nu * Curl A) + Sigma * dA/dt = J
% 算法: Implicit Euler + Newton-Raphson
%
% 更新:
%   1. 自动处理空气域的正则化 (防止 Sigma=0 区域奇异)。
%   2. 适配 v4.2 NonlinearSolver 接口。

    fprintf('==============================================\n');
    fprintf('   Transient Nonlinear Solver (v3.2)          \n');
    fprintf('   Method: Implicit Euler + Newton-Raphson    \n');
    fprintf('==============================================\n');

    dt = TimeParams.dt;
    N_steps = TimeParams.N_steps;
    
    % 1. 初始化求解器
    Solver = NonlinearSolver(Assembler);
    Solver.Tolerance = 1e-4; 
    
    % 配置线性求解器 (MUMPS)
    Solver.LinearSolver.Method = 'Auto';
    Solver.LinearSolver.MumpsICNTL.i14 = 40; 
    
    % 2. 准备质量矩阵
    % M_total: 用于正则化 (几何质量矩阵)
    % M_sigma: 用于物理时间项 (电导率加权)
    fprintf('[Transient] Assembling Mass Matrices...\n');
    M_total = Assembler.assembleMass(Space);
    
    Sigma_val = 0; 
    if isfield(TimeParams, 'Sigma'), Sigma_val = TimeParams.Sigma; end
    
    % 假设全域导电 (简化演示)。实际应传入 SigmaMap。
    M_sigma = M_total * Sigma_val; 
    
    % [正则化策略]
    % 为了防止空气域 (Sigma=0) 在瞬态步进中奇异，我们需要添加微小的 epsilon * M_total
    % 系数选取: 相对于 K_stiffness 的量级
    % 简单的估算: eps ~ 1e-8 / dt (经验值) 或者基于 K 的对角线
    % 这里我们给一个保守的固定正则化，或者让 Assembler 算 K_lin 来估算
    reg_epsilon = 1e-9; 
    M_reg = M_total * reg_epsilon;
    
    % 3. 初始条件
    numDofs = Assembler.DofHandler.NumGlobalDofs;
    A_curr = zeros(numDofs, 1);
    
    % 边界条件
    is_bnd = BoundaryCondition.findOuterBoundaryDofs(Assembler.Mesh, Assembler.DofHandler, Space);
    
    % 4. 时间步进
    Solution.A_history = cell(N_steps, 1);
    Solution.Time = zeros(N_steps, 1);
    
    % 预计算时间步进的附加刚度矩阵
    % K_add = M_sigma / dt + M_reg
    K_add = (M_sigma / dt) + M_reg;
    
    for n = 1:N_steps
        t_curr = n * dt;
        fprintf('\n--- Time Step %d/%d (t=%.4f s) ---\n', n, N_steps, t_curr);
        
        % 源项更新
        if isa(TimeParams.SourceFunc, 'function_handle')
            SourceMap = TimeParams.SourceFunc(t_curr);
        else
            SourceMap = TimeParams.SourceFunc; 
        end
        
        % 附加力向量 (Backward Euler)
        % F_add = (M_sigma / dt) * A_n
        F_add = (M_sigma * A_curr) / dt;
        
        % 配置求解器
        Solver.MatrixK_Add = K_add;
        Solver.VectorF_Add = F_add;
        
        % 求解 (使用上一时刻解作为初值 A_curr)
        [A_next, info] = Solver.solve(Space, MatLibData, SourceMap, is_bnd, A_curr);
        
        if ~info.Converged
            warning('Time step %d failed to converge. Result may be inaccurate.', n);
        end
        
        % 更新状态
        A_curr = A_next;
        
        % 存储
        Solution.A_history{n} = A_curr;
        Solution.Time(n) = t_curr;
    end
    
    Info.Status = 'Completed';
    Info.FinalLinearSolver = Solver.LinearSolver.Method;
end