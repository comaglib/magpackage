function [Solution, Info] = TransientNonlinearSolver(Assembler, TimeParams, MatLibData, Space)
% SOLVE_TRANSIENT_NONLINEAR 非线性瞬态求解器 (v3.4 - Robust Fix)
% 
% 求解方程: Curl(Nu * Curl A) + Sigma * dA/dt = J
% 算法: Implicit Euler + Newton-Raphson
%
% 更新:
%   1. 关闭 MUMPS 符号分析复用 (ReuseAnalysis = false)。

    fprintf('==============================================\n');
    fprintf('   Transient Nonlinear Solver (v3.4)          \n');
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
    
    % [Robustness] 必须关闭
    Solver.LinearSolver.ReuseAnalysis = false;
    
    % 2. 准备质量矩阵
    fprintf('[Transient] Assembling Mass Matrices...\n');
    M_total = Assembler.assembleMass(Space);
    
    Sigma_val = 0; 
    if isfield(TimeParams, 'Sigma'), Sigma_val = TimeParams.Sigma; end
    
    M_sigma = M_total * Sigma_val; 
    
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
    
    K_add = (M_sigma / dt) + M_reg;
    
    for n = 1:N_steps
        t_curr = n * dt;
        fprintf('\n--- Time Step %d/%d (t=%.4f s) ---\n', n, N_steps, t_curr);
        
        if isa(TimeParams.SourceFunc, 'function_handle')
            SourceMap = TimeParams.SourceFunc(t_curr);
        else
            SourceMap = TimeParams.SourceFunc; 
        end
        
        F_add = (M_sigma * A_curr) / dt;
        
        Solver.MatrixK_Add = K_add;
        Solver.VectorF_Add = F_add;
        
        [A_next, info] = Solver.solve(Space, MatLibData, SourceMap, is_bnd, A_curr);
        
        if ~info.Converged
            warning('Time step %d failed to converge. Result may be inaccurate.', n);
        end
        
        A_curr = A_next;
        
        Solution.A_history{n} = A_curr;
        Solution.Time(n) = t_curr;
    end
    
    Info.Status = 'Completed';
    Info.FinalLinearSolver = Solver.LinearSolver.Method;
end