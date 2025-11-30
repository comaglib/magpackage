classdef TransientNonlinearSolver < handle
    % TRANSIENTNONLINEARSOLVER 非线性瞬态求解器 (Class Version v1.0)
    %
    % 功能:
    %   求解方程: Curl(Nu * Curl A) + Sigma * dA/dt = J
    %   算法: Implicit Euler (时间离散) + Newton-Raphson (非线性迭代)
    %
    % 接口变更:
    %   旧(函数): [Sol, Info] = TransientNonlinearSolver(Assembler, TimeParams, MatLibData, Space)
    %   新(类):   solver = TransientNonlinearSolver(Assembler);
    %             [Sol, Info] = solver.solve(space, matLibData, timeParams, x0);
    
    properties
        Assembler
        NonlinearSolver % 内部持有的非线性求解器实例 (Newton-Raphson)
    end
    
    methods
        function obj = TransientNonlinearSolver(assembler)
            % 构造函数
            obj.Assembler = assembler;
            
            % 初始化内部非线性求解器
            obj.NonlinearSolver = NonlinearSolver(assembler);
            
            % 默认配置
            obj.NonlinearSolver.Tolerance = 1e-4;
            
            % 配置底层线性求解器 (MUMPS)
            obj.NonlinearSolver.LinearSolver.Method = 'Auto';
            obj.NonlinearSolver.LinearSolver.MumpsICNTL.i14 = 40;
            
            % [Robustness] 必须关闭符号分析复用，防止稀疏结构变化导致 Crash
            obj.NonlinearSolver.LinearSolver.ReuseAnalysis = false;
        end
        
        function [Solution, Info] = solve(obj, space, matLibData, timeParams, x0)
            % SOLVE 执行瞬态求解
            % 输入:
            %   space      - FunctionSpace 对象
            %   matLibData - 材料库数据
            %   timeParams - 结构体，包含 dt, N_steps, Sigma, SourceFunc
            %   x0         - (可选) 初始条件向量 A_0
            
            fprintf('==============================================\n');
            fprintf('   Transient Nonlinear Solver (Class v1.0)    \n');
            fprintf('   Method: Implicit Euler + Newton-Raphson    \n');
            fprintf('==============================================\n');
            
            dt = timeParams.dt;
            N_steps = timeParams.N_steps;
            
            % 1. 准备质量矩阵 (Mass Matrices)
            fprintf('[Transient] Assembling Mass Matrices...\n');
            M_total = obj.Assembler.assembleMass(space);
            
            Sigma_val = 0;
            if isfield(timeParams, 'Sigma'), Sigma_val = timeParams.Sigma; end
            
            % M_sigma: 用于物理时间导数项 (Sigma * dA/dt)
            M_sigma = M_total * Sigma_val;
            
            % M_reg: 用于正则化 (防止空气域奇异)
            reg_epsilon = 1e-9;
            M_reg = M_total * reg_epsilon;
            
            % 2. 初始条件
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            if nargin < 5 || isempty(x0)
                A_curr = zeros(numDofs, 1);
            else
                A_curr = x0;
            end
            
            % 识别边界条件
            is_bnd = BoundaryCondition.findOuterBoundaryDofs(obj.Assembler.Mesh, obj.Assembler.DofHandler, space);
            
            % 3. 时间步进
            Solution.A_history = cell(N_steps, 1);
            Solution.Time = zeros(N_steps, 1);
            
            % 预计算时间步进的附加刚度矩阵 (Backward Euler)
            % 方程: (K(A) + M_sigma/dt + M_reg) * A_n = J + (M_sigma/dt)*A_{n-1}
            % K_add = M_sigma/dt + M_reg
            K_add = (M_sigma / dt) + M_reg;
            
            for n = 1:N_steps
                t_curr = n * dt;
                fprintf('\n--- Time Step %d/%d (t=%.4f s) ---\n', n, N_steps, t_curr);
                
                % 更新源项
                if isa(timeParams.SourceFunc, 'function_handle')
                    SourceMap = timeParams.SourceFunc(t_curr);
                else
                    SourceMap = timeParams.SourceFunc;
                end
                
                % 计算附加力向量 (来自上一时刻)
                % F_add = (M_sigma / dt) * A_{n-1}
                F_add = (M_sigma * A_curr) / dt;
                
                % 将附加项注入到 NonlinearSolver
                obj.NonlinearSolver.MatrixK_Add = K_add;
                obj.NonlinearSolver.VectorF_Add = F_add;
                
                % 调用 Newton 求解当前步
                % solve(space, matLibData, sourceMap, fixedDofs, x0)
                [A_next, stepInfo] = obj.NonlinearSolver.solve(space, matLibData, SourceMap, is_bnd, A_curr);
                
                if ~stepInfo.Converged
                    warning('Time step %d failed to converge. Result may be inaccurate.', n);
                end
                
                % 更新状态
                A_curr = A_next;
                
                % 记录结果
                Solution.A_history{n} = A_curr;
                Solution.Time(n) = t_curr;
            end
            
            Info.Status = 'Completed';
            Info.FinalLinearSolver = obj.NonlinearSolver.LinearSolver.Method;
        end
    end
end