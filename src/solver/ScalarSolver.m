classdef ScalarSolver < handle
    % SCALARSOLVER 标量场求解器 (v1.0)
    % 
    % 功能:
    %   求解 div(sigma * grad V) = 0
    %   应用: 直流电流传导、静电场、热传导
    
    properties
        Assembler
        LinearSolver
    end
    
    methods
        function obj = ScalarSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 40;
        end
        
        function [V_sol, Info] = solve(obj, space, sigmaMap, bcNodes, bcValues)
            % SOLVE 执行标量场求解
            % 输入:
            %   sigmaMap  - 区域电导率映射 [RegionID -> Sigma]
            %   bcNodes   - Dirichlet 边界节点索引 [N_bc x 1]
            %   bcValues  - Dirichlet 边界值 [N_bc x 1]
            
            fprintf('==============================================\n');
            fprintf('   Scalar Field Solver (Conduction)           \n');
            fprintf('==============================================\n');
            
            % 1. 组装刚度矩阵 S
            S = obj.Assembler.assembleScalarStiffness(space, sigmaMap);
            numDofs = size(S, 1);
            
            % 2. 载荷向量 F (默认为 0，电流守恒)
            % 如果有体电流源，可以在此添加
            F = sparse(numDofs, 1);
            
            % 3. 施加边界条件
            % 构建逻辑掩码
            is_fixed = false(numDofs, 1);
            is_fixed(bcNodes) = true;
            
            % 这里我们需要一个处理非零 Dirichlet 值的函数
            % 之前的 BoundaryCondition.applyDirichlet 主要是为齐次或简单情况设计的
            % 对于非零 BC (V=V0)，通常做法:
            % F = F - S(:, fixed) * V_fixed
            % S(fixed, :) = 0; S(:, fixed) = 0; S(fixed, fixed) = I
            % F(fixed) = V_fixed
            
            % 为了复用，我们手动处理非零 BC 的 RHS 修正
            V_fixed_full = zeros(numDofs, 1);
            V_fixed_full(bcNodes) = bcValues;
            
            % 移项: F_new = F - S * V_known
            F = F - S * V_fixed_full;
            
            % 调用通用 BC 处理 (置1法)
            [S_sys, F_sys] = BoundaryCondition.applyDirichlet(S, F, is_fixed);
            
            % 修正 F_sys 的固定位 (applyDirichlet 默认设为 0，我们需要设为 bcValues - 实际上是 0 因为是增量?)
            % 等等，ApplyDirichlet 的逻辑是: F_sys(fixed) = 0
            % 如果我们解的是增量 dV，这是对的。
            % 但我们要解全量 V。
            % 标准置1法: S_ii=1, F_i = V_bc
            
            F_sys(bcNodes) = bcValues;
            
            % 4. 求解
            tic;
            V_sol = obj.LinearSolver.solve(S_sys, F_sys);
            t_solve = toc;
            
            fprintf('   -> Solved in %.4f s.\n', t_solve);
            fprintf('==============================================\n');
            
            Info.Time = t_solve;
        end
    end
end