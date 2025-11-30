classdef ScalarSolver < handle
    % SCALARSOLVER 标量场求解器 (v1.3 - Symmetric Optimized)
    %
    % 改进:
    %   1. [Optimization] 启用线性求解器对称正定模式 (MumpsSymmetry = 1)。
    
    properties
        Assembler
        LinearSolver
    end
    
    methods
        function obj = ScalarSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 40;
            obj.LinearSolver.ReuseAnalysis = false;
            
            % [New] 启用对称正定模式
            % SYM = 1: Symmetric Positive Definite (SPD)
            % 标量拉普拉斯矩阵是实对称且正定的
            obj.LinearSolver.MumpsSymmetry = 1;
        end
        
        function [V_sol, Info] = solve(obj, space, sigmaMap, bcNodes, bcValues)
            fprintf('==============================================\n');
            fprintf('   Scalar Field Solver (Symmetric Mode)       \n');
            fprintf('==============================================\n');
            
            % 1. 组装刚度矩阵 S
            S = obj.Assembler.assembleScalarLaplacian(space, sigmaMap);
            numDofs = size(S, 1);
            
            % 2. 载荷向量
            F = sparse(numDofs, 1);
            
            % 3. 施加边界条件
            is_fixed = false(numDofs, 1);
            is_fixed(bcNodes) = true;
            
            V_fixed_full = zeros(numDofs, 1);
            V_fixed_full(bcNodes) = bcValues;
            
            % 提升边界条件到右端项
            F = F - S * V_fixed_full;
            
            [S_sys, F_sys] = BoundaryCondition.applyDirichlet(S, F, is_fixed);
            F_sys(bcNodes) = bcValues;
            
            % 4. 求解
            % LinearSolver 将利用 MumpsSymmetry = 1 进行加速
            tic;
            V_sol = obj.LinearSolver.solve(S_sys, F_sys);
            t_solve = toc;
            
            Info.SolveTime = t_solve;
            fprintf('   Solved in %.4f s.\n', t_solve);
        end
    end
end