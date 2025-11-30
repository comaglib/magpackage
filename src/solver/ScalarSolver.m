classdef ScalarSolver < handle
    % SCALARSOLVER 标量场求解器 (v1.2 - Robust Fix)
    
    properties
        Assembler
        LinearSolver
    end
    
    methods
        function obj = ScalarSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 40;
            
            % [Robustness]
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [V_sol, Info] = solve(obj, space, sigmaMap, bcNodes, bcValues)
            fprintf('==============================================\n');
            fprintf('   Scalar Field Solver (Conduction)           \n');
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
            
            F = F - S * V_fixed_full;
            
            [S_sys, F_sys] = BoundaryCondition.applyDirichlet(S, F, is_fixed);
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