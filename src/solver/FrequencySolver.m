classdef FrequencySolver < handle
    % FREQUENCYSOLVER 频域时谐磁场求解器 (v1.3 - Robust Fix)
    
    properties
        Assembler
        LinearSolver
        Frequency % Hz
    end
    
    methods
        function obj = FrequencySolver(assembler, freq)
            obj.Assembler = assembler;
            obj.Frequency = freq;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 60; 
            
            % [Robustness]
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function A_phasor = solve(obj, space, nuMap, sigmaMap, sourceMap, fixedDofs)
            fprintf('==============================================\n');
            fprintf('   Frequency Domain Solver (%.1f Hz)          \n', obj.Frequency);
            fprintf('==============================================\n');
            
            omega = 2 * pi * obj.Frequency;
            
            % 1. 组装刚度矩阵 K
            K = obj.Assembler.assembleStiffness(space, nuMap);
            
            % 2. 组装加权质量矩阵 M_sigma
            if all(sigmaMap == 0)
                fprintf('   [Info] Conductivity is zero everywhere. Solving Magnetostatic AC.\n');
                M_sigma = sparse(size(K,1), size(K,2)); 
            else
                fprintf('   [Info] Assembling Eddy Current terms...\n');
                M_sigma = obj.Assembler.assembleMassWeighted(space, sigmaMap); 
            end
            
            % 3. 构建系统矩阵 S = K + j*w*M
            ref_val = mean(abs(diag(K)));
            eps_reg = ref_val * 1e-9;
            M_reg = obj.Assembler.assembleMass(space); 
            
            % S = K + j*w*M_sigma + eps*M_reg
            SystemMatrix = K + 1j * omega * M_sigma + eps_reg * M_reg;
            
            % 4. 组装载荷 F
            F = obj.Assembler.assembleSource(space, sourceMap);
            
            % 5. 施加边界条件
            [S_sys, F_sys] = BoundaryCondition.applyDirichlet(SystemMatrix, F, fixedDofs);
            
            % 6. 求解
            fprintf('   [Solver] Solving Complex Linear System...\n');
            tic;
            A_phasor = obj.LinearSolver.solve(S_sys, F_sys);
            t_solve = toc;
            
            fprintf('   -> Solved in %.4f s. |A| = %.4e\n', t_solve, norm(A_phasor));
            fprintf('==============================================\n');
        end
    end
end