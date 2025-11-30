classdef FrequencyCoupledSolver < handle
    % FREQUENCYCOUPLEDSOLVER 频域场路耦合求解器 (v1.2 - Robust Fix)
    
    properties
        Assembler
        LinearSolver
        Frequency % Hz
        WindingObj
        CircuitR
    end
    
    methods
        function obj = FrequencyCoupledSolver(assembler, winding, R, freq)
            obj.Assembler = assembler;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            obj.Frequency = freq;
            
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 60; 
            obj.LinearSolver.Symmetric = false; 
            
            % [Robustness]
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [Solution, Info] = solve(obj, space, nuMap, sigmaMap, V_phasor, fixedDofs)
            fprintf('==============================================\n');
            fprintf('   Frequency Coupled Solver (%.1f Hz)         \n', obj.Frequency);
            fprintf('==============================================\n');
            
            omega = 2 * pi * obj.Frequency;
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            
            % 1. 组装磁场矩阵
            fprintf('[Coupled] Assembling Field Matrices...\n');
            K = obj.Assembler.assembleStiffness(space, nuMap);
            
            if all(sigmaMap == 0)
                 M_sigma = sparse(numDofs, numDofs);
            else
                 M_sigma = obj.Assembler.assembleMassWeighted(space, sigmaMap);
            end
            
            % 正则化
            ref_val = mean(abs(diag(K)));
            eps_reg = ref_val * 1e-9;
            M_reg = obj.Assembler.assembleMass(space);
            
            % Z_mag = K + j*w*M_sigma + eps*M_reg
            Z_mag = K + 1j * omega * M_sigma + eps_reg * M_reg;
            
            % 2. 组装耦合向量
            C_vec = obj.Assembler.assembleWinding(space, obj.WindingObj);
            
            % 3. 组装系统矩阵
            sys_size = numDofs + 1;
            top_right = -C_vec;
            bot_left = 1j * omega * C_vec';
            bot_right = sparse(obj.CircuitR);
            
            SystemMatrix = [Z_mag, top_right; bot_left, bot_right];
            
            % 4. 组装载荷向量
            RHS = zeros(sys_size, 1);
            RHS(end) = V_phasor; 
            
            % 5. 处理边界条件
            is_fixed_sys = [fixedDofs; false];
            [S_sys, F_sys] = BoundaryCondition.applyDirichlet(SystemMatrix, RHS, is_fixed_sys);
            
            % 6. 求解
            fprintf('[Solver] Solving Complex Coupled System...\n');
            tic;
            X = obj.LinearSolver.solve(S_sys, F_sys);
            t_solve = toc;
            
            % 7. 结果提取
            A_phasor = X(1:numDofs);
            I_phasor = X(end);
            
            fprintf('   -> Solved in %.4f s.\n', t_solve);
            fprintf('   -> Current I = %.4f + %.4fi A (Mag: %.4f)\n', ...
                real(I_phasor), imag(I_phasor), abs(I_phasor));
            
            Solution.A = A_phasor;
            Solution.I = I_phasor;
            Info.SystemSize = sys_size;
        end
    end
end