classdef NonlinearSolver < handle
    % NONLINEARSOLVER 牛顿-拉夫逊非线性求解器 (v4.6 - Robust Fix)
    % 
    % 更新:
    %   1. [Fix] 关闭 ReuseAnalysis。MATLAB sparse() 会自动丢弃零元，
    %      导致 Jacobian 结构动态变化，MUMPS 重用分析会报错 -53。
    
    properties
        Assembler
        Config
        LinearSolver
        
        MaxIter = 100
        Tolerance = 1e-5
        
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        MatrixK_Add
        VectorF_Add
    end
    
    methods
        function obj = NonlinearSolver(assembler)
            obj.Assembler = assembler;
            obj.Config = assembler.Config;
            obj.LinearSolver = LinearSolver('Auto');
            
            % [Robustness] 必须关闭分析复用，因为 sparse() 不保证结构恒定
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [x, info] = solve(obj, space, matLibData, sourceMap, fixedDofs, x0)
            fprintf('==============================================\n');
            fprintf('   Nonlinear Solver (Fast Newton)             \n');
            fprintf('==============================================\n');
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            
            if nargin < 6 || isempty(x0), x = zeros(numDofs, 1); else, x = x0; end
            
            if isempty(sourceMap)
                F_ext = sparse(numDofs, 1);
            else
                F_ext = obj.Assembler.assembleSource(space, sourceMap);
            end
            
            has_add_K = ~isempty(obj.MatrixK_Add);
            
            % Initial Residual (Only need R)
            [~, R_mag_0] = obj.Assembler.assembleJacobian(space, x, matLibData, false);
            
            R_0 = R_mag_0 - F_ext;
            if has_add_K, R_0 = R_0 + (obj.MatrixK_Add * x); end
            if ~isempty(obj.VectorF_Add), R_0 = R_0 - obj.VectorF_Add; end
            
            free_mask = ~fixedDofs;
            current_res_norm = norm(R_0(free_mask));
            norm_base = norm(F_ext); 
            if has_add_K && ~isempty(obj.VectorF_Add), norm_base = max(norm_base, norm(obj.VectorF_Add)); end
            if norm_base < 1e-10, norm_base = 1.0; end
            
            fprintf('   Iter  0: Res=%.4e (Rel=%.4e)\n', current_res_norm, current_res_norm/norm_base);
            
            res_history = current_res_norm;
            converged = false;
            
            for iter = 1:obj.MaxIter
                % A. 组装 Jacobian (Need J)
                [J_mag, R_mag] = obj.Assembler.assembleJacobian(space, x, matLibData, true);
                
                % B. 系统组装
                if has_add_K
                    J_total = J_mag + obj.MatrixK_Add;
                    R_total = R_mag + (obj.MatrixK_Add * x) - F_ext;
                else
                    J_total = J_mag;
                    R_total = R_mag - F_ext;
                end
                if ~isempty(obj.VectorF_Add), R_total = R_total - obj.VectorF_Add; end
                
                % C. 求解
                RHS = -R_total;
                [J_sys, RHS_sys] = BoundaryCondition.applyDirichlet(J_total, RHS, fixedDofs);
                dx = obj.LinearSolver.solve(J_sys, RHS_sys);
                
                % D. 线搜索 (Fast: No Jacobian)
                alpha = 1.0;
                accepted = false;
                
                for bt = 0:obj.MaxBacktracks
                    x_trial = x + alpha * dx;
                    
                    % [Fast Calc] 仅计算残差
                    [~, R_mag_trial] = obj.Assembler.assembleJacobian(space, x_trial, matLibData, false);
                    
                    if has_add_K
                        R_trial = R_mag_trial + (obj.MatrixK_Add * x_trial) - F_ext;
                    else
                        R_trial = R_mag_trial - F_ext;
                    end
                    if ~isempty(obj.VectorF_Add), R_trial = R_trial - obj.VectorF_Add; end
                    
                    new_res_norm = norm(R_trial(free_mask));
                    
                    if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                        x = x_trial;
                        current_res_norm = new_res_norm;
                        accepted = true;
                        
                        marker = ''; if bt>0, marker=sprintf('[BT %d]', bt); end
                        fprintf('   Iter %2d: Res=%.4e (Rel=%.4e) Step=%.4f %s\n', ...
                            iter, current_res_norm, current_res_norm/norm_base, alpha, marker);
                        break;
                    else
                        alpha = alpha * obj.BacktrackFactor;
                    end
                end
                
                if ~accepted
                    warning('Line search failed.');
                    break;
                end
                
                res_history(end+1) = current_res_norm; %#ok<AGROW>
                
                if current_res_norm / norm_base < obj.Tolerance
                    fprintf('   -> Converged!\n');
                    converged = true;
                    break;
                end
            end
            
            if ~converged
                warning('Nonlinear solver reached max iterations without convergence.');
            end
            
            info.Iterations = iter;
            info.Residuals = res_history;
            info.Converged = converged;
            fprintf('==============================================\n');
        end
    end
end