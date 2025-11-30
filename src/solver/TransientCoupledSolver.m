classdef TransientCoupledSolver < handle
    % TRANSIENTCOUPLEDSOLVER 场路耦合瞬态求解器 (v2.3 - Robust Fix)
    % 
    % 更新:
    %   1. 关闭 LinearSolver.ReuseAnalysis。
    
    properties
        Assembler
        LinearSolver
        
        MaxIter = 50
        Tolerance = 1e-5
        
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        WindingObj
        CircuitR
    end
    
    methods
        function obj = TransientCoupledSolver(assembler, winding, R)
            obj.Assembler = assembler;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            
            obj.LinearSolver = LinearSolver('Auto');
            % 耦合矩阵是非对称的
            obj.LinearSolver.Symmetric = false; 
            obj.LinearSolver.MumpsICNTL.i14 = 50; 
            
            % [Robustness] 关闭重用以防止 -53 错误
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [Solution, Info] = solve(obj, space, matLibData, timeParams, x0)
            fprintf('==============================================\n');
            fprintf('   Transient Coupled Solver (Robust)          \n');
            fprintf('==============================================\n');
            
            dt = timeParams.dt;
            N_steps = timeParams.N_steps;
            V_func = timeParams.V_func;
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            
            if nargin < 5 || isempty(x0)
                A_curr = zeros(numDofs, 1);
                I_curr = 0;
            else
                A_curr = x0(1:numDofs);
                I_curr = x0(end);
            end
            
            % 预计算耦合向量
            C_vec = obj.Assembler.assembleWinding(space, obj.WindingObj);
            is_bnd_dof = BoundaryCondition.findOuterBoundaryDofs(obj.Assembler.Mesh, obj.Assembler.DofHandler, space);
            is_fixed_sys = [is_bnd_dof; false];
            free_mask = ~is_fixed_sys;
            
            Solution.Time = zeros(N_steps, 1);
            Solution.I = zeros(N_steps, 1);
            Solution.V = zeros(N_steps, 1);
            Solution.A_history = cell(N_steps, 1);
            
            A_prev = A_curr;
            
            for n = 1:N_steps
                t_curr = n * dt;
                V_val = V_func(t_curr);
                
                % 初始残差 (Fast)
                [~, K_A_0] = obj.Assembler.assembleJacobian(space, A_curr, matLibData, false); 
                
                Res_A = K_A_0 - C_vec * I_curr;
                Psi_curr = C_vec' * A_curr;
                Psi_prev = C_vec' * A_prev;
                Res_I = (Psi_curr - Psi_prev) + dt * (obj.CircuitR * I_curr - V_val);
                
                % Auto-Scaling
                mag_scale = mean(abs(Res_A)); if mag_scale < 1e-5, mag_scale = 1.0; end
                cir_scale_base = abs(dt * obj.CircuitR); if cir_scale_base < 1e-10, cir_scale_base = 1e-10; end
                row_scale = mag_scale / cir_scale_base;
                row_scale = min(max(row_scale, 1e3), 1e9);
                
                Total_Res = [Res_A; Res_I * row_scale];
                current_res_norm = norm(Total_Res(free_mask));
                norm_base = current_res_norm; if norm_base < 1e-5, norm_base = 1.0; end
                
                fprintf('\n--- Step %d (t=%.4f) InitRes=%.2e ---\n', n, t_curr, current_res_norm);

                converged = false;
                
                for iter = 1:obj.MaxIter
                    % A. 组装 Jacobian
                    [J_mag, K_A] = obj.Assembler.assembleJacobian(space, A_curr, matLibData, true);
                    
                    % B. 残差
                    Res_A = K_A - C_vec * I_curr;
                    Psi_curr = C_vec' * A_curr;
                    Res_I = (Psi_curr - Psi_prev) + dt * (obj.CircuitR * I_curr - V_val);
                    
                    % C. 系统矩阵 (非对称)
                    % J = [ K     -C ]
                    %     [ C'*s   R'*s ]
                    J_sys = [J_mag, -C_vec; 
                             C_vec' * row_scale, sparse(dt * obj.CircuitR * row_scale)];
                         
                    Total_Res = [Res_A; Res_I * row_scale];
                    
                    % D. 求解
                    RHS = -Total_Res;
                    [J_solve, R_solve] = BoundaryCondition.applyDirichlet(J_sys, RHS, is_fixed_sys);
                    
                    dX_sys = obj.LinearSolver.solve(J_solve, R_solve);
                    
                    dA = dX_sys(1:end-1);
                    dI = dX_sys(end);
                    
                    % E. 线搜索
                    alpha = 1.0;
                    accepted = false;
                    
                    for bt = 0:obj.MaxBacktracks
                        A_try = A_curr + alpha * dA;
                        I_try = I_curr + alpha * dI;
                        
                        % Fast Residual
                        [~, K_A_try] = obj.Assembler.assembleJacobian(space, A_try, matLibData, false);
                        
                        R_A_try = K_A_try - C_vec * I_try;
                        Psi_try = C_vec' * A_try;
                        R_I_try = (Psi_try - Psi_prev) + dt * (obj.CircuitR * I_try - V_val);
                        
                        Total_Res_try = [R_A_try; R_I_try * row_scale];
                        new_res_norm = norm(Total_Res_try(free_mask));
                        
                        if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                            A_curr = A_try;
                            I_curr = I_try;
                            current_res_norm = new_res_norm;
                            accepted = true;
                            
                            marker = ''; if bt>0, marker=sprintf('[BT %d]', bt); end
                            fprintf('   Iter %d: Res=%.2e (Rel=%.2e) Step=%.4f %s\n', ...
                                iter, new_res_norm, new_res_norm/norm_base, alpha, marker);
                            break;
                        else
                            alpha = alpha * obj.BacktrackFactor;
                        end
                    end
                    
                    if ~accepted
                        warning('Line search failed.');
                        A_curr = A_curr + 1e-4 * dA; % 强制小步
                    end
                    
                    if current_res_norm / norm_base < obj.Tolerance
                        converged = true;
                        break;
                    end
                end
                
                if ~converged, warning('Step %d did not converge.', n); end
                
                A_prev = A_curr;
                Solution.A_history{n} = A_curr;
                Solution.I(n) = I_curr;
                Solution.V(n) = V_val;
                Solution.Time(n) = t_curr;
            end
            Info.Status = 'Done';
        end
    end
end