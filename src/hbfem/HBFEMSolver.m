classdef HBFEMSolver < handle
    % HBFEMSOLVER 谐波平衡有限元求解器 (v3.1 - Optimized Line Search)
    % 
    % 更新:
    %   1. 在线搜索阶段显式关闭雅可比矩阵计算 (flag=false)，大幅提升速度。
    
    properties
        Assembler
        AFT
        LinearSolver
        
        MaxIter = 50
        Tolerance = 1e-4
        
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        MatrixK_Add 
    end
    
    methods
        function obj = HBFEMSolver(assembler, aft)
            obj.Assembler = assembler;
            obj.AFT = aft;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
        end
        
        function [X_harmonics, info] = solve(obj, space, matLibData, sourceMaps, fixedDofs, x0)
            fprintf('==============================================\n');
            fprintf('   HBFEM Solver (Optimized Line Search)       \n');
            fprintf('==============================================\n');
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            numHarm = obj.AFT.NumHarmonics;
            totalDofs = numDofs * numHarm;
            
            if nargin < 6 || isempty(x0), X_mat = zeros(numDofs, numHarm); else, X_mat = x0; end
            
            % 1. 组装全局载荷 F
            F_global = zeros(totalDofs, 1);
            for k = 1:numHarm
                if ~isempty(sourceMaps{k})
                    F_k = obj.Assembler.assembleSource(space, sourceMaps{k});
                    idx_start = (k-1)*numDofs + 1;
                    idx_end = k*numDofs;
                    F_global(idx_start:idx_end) = F_k;
                end
            end
            
            % 2. 正则化矩阵
            K_add_global = sparse(totalDofs, totalDofs);
            has_reg = ~isempty(obj.MatrixK_Add);
            
            if has_reg
                fprintf('   [Info] Applying Regularization to all harmonics.\n');
                M_small = obj.MatrixK_Add;
                [mi, mj, mv] = find(M_small);
                n_nz = length(mv);
                I_reg = zeros(n_nz * numHarm, 1);
                J_reg = zeros(n_nz * numHarm, 1);
                V_reg = zeros(n_nz * numHarm, 1);
                
                for h = 1:numHarm
                    offset = (h-1) * numDofs;
                    idx_range = (h-1)*n_nz + 1 : h*n_nz;
                    I_reg(idx_range) = mi + offset;
                    J_reg(idx_range) = mj + offset;
                    V_reg(idx_range) = mv;
                end
                K_add_global = sparse(I_reg, J_reg, V_reg, totalDofs, totalDofs);
            end
            
            Fixed_global = repmat(fixedDofs(:), numHarm, 1);
            
            % 3. 初始残差 (Flag = false, 只算 R)
            [~, R_mat_0] = obj.Assembler.assembleHBFEM(space, X_mat, obj.AFT, matLibData, false);
            R_vec_0 = R_mat_0(:);
            Total_Res_0 = R_vec_0 - F_global;
            if has_reg
                Total_Res_0 = Total_Res_0 + K_add_global * X_mat(:);
            end
            
            free_mask = ~Fixed_global;
            current_res_norm = norm(Total_Res_0(free_mask));
            norm_base = norm(F_global); if norm_base < 1e-10, norm_base = 1.0; end
            
            fprintf('   Iter  0: Res=%.4e (Rel=%.4e)\n', current_res_norm, current_res_norm/norm_base);
            
            res_history = current_res_norm;
            converged = false;
            
            for iter = 1:obj.MaxIter
                % A. 组装物理 Jacobian (Flag = true, 需要 J)
                [J_tri, R_mat] = obj.Assembler.assembleHBFEM(space, X_mat, obj.AFT, matLibData, true);
                
                R_vec = R_mat(:);
                J_phys = sparse(J_tri.I, J_tri.J, J_tri.V, totalDofs, totalDofs);
                
                % B. 添加正则化
                if has_reg
                    X_vec = X_mat(:);
                    Total_Res = R_vec + (K_add_global * X_vec) - F_global;
                    J_total = J_phys + K_add_global;
                else
                    Total_Res = R_vec - F_global;
                    J_total = J_phys;
                end
                
                % C. 求解
                RHS = -Total_Res;
                [J_sys, RHS_sys] = BoundaryCondition.applyDirichlet(J_total, RHS, Fixed_global);
                
                try
                    dX_vec = obj.LinearSolver.solve(J_sys, RHS_sys);
                catch
                    warning('Linear solver failed. Retrying with identity perturbation.');
                    J_sys = J_sys + 1e-3 * speye(size(J_sys));
                    dX_vec = obj.LinearSolver.solve(J_sys, RHS_sys);
                end
                
                dX_mat = reshape(dX_vec, numDofs, numHarm);
                
                % D. 线搜索 (优化点: Flag = false)
                alpha = 1.0;
                accepted = false;
                
                for bt = 0:obj.MaxBacktracks
                    X_trial = X_mat + alpha * dX_mat;
                    
                    % [Optimization] 仅计算残差，跳过 Jacobian
                    [~, R_mat_trial] = obj.Assembler.assembleHBFEM(space, X_trial, obj.AFT, matLibData, false);
                    R_vec_trial = R_mat_trial(:);
                    
                    if has_reg
                        X_vec_trial = X_trial(:);
                        Res_trial = R_vec_trial + (K_add_global * X_vec_trial) - F_global;
                    else
                        Res_trial = R_vec_trial - F_global;
                    end
                    
                    if any(isnan(Res_trial)) || any(isinf(Res_trial))
                        new_res_norm = inf;
                    else
                        new_res_norm = norm(Res_trial(free_mask));
                    end
                    
                    if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                        X_mat = X_trial;
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
                    warning('Line search stuck. Forcing small step.');
                    X_mat = X_mat + 1e-3 * dX_mat; 
                    [~, R_mat_next] = obj.Assembler.assembleHBFEM(space, X_mat, obj.AFT, matLibData, false);
                    current_res_norm = norm(R_mat_next(:) - F_global); 
                end
                
                res_history(end+1) = current_res_norm; %#ok<AGROW>
                
                if current_res_norm / norm_base < obj.Tolerance
                    fprintf('   -> Converged!\n');
                    converged = true;
                    break;
                end
            end
            
            X_harmonics = X_mat;
            info.Iterations = iter;
            info.Residuals = res_history;
            info.Converged = converged;
            fprintf('==============================================\n');
        end
    end
end