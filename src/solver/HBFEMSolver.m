classdef HBFEMSolver < handle
    % HBFEMSOLVER 谐波平衡有限元求解器 (v4.1 - Sparse Print Fix)
    % 
    % 修复日志:
    %   1. [Print Error] 修复 fprintf 不支持稀疏标量输入的问题 (使用 full() 转换)。
    %   2. [Dimension] 保持 reshape 逻辑以匹配 Assembler。
    
    properties
        Assembler
        AFT
        LinearSolver
        
        % 牛顿迭代参数
        MaxIter = 50
        Tolerance = 1e-4
        
        % 线搜索参数
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        % 正则化/规范固定
        AutoRegularization = true 
        RegScaleRel = 1e-6  % 相对缩放因子
        MatrixK_Add = []    % 附加矩阵
    end
    
    methods
        function obj = HBFEMSolver(assembler, aft)
            obj.Assembler = assembler;
            obj.AFT = aft;
            
            % 初始化线性求解器
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 200; 
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [X_harmonics, info] = solve(obj, space, matLibData, sourceMaps, fixedDofs, x0)
            fprintf('==============================================\n');
            fprintf('   HBFEM Solver (v4.1 - Sparse Print Fix)\n');
            fprintf('==============================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            numDofs = dofHandler.NumGlobalDofs;
            numHarm = obj.AFT.NumHarmonics;
            totalDofs = numDofs * numHarm;
            
            % --- 1. 预计算正则化矩阵 (Gauge Fixing) ---
            M_reg_full = sparse(totalDofs, totalDofs);
            if obj.AutoRegularization
                fprintf('   [Init] Assembling mass matrix for regularization...\n');
                M_local = obj.Assembler.assembleMass(space);
                M_reg_full = kron(sparse(eye(numHarm)), M_local);
            end
            
            % --- 2. 组装全局源项 F (RHS) ---
            fprintf('   [Init] Assembling source vectors...\n');
            F_global = zeros(totalDofs, 1);
            
            if nargin >= 4 && ~isempty(sourceMaps)
                for h = 1:numHarm
                    map_h = [];
                    if iscell(sourceMaps)
                        if length(sourceMaps) >= h, map_h = sourceMaps{h}; end
                    elseif isa(sourceMaps, 'containers.Map')
                        if h == 1, map_h = sourceMaps; end
                    end
                    
                    if ~isempty(map_h)
                        F_h_local = obj.Assembler.assembleSource(space, map_h);
                        if length(F_h_local) < numDofs
                            F_h_local(numDofs, 1) = 0;
                        end
                        idx_start = (h-1) * numDofs + 1;
                        idx_end   = h * numDofs;
                        F_global(idx_start : idx_end) = F_h_local;
                    end
                end
            end
            
            % --- 3. 边界条件 ---
            if islogical(fixedDofs), idx_A = find(fixedDofs); else, idx_A = fixedDofs(:); end
            fixed_all = [];
            for h = 0:numHarm-1
                fixed_all = [fixed_all; idx_A + h*numDofs];
            end
            
            K_add = sparse(totalDofs, totalDofs);
            if ~isempty(obj.MatrixK_Add)
                K_add = obj.MatrixK_Add;
            end
            
            % --- 4. 初始解 ---
            if nargin < 6 || isempty(x0)
                X_mat = zeros(totalDofs, 1);
            else
                X_mat = x0;
            end
            
            converged = false;
            res_history = [];
            RegScaleAbs = 0;
            norm_base = 1.0;
            
            % --- 5. 牛顿迭代 ---
            for iter = 1:obj.MaxIter
                
                % 5.1 维度重塑 (Vector -> Matrix)
                X_input_matrix = reshape(X_mat, numDofs, numHarm);
                
                % 5.2 组装 (Jacobian & Residual)
                [J_triplets, R_mat] = obj.Assembler.assembleHBFEM(space, X_input_matrix, obj.AFT, matLibData, true);
                
                % 5.3 扁平化残差
                R_vec = R_mat(:); 
                
                % 稀疏 Jacobian
                H_mat = sparse(J_triplets.I, J_triplets.J, J_triplets.V, totalDofs, totalDofs);
                
                % 5.4 物理残差
                R_vec = R_vec - F_global;
                
                % 5.5 自动正则化 (Fixing Gauge)
                if obj.AutoRegularization
                    if RegScaleAbs == 0
                        diag_H = abs(diag(H_mat));
                        avg_diag = mean(diag_H(diag_H > 1e-12)); 
                        if isnan(avg_diag) || avg_diag == 0, avg_diag = 1.0; end
                        RegScaleAbs = avg_diag * obj.RegScaleRel;
                        
                        % [修复] 使用 full() 确保 fprintf 不报错
                        fprintf('   [Auto-Reg] Applied scale: %.2e\n', full(RegScaleAbs));
                    end
                    
                    H_mat = H_mat + RegScaleAbs * M_reg_full;
                    R_vec = R_vec + (RegScaleAbs * M_reg_full) * X_mat;
                end
                
                % 附加项
                H_mat = H_mat + K_add;
                R_vec = R_vec + K_add * X_mat;
                
                % 5.6 边界条件
                [H_bc, R_bc] = BoundaryCondition.applyDirichlet(H_mat, R_vec, fixed_all);
                
                current_res_norm = norm(R_bc);
                if iter == 1, norm_base = max(current_res_norm, 1e-6); end
                
                % 5.7 线性求解
                dX_mat = obj.LinearSolver.solve(H_bc, -R_bc);
                
                if isempty(dX_mat) || any(isnan(dX_mat))
                    warning('Linear solver failed (NaN). Aborting.');
                    break;
                end
                
                % 5.8 线搜索
                alpha = 1.0;
                accepted = false;
                
                for bt = 0:obj.MaxBacktracks
                    X_try = X_mat + alpha * dX_mat;
                    
                    % 快速计算残差
                    X_try_matrix = reshape(X_try, numDofs, numHarm);
                    [~, R_mat_try] = obj.Assembler.assembleHBFEM(space, X_try_matrix, obj.AFT, matLibData, false);
                    R_try_vec = R_mat_try(:);
                    
                    R_try_total = R_try_vec - F_global; 
                    if obj.AutoRegularization
                        R_try_total = R_try_total + (RegScaleAbs * M_reg_full) * X_try;
                    end
                    R_try_total = R_try_total + K_add * X_try;
                    
                    R_try_total(fixed_all) = 0;
                    new_res_norm = norm(R_try_total);
                    
                    if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                        X_mat = X_try;
                        current_res_norm = new_res_norm;
                        accepted = true;
                        
                        tag = ''; if bt>0, tag=sprintf('[BT %d]', bt); end
                        fprintf('   Iter %2d: Res=%.4e (Rel=%.4e) Step=%.4f %s\n', ...
                            iter, current_res_norm, current_res_norm/norm_base, alpha, tag);
                        break;
                    else
                        alpha = alpha * obj.BacktrackFactor;
                    end
                end
                
                if ~accepted
                    warning('Line search stuck. Forcing small step.');
                    X_mat = X_mat + 1e-3 * dX_mat; 
                end
                
                res_history(end+1) = current_res_norm;
                
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