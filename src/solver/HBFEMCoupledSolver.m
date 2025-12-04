classdef HBFEMCoupledSolver < handle
    % HBFEMCOUPLEDSOLVER 场路耦合谐波平衡求解器 (v1.0)
    % 
    % 功能:
    %   求解电压驱动的非线性稳态磁场问题 (HBFEM + Circuit)。
    %   支持任意次谐波电压激励，自动计算各次谐波电流。
    
    properties
        Assembler
        AFT
        LinearSolver
        
        WindingObj
        CircuitR
        
        MaxIter = 50
        Tolerance = 1e-4
        
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        MatrixK_Add
    end
    
    methods
        function obj = HBFEMCoupledSolver(assembler, aft, winding, R)
            obj.Assembler = assembler;
            obj.AFT = aft;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
        end
        
        function [X_sol, I_sol, info] = solve(obj, space, matLibData, V_harmonics, fixedDofs, x0)
            % SOLVE 执行耦合求解
            % 输入:
            %   V_harmonics - [NumHarm x 1] 各次谐波电压相量
            %   x0 - [NumDofs x NumHarm] 初值 (可选)
            % 输出:
            %   X_sol - [NumDofs x NumHarm] 磁矢位解
            %   I_sol - [NumHarm x 1] 电流解
            
            fprintf('==============================================\n');
            fprintf('   HBFEM Coupled Solver (Voltage Driven)      \n');
            fprintf('==============================================\n');
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            numHarm = obj.AFT.NumHarmonics;
            harmonics = obj.AFT.Harmonics;
            baseFreq = obj.AFT.BaseFreq;
            
            totalFieldDofs = numDofs * numHarm;
            totalCircuitDofs = numHarm;
            totalSystemSize = totalFieldDofs + totalCircuitDofs;
            
            % 1. 初始化
            if nargin < 6 || isempty(x0)
                A_mat = zeros(numDofs, numHarm);
            else
                A_mat = x0;
            end
            I_vec = zeros(numHarm, 1); % 电流谐波初值
            
            % 2. 预计算耦合向量 C
            C_vec = obj.Assembler.assembleWinding(space, obj.WindingObj);
            
            % 3. 准备正则化 (如果有)
            K_add_global = sparse(totalFieldDofs, totalFieldDofs);
            has_reg = ~isempty(obj.MatrixK_Add);
            if has_reg
                fprintf('   [Info] Applying Regularization.\n');
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
                K_add_global = sparse(I_reg, J_reg, V_reg, totalFieldDofs, totalFieldDofs);
            end
            
            % 4. 扩展边界条件
            % [A_harm1; ... A_harmK; I_1; ... I_K]
            % 电流自由度通常不固定
            Fixed_field = repmat(fixedDofs(:), numHarm, 1);
            Fixed_circuit = false(numHarm, 1); 
            Fixed_sys = [Fixed_field; Fixed_circuit];
            free_mask = ~Fixed_sys;
            
            % 5. 初始残差计算
            % 场方程: K(A)*A_k - C*I_k = 0
            % 路方程: j*w_k*C'*A_k + R*I_k - V_k = 0
            
            [~, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, false);
            R_mag_vec = R_mag_mat(:);
            
            % 场方程残差需要减去电流源项 C*I
            % 我们需要构建一个包含 I 贡献的向量
            F_current = zeros(totalFieldDofs, 1);
            for k = 1:numHarm
                idx_s = (k-1)*numDofs + 1;
                idx_e = k*numDofs;
                F_current(idx_s:idx_e) = C_vec * I_vec(k);
            end
            
            Res_Field = R_mag_vec - F_current;
            if has_reg, Res_Field = Res_Field + K_add_global * A_mat(:); end
            
            Res_Circuit = zeros(numHarm, 1);
            for k = 1:numHarm
                h = harmonics(k);
                w = 2 * pi * baseFreq * h;
                Psi_k = C_vec' * A_mat(:, k);
                % V = RI + jwPsi
                Res_Circuit(k) = obj.CircuitR * I_vec(k) + 1j * w * Psi_k - V_harmonics(k);
            end
            
            % --- Auto-Scaling ---
            mag_scale = mean(abs(Res_Field)); if mag_scale < 1e-5, mag_scale = 1.0; end
            cir_scale_base = abs(obj.CircuitR); if cir_scale_base < 1e-10, cir_scale_base = 1.0; end
            row_scale = mag_scale / cir_scale_base;
            row_scale = min(max(row_scale, 1e3), 1e9); % 限制范围
            
            fprintf('   [Scaling] Circuit Row Scale Factor: %.2e\n', row_scale);
            
            Total_Res = [Res_Field; Res_Circuit * row_scale];
            current_res_norm = norm(Total_Res(free_mask));
            norm_base = current_res_norm; if norm_base < 1e-5, norm_base = 1.0; end
            
            fprintf('   Iter  0: Res=%.4e (Rel=%.4e)\n', current_res_norm, current_res_norm/norm_base);
            
            res_history = current_res_norm;
            converged = false;
            
            % 6. 迭代循环
            for iter = 1:obj.MaxIter
                % A. 组装 Jacobian (Field Part)
                [J_tri, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, true);
                J_field = sparse(J_tri.I, J_tri.J, J_tri.V, totalFieldDofs, totalFieldDofs);
                
                if has_reg, J_field = J_field + K_add_global; end
                
                % B. 构建耦合矩阵块
                % Top-Right: -C (Diagonal blocks)
                % Bot-Left:  jwC' (Diagonal blocks)
                % Bot-Right: R (Diagonal)
                
                % 使用稀疏三元组构建耦合块
                % Top-Right (-C)
                [ci, ~, cv] = find(C_vec);
                n_c = length(cv);
                TR_I = zeros(n_c * numHarm, 1);
                TR_J = zeros(n_c * numHarm, 1);
                TR_V = zeros(n_c * numHarm, 1);
                
                % Bot-Left (jwC')
                BL_I = zeros(n_c * numHarm, 1);
                BL_J = zeros(n_c * numHarm, 1);
                BL_V = zeros(n_c * numHarm, 1);
                
                % Bot-Right (R)
                BR_I = (1:numHarm)';
                BR_J = (1:numHarm)';
                BR_V = ones(numHarm, 1) * obj.CircuitR * row_scale;
                
                for k = 1:numHarm
                    h = harmonics(k);
                    w = 2 * pi * baseFreq * h;
                    
                    offset_field_row = (k-1)*numDofs;
                    col_idx_circuit = k; % Relative to circuit block start
                    
                    idx_rng = (k-1)*n_c + 1 : k*n_c;
                    
                    % Top-Right: Row=Field, Col=Circuit. Val = -C
                    TR_I(idx_rng) = ci + offset_field_row;
                    TR_J(idx_rng) = col_idx_circuit;
                    TR_V(idx_rng) = -cv;
                    
                    % Bot-Left: Row=Circuit, Col=Field. Val = jwC' * scale
                    BL_I(idx_rng) = col_idx_circuit;
                    BL_J(idx_rng) = ci + offset_field_row;
                    BL_V(idx_rng) = 1j * w * cv * row_scale;
                end
                
                % 拼接大矩阵
                J_TR = sparse(TR_I, TR_J, TR_V, totalFieldDofs, totalCircuitDofs);
                J_BL = sparse(BL_I, BL_J, BL_V, totalCircuitDofs, totalFieldDofs);
                J_BR = sparse(BR_I, BR_J, BR_V, totalCircuitDofs, totalCircuitDofs);
                
                J_sys = [J_field, J_TR; J_BL, J_BR];
                
                % C. 计算残差
                R_mag_vec = R_mag_mat(:);
                
                % 更新 F_current (基于当前 I_vec)
                F_current = zeros(totalFieldDofs, 1);
                for k = 1:numHarm
                    idx_s = (k-1)*numDofs + 1;
                    idx_e = k*numDofs;
                    F_current(idx_s:idx_e) = C_vec * I_vec(k);
                end
                
                Res_Field = R_mag_vec - F_current;
                if has_reg, Res_Field = Res_Field + K_add_global * A_mat(:); end
                
                Res_Circuit = zeros(numHarm, 1);
                for k = 1:numHarm
                    h = harmonics(k);
                    w = 2 * pi * baseFreq * h;
                    Psi_k = C_vec' * A_mat(:, k);
                    Res_Circuit(k) = obj.CircuitR * I_vec(k) + 1j * w * Psi_k - V_harmonics(k);
                end
                
                Total_Res = [Res_Field; Res_Circuit * row_scale];
                
                % D. 求解
                RHS = -Total_Res;
                [J_solve, R_solve] = BoundaryCondition.applyDirichlet(J_sys, RHS, Fixed_sys);
                
                dX_sys = obj.LinearSolver.solve(J_solve, R_solve);
                
                dA_vec = dX_sys(1:totalFieldDofs);
                dI_vec = dX_sys(totalFieldDofs+1:end);
                
                dA_mat = reshape(dA_vec, numDofs, numHarm);
                
                % E. 线搜索
                alpha = 1.0;
                accepted = false;
                
                for bt = 0:obj.MaxBacktracks
                    A_try = A_mat + alpha * dA_mat;
                    I_try = I_vec + alpha * dI_vec;
                    
                    % Fast Residual
                    [~, R_mag_try] = obj.Assembler.assembleHBFEM(space, A_try, obj.AFT, matLibData, false);
                    R_vec_try = R_mag_try(:);
                    
                    F_curr_try = zeros(totalFieldDofs, 1);
                    for k = 1:numHarm
                        idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                        F_curr_try(idx_s:idx_e) = C_vec * I_try(k);
                    end
                    
                    Res_F_try = R_vec_try - F_curr_try;
                    if has_reg, Res_F_try = Res_F_try + K_add_global * A_try(:); end
                    
                    Res_C_try = zeros(numHarm, 1);
                    for k = 1:numHarm
                        h = harmonics(k); w = 2*pi*baseFreq*h;
                        Psi_k = C_vec' * A_try(:, k);
                        Res_C_try(k) = obj.CircuitR * I_try(k) + 1j * w * Psi_k - V_harmonics(k);
                    end
                    
                    Total_Res_try = [Res_F_try; Res_C_try * row_scale];
                    new_res_norm = norm(Total_Res_try(free_mask));
                    
                    if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                        A_mat = A_try;
                        I_vec = I_try;
                        current_res_norm = new_res_norm;
                        accepted = true;
                        
                        marker = ''; if bt>0, marker=sprintf('[BT %d]', bt); end
                        fprintf('   Iter %d: Res=%.2e (Rel=%.2e) Step=%.4f %s\n', ...
                            iter, current_res_norm, current_res_norm/norm_base, alpha, marker);
                        break;
                    else
                        alpha = alpha * obj.BacktrackFactor;
                    end
                end
                
                if ~accepted
                    warning('Line search stuck. Forcing step.');
                    A_mat = A_mat + 1e-3 * dA_mat;
                    I_vec = I_vec + 1e-3 * dI_vec;
                end
                
                if current_res_norm / norm_base < obj.Tolerance
                    converged = true;
                    break;
                end
            end
            
            X_sol = A_mat;
            I_sol = I_vec;
            info.Iterations = iter;
            info.Converged = converged;
            fprintf('==============================================\n');
        end
    end
end