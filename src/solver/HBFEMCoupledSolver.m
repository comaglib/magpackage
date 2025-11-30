classdef HBFEMCoupledSolver < handle
    % HBFEMCOUPLEDSOLVER 场路耦合谐波平衡求解器 (v2.0 - Auto Reg & Inductance)
    % 
    % 更新日志:
    %   1. [Feature] 增加电路漏感 L (CircuitL) 支持。
    %   2. [Robustness] 内置自动正则化 (Gauge Fixing)，无需外部手动注入。
    
    properties
        Assembler
        AFT
        LinearSolver
        
        WindingObj
        CircuitR
        CircuitL % [New] 外部电感
        
        MaxIter = 50
        Tolerance = 1e-4
        
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        % 正则化配置
        AutoRegularization = true 
        RegScaleRel = 1e-6
        MatrixK_Add % 保留用于额外的特殊刚度修正
    end
    
    methods
        function obj = HBFEMCoupledSolver(assembler, aft, winding, R, L)
            obj.Assembler = assembler;
            obj.AFT = aft;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            if nargin < 5, L = 0; end
            obj.CircuitL = L;
            
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
            obj.LinearSolver.ReuseAnalysis = false;
        end
        
        function [X_sol, I_sol, info] = solve(obj, space, matLibData, V_harmonics, fixedDofs, x0)
            % SOLVE 执行耦合求解
            
            fprintf('==============================================\n');
            fprintf('   HBFEM Coupled Solver (v2.0 - R-L Circuit)  \n');
            fprintf('==============================================\n');
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            numHarm = obj.AFT.NumHarmonics;
            harmonics = obj.AFT.Harmonics;
            baseFreq = obj.AFT.BaseFreq;
            
            totalFieldDofs = numDofs * numHarm;
            totalCircuitDofs = numHarm;
            
            % 1. 初始化
            if nargin < 6 || isempty(x0)
                A_mat = zeros(numDofs, numHarm);
            else
                A_mat = x0;
            end
            I_vec = zeros(numHarm, 1); 
            
            % 2. 预计算耦合向量 C
            C_vec = obj.Assembler.assembleWinding(space, obj.WindingObj);
            
            % 3. 准备正则化矩阵 (M_reg_full)
            M_reg_full = sparse(totalFieldDofs, totalFieldDofs);
            if obj.AutoRegularization
                fprintf('   [Init] Assembling mass matrix for regularization...\n');
                M_local = obj.Assembler.assembleMass(space);
                % 块对角矩阵 kron(I, M)
                M_reg_full = kron(sparse(eye(numHarm)), M_local);
            end
            
            % 准备额外的 K_Add (如果用户手动指定)
            K_add_global = sparse(totalFieldDofs, totalFieldDofs);
            if ~isempty(obj.MatrixK_Add)
                M_small = obj.MatrixK_Add;
                [mi, mj, mv] = find(M_small);
                n_nz = length(mv);
                I_reg = zeros(n_nz * numHarm, 1); J_reg = zeros(n_nz * numHarm, 1); V_reg = zeros(n_nz * numHarm, 1);
                for h = 1:numHarm
                    offset = (h-1) * numDofs; idx_range = (h-1)*n_nz + 1 : h*n_nz;
                    I_reg(idx_range) = mi + offset; J_reg(idx_range) = mj + offset; V_reg(idx_range) = mv;
                end
                K_add_global = sparse(I_reg, J_reg, V_reg, totalFieldDofs, totalFieldDofs);
            end
            
            % 4. 扩展边界条件
            Fixed_field = repmat(fixedDofs(:), numHarm, 1);
            Fixed_circuit = false(numHarm, 1); 
            Fixed_sys = [Fixed_field; Fixed_circuit];
            free_mask = ~Fixed_sys;
            
            % 5. 初始残差计算
            [~, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, false);
            R_mag_vec = R_mag_mat(:);
            
            F_current = zeros(totalFieldDofs, 1);
            for k = 1:numHarm
                idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                F_current(idx_s:idx_e) = C_vec * I_vec(k);
            end
            
            Res_Field = R_mag_vec - F_current + K_add_global * A_mat(:);
            
            Res_Circuit = zeros(numHarm, 1);
            for k = 1:numHarm
                h = harmonics(k); w = 2 * pi * baseFreq * h;
                Psi_k = C_vec' * A_mat(:, k);
                % [Update] V = Z*I + jwPsi, Z = R + jwL
                Z_k = obj.CircuitR + 1j * w * obj.CircuitL;
                Res_Circuit(k) = Z_k * I_vec(k) + 1j * w * Psi_k - V_harmonics(k);
            end
            
            % --- Auto-Scaling ---
            mag_scale = mean(abs(Res_Field)); if mag_scale < 1e-5, mag_scale = 1.0; end
            cir_scale_base = abs(obj.CircuitR) + abs(2*pi*baseFreq*obj.CircuitL); 
            if cir_scale_base < 1e-10, cir_scale_base = 1.0; end
            row_scale = mag_scale / cir_scale_base;
            row_scale = min(max(row_scale, 1e3), 1e9); 
            
            fprintf('   [Scaling] Circuit Row Scale Factor: %.2e\n', row_scale);
            
            Total_Res = [Res_Field; Res_Circuit * row_scale];
            current_res_norm = norm(Total_Res(free_mask));
            norm_base = current_res_norm; if norm_base < 1e-5, norm_base = 1.0; end
            
            fprintf('   Iter  0: Res=%.4e (Rel=%.4e)\n', current_res_norm, current_res_norm/norm_base);
            
            RegScaleAbs = 0; % 自动正则化系数
            converged = false;
            
            % 6. 迭代循环
            for iter = 1:obj.MaxIter
                % A. 组装 Jacobian (Field Part)
                [J_tri, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, true);
                J_field = sparse(J_tri.I, J_tri.J, J_tri.V, totalFieldDofs, totalFieldDofs);
                J_field = J_field + K_add_global;
                
                % [Auto-Reg] 计算并应用正则化
                if obj.AutoRegularization
                    if RegScaleAbs == 0
                        diag_J = abs(diag(J_field));
                        avg_diag = mean(diag_J(diag_J > 1e-12)); 
                        if isnan(avg_diag) || avg_diag == 0, avg_diag = 1.0; end
                        RegScaleAbs = avg_diag * obj.RegScaleRel;
                        fprintf('   [Auto-Reg] Applied scale: %.2e\n', full(RegScaleAbs));
                    end
                    J_field = J_field + RegScaleAbs * M_reg_full;
                end
                
                % B. 构建耦合矩阵块
                [ci, ~, cv] = find(C_vec); n_c = length(cv);
                
                % Top-Right (-C)
                TR_I = zeros(n_c * numHarm, 1); TR_J = zeros(n_c * numHarm, 1); TR_V = zeros(n_c * numHarm, 1);
                % Bot-Left (jwC')
                BL_I = zeros(n_c * numHarm, 1); BL_J = zeros(n_c * numHarm, 1); BL_V = zeros(n_c * numHarm, 1);
                % Bot-Right (Z)
                BR_I = (1:numHarm)'; BR_J = (1:numHarm)'; BR_V = zeros(numHarm, 1);
                
                for k = 1:numHarm
                    h = harmonics(k); w = 2 * pi * baseFreq * h;
                    offset_field_row = (k-1)*numDofs;
                    
                    idx_rng = (k-1)*n_c + 1 : k*n_c;
                    
                    % Top-Right: -C
                    TR_I(idx_rng) = ci + offset_field_row; TR_J(idx_rng) = k; TR_V(idx_rng) = -cv;
                    
                    % Bot-Left: jwC' * scale
                    BL_I(idx_rng) = k; BL_J(idx_rng) = ci + offset_field_row; BL_V(idx_rng) = 1j * w * cv * row_scale;
                    
                    % Bot-Right: Z * scale
                    Z_k = obj.CircuitR + 1j * w * obj.CircuitL;
                    BR_V(k) = Z_k * row_scale;
                end
                
                J_TR = sparse(TR_I, TR_J, TR_V, totalFieldDofs, totalCircuitDofs);
                J_BL = sparse(BL_I, BL_J, BL_V, totalCircuitDofs, totalFieldDofs);
                J_BR = sparse(BR_I, BR_J, BR_V, totalCircuitDofs, totalCircuitDofs);
                J_sys = [J_field, J_TR; J_BL, J_BR];
                
                % C. 计算残差
                R_mag_vec = R_mag_mat(:);
                F_current = zeros(totalFieldDofs, 1);
                for k = 1:numHarm
                    idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                    F_current(idx_s:idx_e) = C_vec * I_vec(k);
                end
                
                Res_Field = R_mag_vec - F_current + K_add_global * A_mat(:);
                if obj.AutoRegularization
                    Res_Field = Res_Field + (RegScaleAbs * M_reg_full) * A_mat(:);
                end
                
                Res_Circuit = zeros(numHarm, 1);
                for k = 1:numHarm
                    h = harmonics(k); w = 2 * pi * baseFreq * h;
                    Psi_k = C_vec' * A_mat(:, k);
                    Z_k = obj.CircuitR + 1j * w * obj.CircuitL;
                    Res_Circuit(k) = Z_k * I_vec(k) + 1j * w * Psi_k - V_harmonics(k);
                end
                Total_Res = [Res_Field; Res_Circuit * row_scale];
                
                % D. 求解
                [J_solve, R_solve] = BoundaryCondition.applyDirichlet(J_sys, -Total_Res, Fixed_sys);
                dX_sys = obj.LinearSolver.solve(J_solve, R_solve);
                if isempty(dX_sys), error('Linear Solver Failed'); end
                
                dA_vec = dX_sys(1:totalFieldDofs);
                dI_vec = dX_sys(totalFieldDofs+1:end);
                dA_mat = reshape(dA_vec, numDofs, numHarm);
                
                % E. 线搜索
                alpha = 1.0; accepted = false;
                for bt = 0:obj.MaxBacktracks
                    A_try = A_mat + alpha * dA_mat;
                    I_try = I_vec + alpha * dI_vec;
                    
                    [~, R_mag_try] = obj.Assembler.assembleHBFEM(space, A_try, obj.AFT, matLibData, false);
                    R_vec_try = R_mag_try(:);
                    
                    F_curr_try = zeros(totalFieldDofs, 1);
                    for k = 1:numHarm
                        idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                        F_curr_try(idx_s:idx_e) = C_vec * I_try(k);
                    end
                    
                    Res_F_try = R_vec_try - F_curr_try + K_add_global * A_try(:);
                    if obj.AutoRegularization
                        Res_F_try = Res_F_try + (RegScaleAbs * M_reg_full) * A_try(:);
                    end
                    
                    Res_C_try = zeros(numHarm, 1);
                    for k = 1:numHarm
                        h = harmonics(k); w = 2*pi*baseFreq*h;
                        Psi_k = C_vec' * A_try(:, k);
                        Z_k = obj.CircuitR + 1j * w * obj.CircuitL;
                        Res_C_try(k) = Z_k * I_try(k) + 1j * w * Psi_k - V_harmonics(k);
                    end
                    
                    Total_Res_try = [Res_F_try; Res_C_try * row_scale];
                    new_res_norm = norm(Total_Res_try(free_mask));
                    
                    if new_res_norm < current_res_norm || alpha < obj.MinStepSize
                        A_mat = A_try; I_vec = I_try;
                        current_res_norm = new_res_norm; accepted = true;
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
                    A_mat = A_mat + 1e-3 * dA_mat; I_vec = I_vec + 1e-3 * dI_vec;
                end
                
                if current_res_norm / norm_base < obj.Tolerance, converged = true; break; end
            end
            
            X_sol = A_mat; I_sol = I_vec;
            info.Iterations = iter; info.Converged = converged;
            fprintf('==============================================\n');
        end
    end
end