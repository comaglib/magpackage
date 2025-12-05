classdef DFTFEMCoupledSolver < handle
    % DFTFEMCoupledSolver 场路耦合离散傅里叶变换求解器
    % 
    % 描述:
    %   基于 HBFEMCoupledSolver 开发。
    %   专用于全频谱 (Full Spectrum) 求解，即时域采样点数与频域自由度数匹配的情况。
    %   此时变换矩阵是方阵 (或满秩矩形阵)，等效于在每个时间步满足方程。
    %
    %   适用场景:
    %   - 波形畸变极其严重，难以预知主要谐波次数的情况。
    %   - 需要极高精度的非线性时域重构。
    %
    %   注意: 
    %   相比稀疏 HBFEM，DFTFEM 的系统矩阵更稠密，内存消耗显著增加。
    
    properties
        Assembler       % 有限元组装器
        AFT             % 谐波-时域转换器 (应配置为全谐波模式)
        LinearSolver    % 线性求解器
        
        WindingObj      % 绕组对象
        CircuitR        % 外电路电阻
        
        % --- 迭代控制 ---
        MaxIter = 50
        Tolerance = 1e-4
        
        % --- 线搜索 ---
        MaxBacktracks = 20
        BacktrackFactor = 0.5
        MinStepSize = 1e-6
        
        % --- 正则化 ---
        MatrixK_Add     % 用于处理 Sigma=0 区域的奇异性
    end
    
    methods
        function obj = DFTFEMCoupledSolver(assembler, aft, winding, R)
            % 构造函数
            obj.Assembler = assembler;
            obj.AFT = aft;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            
            % 初始化线性求解器
            % DFTFEM 产生的矩阵规模较大且非对称，需优化配置
            obj.LinearSolver = LinearSolver('Auto');
            % 必须设置为非对称模式 (0)
            obj.LinearSolver.MumpsSymmetry = 0; 
            % DFTFEM 填充率高，大幅增加内存预估 (400%)
            obj.LinearSolver.MumpsICNTL.i14 = 400; 
        end
        
        function [X_sol, I_sol, info] = solve(obj, space, matLibData, V_spectrum, fixedDofs, x0)
            % SOLVE 执行 DFT-FEM 求解
            %
            % 输入:
            %   V_spectrum: [NumHarm x 1] 电压源的频谱分量 (全频谱)
            
            fprintf('==============================================\n');
            fprintf('   DFT-FEM Coupled Solver (Full Spectrum)     \n');
            fprintf('==============================================\n');
            
            % 1. 检查 DFTFEM 条件 (时频点数匹配)
            % -----------------------------------------------------------
            numHarm = obj.AFT.NumHarmonics;
            numTimeSteps = obj.AFT.NumTimeSteps;
            maxHarm = max(obj.AFT.Harmonics);
            
            % 理论上，实数 DFT 需要 N_t >= 2*max(H) + 1
            % 如果由 AFT 自动生成，通常是满足的。这里仅做提示。
            fprintf('   [Check] Harmonics: %d, TimeSteps: %d\n', numHarm, numTimeSteps);
            if numTimeSteps > 4 * numHarm
                fprintf(['   DFTFEM: Oversampling: ' ...
                    '时域采样点远多于频域未知数，这更像是 HBFEM 而非标准 DFTFEM。建议检查 AFT 配置。\n']);
            end
            
            % 2. 基础初始化
            % -----------------------------------------------------------
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            totalFieldDofs = numDofs * numHarm;
            totalCircuitDofs = numHarm;
            
            if nargin < 6 || isempty(x0)
                A_mat = zeros(numDofs, numHarm);
            else
                A_mat = x0;
            end
            I_vec = zeros(numHarm, 1);
            
            % 3. 预计算常数向量 C (线圈几何)
            C_vec = obj.Assembler.assembleWinding(space, obj.WindingObj);
            
            % 4. 准备正则化矩阵 (Robust Regularization)
            % -----------------------------------------------------------
            K_add_global = sparse(totalFieldDofs, totalFieldDofs);
            has_reg = ~isempty(obj.MatrixK_Add);
            if has_reg
                fprintf('   [Info] Applying Regularization to eliminate singularity.\n');
                M_small = obj.MatrixK_Add;
                [mi, mj, mv] = find(M_small);
                n_nz = length(mv);
                
                % 扩展到所有频率块对角线
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
            
            % 5. 扩展边界条件
            Fixed_field = repmat(fixedDofs(:), numHarm, 1);
            Fixed_circuit = false(numHarm, 1); 
            Fixed_sys = [Fixed_field; Fixed_circuit];
            free_mask = ~Fixed_sys;
            
            % 6. 初始残差与自动缩放 (Auto-Scaling)
            % -----------------------------------------------------------
            % 场方程残差
            [~, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, false);
            R_mag_vec = R_mag_mat(:);
            
            F_current = zeros(totalFieldDofs, 1);
            for k = 1:numHarm
                idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                F_current(idx_s:idx_e) = C_vec * I_vec(k);
            end
            
            Res_Field = R_mag_vec - F_current;
            if has_reg, Res_Field = Res_Field + K_add_global * A_mat(:); end
            
            % 电路方程残差
            Res_Circuit = zeros(numHarm, 1);
            harmonics_list = obj.AFT.Harmonics;
            baseFreq = obj.AFT.BaseFreq;
            
            for k = 1:numHarm
                h = harmonics_list(k);
                w = 2 * pi * baseFreq * h;
                Psi_k = C_vec' * A_mat(:, k);
                Res_Circuit(k) = obj.CircuitR * I_vec(k) + 1j * w * Psi_k - V_spectrum(k);
            end
            
            % 计算缩放因子
            mag_scale = mean(abs(Res_Field)); if mag_scale < 1e-5, mag_scale = 1.0; end
            cir_scale_base = abs(obj.CircuitR); if cir_scale_base < 1e-10, cir_scale_base = 1.0; end
            row_scale = mag_scale / cir_scale_base;
            row_scale = min(max(row_scale, 1e3), 1e9); 
            
            fprintf('   [Scaling] Circuit Row Scale Factor: %.2e\n', row_scale);
            
            Total_Res = [Res_Field; Res_Circuit * row_scale];
            current_res_norm = norm(Total_Res(free_mask));
            norm_base = current_res_norm; if norm_base < 1e-5, norm_base = 1.0; end
            
            fprintf('   Iter  0: Res=%.4e (Rel=%.4e)\n', current_res_norm, current_res_norm/norm_base);
            
            converged = false;
            
            % 7. 牛顿迭代循环
            % -----------------------------------------------------------
            for iter = 1:obj.MaxIter
                % A. 组装场雅可比 (Dense Blocks in Frequency Domain)
                [J_tri, R_mag_mat] = obj.Assembler.assembleHBFEM(space, A_mat, obj.AFT, matLibData, true);
                J_field = sparse(J_tri.I, J_tri.J, J_tri.V, totalFieldDofs, totalFieldDofs);
                
                if has_reg, J_field = J_field + K_add_global; end
                
                % B. 构建耦合块
                [ci, ~, cv] = find(C_vec);
                n_c = length(cv);
                
                TR_I = zeros(n_c * numHarm, 1); TR_J = zeros(n_c * numHarm, 1); TR_V = zeros(n_c * numHarm, 1);
                BL_I = zeros(n_c * numHarm, 1); BL_J = zeros(n_c * numHarm, 1); BL_V = zeros(n_c * numHarm, 1);
                BR_I = (1:numHarm)'; BR_J = (1:numHarm)'; BR_V = ones(numHarm, 1) * obj.CircuitR * row_scale;
                
                for k = 1:numHarm
                    h = harmonics_list(k);
                    w = 2 * pi * baseFreq * h;
                    
                    offset_field = (k-1)*numDofs;
                    col_idx = k;
                    idx_rng = (k-1)*n_c + 1 : k*n_c;
                    
                    % Top-Right: -C
                    TR_I(idx_rng) = ci + offset_field; TR_J(idx_rng) = col_idx; TR_V(idx_rng) = -cv;
                    % Bot-Left: jwC' * scale
                    BL_I(idx_rng) = col_idx; BL_J(idx_rng) = ci + offset_field; BL_V(idx_rng) = 1j * w * cv * row_scale;
                end
                
                J_TR = sparse(TR_I, TR_J, TR_V, totalFieldDofs, totalCircuitDofs);
                J_BL = sparse(BL_I, BL_J, BL_V, totalCircuitDofs, totalFieldDofs);
                J_BR = sparse(BR_I, BR_J, BR_V, totalCircuitDofs, totalCircuitDofs);
                
                J_sys = [J_field, J_TR; J_BL, J_BR];
                
                % C. 更新残差
                R_mag_vec = R_mag_mat(:);
                F_current = zeros(totalFieldDofs, 1);
                for k = 1:numHarm
                    idx_s = (k-1)*numDofs + 1; idx_e = k*numDofs;
                    F_current(idx_s:idx_e) = C_vec * I_vec(k);
                end
                Res_Field = R_mag_vec - F_current;
                if has_reg, Res_Field = Res_Field + K_add_global * A_mat(:); end
                
                Res_Circuit = zeros(numHarm, 1);
                for k = 1:numHarm
                    h = harmonics_list(k); w = 2 * pi * baseFreq * h;
                    Psi_k = C_vec' * A_mat(:, k);
                    Res_Circuit(k) = obj.CircuitR * I_vec(k) + 1j * w * Psi_k - V_spectrum(k);
                end
                Total_Res = [Res_Field; Res_Circuit * row_scale];
                
                % D. 求解
                RHS = -Total_Res;
                [J_solve, R_solve] = BoundaryCondition.applyDirichlet(J_sys, RHS, Fixed_sys);
                dX_sys = obj.LinearSolver.solve(J_solve, R_solve);
                
                dA_vec = dX_sys(1:totalFieldDofs);
                dI_vec = dX_sys(totalFieldDofs+1:end);
                dA_mat = reshape(dA_vec, numDofs, numHarm);
                
                % E. 线搜索 (Line Search)
                alpha = 1.0;
                accepted = false;
                
                for bt = 0:obj.MaxBacktracks
                    A_try = A_mat + alpha * dA_mat;
                    I_try = I_vec + alpha * dI_vec;
                    
                    % Fast Residual Check
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
                        h = harmonics_list(k); w = 2*pi*baseFreq*h;
                        Psi_k = C_vec' * A_try(:, k);
                        Res_C_try(k) = obj.CircuitR * I_try(k) + 1j * w * Psi_k - V_spectrum(k);
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