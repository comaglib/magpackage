classdef EnvelopeCoupledSolver < handle
    % ENVELOPECOUPLEDSOLVER 包络有限元场路耦合求解器 (v2.1 - Verbose Trace)
    %
    % 核心升级: 实现了自适应 Levenberg-Marquardt (LM) 算法。
    % 功能增强: 添加了详细的每一步运行日志 (Verbose Logging)。
    %
    % 修复记录 (v2.1): 
    %   1. 修复了 assembleSystem 函数中 row_scale 输出参数未赋值的 Bug。
    %   2. 确保缩放因子在残差R和雅可比J中都保持一致。
    
    properties
        Assembler       
        AFT             
        LinearSolver    
        
        WindingObj      
        CircuitR        
        
        % --- 迭代控制 ---
        MaxIter = 25        
        Tolerance = 1e-4    
        
        % --- LM 算法参数 ---
        LambdaInit = 1e-2   
        LambdaMin = 1e-8    
        LambdaMax = 1e10    
        ScaleUp = 10.0      
        ScaleDown = 0.1     
        
        % --- 质量矩阵 (必须) ---
        MatrixM         
    end
    
    methods
        function obj = EnvelopeCoupledSolver(assembler, aft, winding, R)
            obj.Assembler = assembler;
            obj.AFT = aft;
            obj.WindingObj = winding;
            obj.CircuitR = R;
            
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 200; 
        end
        
        function [Results, info] = solve(obj, space, matLibData, timePoints, V_func, fixedDofs, x_init)
            fprintf('==============================================\n');
            fprintf('   Envelope FEM Solver (Levenberg-Marquardt)  \n');
            fprintf('==============================================\n');
            
            if isempty(obj.MatrixM), error('EnvelopeSolver:NoMass', 'Missing MatrixM.'); end
            
            % --- 1. 维度与初始化 ---
            fprintf('[Step 0] Initializing Solver Context...\n');
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            numHarm = obj.AFT.NumHarmonics;
            harmonics = obj.AFT.Harmonics;
            baseFreq = obj.AFT.BaseFreq;
            
            totalFieldDofs = numDofs * numHarm;
            fprintf('   -> System Size: Field=%d, Circuit=%d, Total=%d\n', ...
                totalFieldDofs, numHarm, totalFieldDofs+numHarm);
            
            if nargin < 7 || isempty(x_init)
                X_curr_mat = zeros(numDofs, numHarm);
                I_curr_vec = zeros(numHarm, 1);
            else
                X_curr_mat = x_init; 
                I_curr_vec = zeros(numHarm, 1); 
            end
            X_curr_flat = X_curr_mat(:);
            
            % --- 2. 预计算常数矩阵 ---
            fprintf('[Step 1] Pre-assembling Invariant Matrices...\n');
            
            fprintf('   -> Assembling Winding Coupling Vector (C)...\n');
            C_vec_small = obj.Assembler.assembleWinding(space, obj.WindingObj);
            
            fprintf('   -> Building Block-Diagonal Mass Matrix (M_sys)...\n');
            M_small = obj.MatrixM;
            [mi, mj, mv] = find(M_small);
            n_nz = length(mv);
            I_reg = zeros(n_nz * numHarm, 1); J_reg = zeros(n_nz * numHarm, 1); V_reg = zeros(n_nz * numHarm, 1);
            for h = 1:numHarm
                offset = (h-1) * numDofs;
                idx = (h-1)*n_nz + 1 : h*n_nz;
                I_reg(idx) = mi + offset; J_reg(idx) = mj + offset; V_reg(idx) = mv;
            end
            M_sys = sparse(I_reg, J_reg, V_reg, totalFieldDofs, totalFieldDofs);
            
            fprintf('   -> Building Eddy Damping Matrix (jw*M)...\n');
            V_jed = zeros(n_nz * numHarm, 1);
            for k = 1:numHarm
                w = 2 * pi * baseFreq * harmonics(k);
                idx = (k-1)*n_nz + 1 : k*n_nz;
                V_jed(idx) = 1j * w * mv;
            end
            J_eddy_sys = sparse(I_reg, J_reg, V_jed, totalFieldDofs, totalFieldDofs);
            
            % --- 3. 边界条件 ---
            Fixed_field = repmat(fixedDofs(:), numHarm, 1);
            Fixed_circuit = false(numHarm, 1);
            Fixed_sys = [Fixed_field; Fixed_circuit];
            free_mask = ~Fixed_sys;
            
            % --- 4. 时间步进循环 ---
            Results.Time = timePoints;
            Results.I_history = zeros(length(timePoints), numHarm); 
            X_prev_flat = X_curr_flat;
            
            % LM 状态变量
            lambda = obj.LambdaInit; 
            
            fprintf('[Step 2] Starting Simulation (%d steps)...\n', length(timePoints));
            
            for t_idx = 1:length(timePoints)
                t_curr = timePoints(t_idx);
                if t_idx == 1, dt = timePoints(2)-timePoints(1); else, dt = t_curr-timePoints(t_idx-1); end
                
                V_curr = V_func(t_curr);
                
                fprintf('\n--- Time Step %d/%d (t=%.4fs, dt=%.4fs) ---\n', t_idx, length(timePoints), t_curr, dt);
                
                % ==================================================
                % Levenberg-Marquardt 迭代循环
                % ==================================================
                
                % 初始残差计算 (需要 J 和 scale)
                % fprintf('   [Init] Calculating initial residual...\n');
                [R_curr, J_curr, ~] = obj.assembleSystem(space, matLibData, ...
                    X_curr_mat, I_curr_vec, X_prev_flat, V_curr, dt, ...
                    C_vec_small, M_sys, J_eddy_sys, numDofs, numHarm, harmonics, baseFreq, true);
                
                res_norm = norm(R_curr(free_mask));
                if t_idx==1, norm_base = res_norm; if norm_base<1e-6, norm_base=1; end; end
                
                fprintf('   [Start] Init Res: %.4e (Rel: %.4e)\n', res_norm, res_norm/norm_base);
                
                for iter = 1:obj.MaxIter
                    if res_norm / norm_base < obj.Tolerance
                        fprintf('      -> Converged (Tol met).\n');
                        break;
                    end
                    
                    % --- LM 内循环 (寻找可接受步长) ---
                    step_accepted = false;
                    
                    % 提取对角线用于缩放阻尼 (LM核心)
                    diag_J = sparse(1:size(J_curr,1), 1:size(J_curr,1), abs(diag(J_curr)));
                    
                    inner_iter = 0;
                    while ~step_accepted
                        inner_iter = inner_iter + 1;
                        
                        % 1. 构造阻尼系统: (J + lambda * diag(J)) * dx = -R
                        J_damp = J_curr + lambda * diag_J;
                        
                        % 2. 求解更新方向
                        [J_solve, R_solve] = BoundaryCondition.applyDirichlet(J_damp, -R_curr, Fixed_sys);
                        dX_sol = obj.LinearSolver.solve(J_solve, R_solve);
                        
                        % 3. 试探更新
                        dA = dX_sol(1:totalFieldDofs);
                        dI = dX_sol(totalFieldDofs+1:end);
                        
                        X_try_flat = X_curr_flat + dA;
                        I_try_vec  = I_curr_vec  + dI;
                        X_try_mat  = reshape(X_try_flat, numDofs, numHarm);
                        
                        % 4. 计算新残差 (不组装J)
                        [R_try, ~, ~] = obj.assembleSystem(space, matLibData, ...
                            X_try_mat, I_try_vec, X_prev_flat, V_curr, dt, ...
                            C_vec_small, M_sys, J_eddy_sys, numDofs, numHarm, harmonics, baseFreq, false);
                        
                        res_try = norm(R_try(free_mask));
                        
                        % 5. 评估: 残差是否下降?
                        if res_try < res_norm
                            % [成功] 接受步长
                            X_curr_flat = X_try_flat;
                            X_curr_mat  = X_try_mat;
                            I_curr_vec  = I_try_vec;
                            
                            % 减小阻尼 (更激进，趋向牛顿法)
                            lambda_old = lambda;
                            lambda = max(lambda * obj.ScaleDown, obj.LambdaMin);
                            
                            % 更新当前点的残差和雅可比 (为下一次迭代做准备)
                            % 注意：这里需要重新计算 J_curr，因为位置变了
                            % fprintf('      [Update] Re-assembling Jacobian at new state...\n');
                            [R_curr, J_curr, ~] = obj.assembleSystem(space, matLibData, ...
                                X_curr_mat, I_curr_vec, X_prev_flat, V_curr, dt, ...
                                C_vec_small, M_sys, J_eddy_sys, numDofs, numHarm, harmonics, baseFreq, true);
                            
                            res_prev = res_norm;
                            res_norm = res_try;
                            step_accepted = true;
                            
                            fprintf('   Iter %d: Res=%.4e RelRes:%.4e | L=%.1e -> %.1e | I_fund=%.2f A\n', ...
                                iter, res_norm, res_norm / norm_base, lambda_old, lambda, abs(I_curr_vec(2)));
                        else
                            % [失败] 拒绝步长
                            % 增大阻尼 (更保守，趋向梯度下降)
                            lambda_old = lambda;
                            lambda = lambda * obj.ScaleUp;
                            
                            % fprintf('      [LM-Rej] Res=%.2e > %.2e. Increasing Damping L: %.1e -> %.1e\n', ...
                            %    res_try, res_norm, lambda_old, lambda);
                            
                            if lambda > obj.LambdaMax
                                warning('LM Damping overflow. Solver stuck.');
                                step_accepted = true; % 强制退出内循环，让外循环结束
                            end
                        end
                    end
                    
                    if lambda > obj.LambdaMax, break; end
                end
                
                X_prev_flat = X_curr_flat;
                Results.I_history(t_idx, :) = I_curr_vec.';
            end
            
            Results.X_final = X_curr_mat;
            info = [];
            fprintf('[Step 3] Solution Complete.\n');
            fprintf('==============================================\n');
        end
        
        function [R_sys, J_sys, row_scale] = assembleSystem(obj, space, matLibData, ...
                X_mat, I_vec, X_prev, V_src, dt, ...
                C_vec, M_sys, J_eddy, nDof, nHarm, harms, freq, calcJ)
            
            if nargin < 16, calcJ = true; end
            X_flat = X_mat(:);
            nField = nDof * nHarm;
            nCirc  = nHarm;
            
            % A. 场方程组装
            [J_tri, R_mag_mat] = obj.Assembler.assembleHBFEM(space, X_mat, obj.AFT, matLibData, calcJ);
            R_field = R_mag_mat(:) + (M_sys/dt)*(X_flat - X_prev) + J_eddy*X_flat;
            
            % 耦合项: -C*I
            F_curr = zeros(nField, 1);
            for k=1:nHarm
                is=(k-1)*nDof+1; ie=k*nDof;
                F_curr(is:ie) = C_vec * I_vec(k);
            end
            R_field = R_field - F_curr;
            
            % B. 电路方程组装
            R_circ = zeros(nHarm, 1);
            for k=1:nHarm
                h = harms(k); w = 2*pi*freq*h;
                Psi = C_vec' * X_mat(:,k);
                Psi_old = C_vec' * X_prev((k-1)*nDof+1:k*nDof);
                
                % V = R*I + jw*Psi + dPsi/dt
                R_circ(k) = obj.CircuitR*I_vec(k) + (Psi-Psi_old)/dt + 1j*w*Psi - V_src(k);
            end
            
            % 缩放 (修复: 确保输出参数 row_scale 被赋值)
            mag_scale = mean(abs(R_field)); if mag_scale<1e-5, mag_scale=1; end
            cir_base = abs(obj.CircuitR); if cir_base<1e-10, cir_base=1; end
            
            row_scale = mag_scale / cir_base; % [Fixed v2.1]
            row_scale = min(max(row_scale, 1e3), 1e9); % 限制范围
            
            R_sys = [R_field; R_circ * row_scale];
            
            % C. 雅可比组装
            if calcJ
                J_field = sparse(J_tri.I, J_tri.J, J_tri.V, nField, nField) + M_sys/dt + J_eddy;
                
                [ci, ~, cv] = find(C_vec); nc = length(cv);
                TR_I=[]; TR_J=[]; TR_V=[];
                BL_I=[]; BL_J=[]; BL_V=[];
                BR_I=(1:nHarm)'; BR_J=(1:nHarm)'; BR_V=ones(nHarm,1)*obj.CircuitR*row_scale;
                
                for k=1:nHarm
                    h=harms(k); w=2*pi*freq*h; offset=(k-1)*nDof;
                    idx=(k-1)*nc+1:k*nc;
                    
                    % Top-Right: -C
                    TR_I(idx)=ci+offset; TR_J(idx)=k; TR_V(idx)=-cv;
                    % Bot-Left: (jw + 1/dt)*C'*scale
                    fact = (1j*w + 1/dt)*row_scale; 
                    BL_I(idx)=k; BL_J(idx)=ci+offset; BL_V(idx)=fact*cv;
                end
                
                J_TR = sparse(TR_I, TR_J, TR_V, nField, nCirc);
                J_BL = sparse(BL_I, BL_J, BL_V, nCirc, nField);
                J_BR = sparse(BR_I, BR_J, BR_V, nCirc, nCirc);
                
                J_sys = [J_field, J_TR; J_BL, J_BR];
            else
                J_sys = [];
            end
        end
    end
end