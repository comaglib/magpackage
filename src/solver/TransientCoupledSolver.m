classdef TransientCoupledSolver < handle
    % TRANSIENTCOUPLEDSOLVER 场路耦合非线性瞬态求解器 (v4.0 - Regularization)
    % 
    % 描述:
    %   回归正则化方法 (Penalty/Regularization)，放弃拉格朗日乘子。
    %   通过在刚度矩阵中添加微小的质量矩阵项 (epsilon * M) 来消除 A 的零空间奇异性。
    %   这种方法在非线性高反差问题中比 Lagrange 乘子法更鲁棒。
    %
    % 更新日志:
    %   v4.0: [Revert] 移除 space_P (Lagrange Multiplier)。
    %         [Add] 引入 Mass-Matrix Regularization。
    
    properties
        Assembler
        LinearSolver
        
        Dt = 0.01
        
        MaxIterations = 50
        Tolerance = 1e-3
        RelTolerance = 1e-2
        UseLineSearch = true
        MaxLineSearchIters = 10
        Relaxation = 1.0
        
        RegularizationFactor = 8e-2; 
        CircuitScale = 1.0;
    end
    
    methods
        function obj = TransientCoupledSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 400; 
            obj.LinearSolver.MumpsSymmetry = 0; % 非对称模式 (为了通用性)
        end
        
        function [solutionResults, info] = solve(obj, space_A, matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, fixedDofs_A, init_State)
            
            fprintf('==================================================\n');
            fprintf('   Transient Coupled Solver (v4.0 Regularization) \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % 1. 自由度分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            numFemDofs = dofHandler.NumGlobalDofs;
            numTotalDofs = numFemDofs + 1; % +1 Circuit Dof
            
            % 2. 预计算
            ctx.num_A = num_A; 
            ctx.numTotalDofs = numTotalDofs;
            ctx.space_A = space_A; ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % M_sigma: 导电区域的涡流项 (sigma * dA/dt)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, numTotalDofs, numTotalDofs); 
            
            % M_reg: 正则化质量矩阵 (全域 Mass)
            % 用于填充 A 的零空间 (Null Space)，防止奇异
            fprintf('   [Reg] Assembling regularization mass matrix...\n');
            M_reg_local = obj.Assembler.assembleMass(space_A);
            [ri, rj, rv] = find(M_reg_local);
            ctx.M_reg = sparse(ri, rj, rv, numTotalDofs, numTotalDofs);
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            ctx.CircuitRowScale = obj.CircuitScale;
            
            % 3. 初始化与BC
            x_prev = zeros(numTotalDofs, 1);
            if nargin >= 9 && ~isempty(init_State)
                n_c = min(length(init_State), numTotalDofs);
                x_prev(1:n_c) = init_State(1:n_c);
            end
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed = idx_A; % 只有 A 的边界，没有 P
            all_fixed = all_fixed(all_fixed <= numTotalDofs);
            
            isLinearSystem = obj.checkLinearity(matLib);
            numSteps = length(timeSteps);
            solutionResults = []; currentHistory = zeros(numSteps, 1); currentTime = 0;
            
            % 4. 求解循环
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f) ---\n', t_idx, numSteps, dt, currentTime);
                ctx.dt = dt; ctx.x_prev = x_prev;
                ctx.V_source_val = circuitProps.V_source_func(currentTime);
                
                x_curr = x_prev; converged = false;
                
                % 首次评估 (计算初始残差)
                [J_sys, R_sys, eps_scale] = obj.evaluateSystem(x_curr, ctx, true, 0);
                
                % 确定正则化强度 (基于刚度矩阵对角线)
                if t_idx == 1
                     avg_diag = mean(abs(diag(J_sys(1:num_A, 1:num_A))));
                     if avg_diag == 0, avg_diag = 1.0; end
                     ctx.eps_val = avg_diag * obj.RegularizationFactor;
                     fprintf('   [Reg] Regularization Strength: %.2e (Factor=%.1e)\n', full(ctx.eps_val), obj.RegularizationFactor);
                     
                     % 重新评估带正则化的系统
                     [J_sys, R_sys] = obj.evaluateSystem(x_curr, ctx, true, ctx.eps_val);
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                norm_res = norm(Res_bc);
                initial_res = norm_res; if initial_res < 1e-20, initial_res = 1.0; end 
                
                for iter = 1:obj.MaxIterations
                    rel_res = norm_res / initial_res;
                    
                    if norm_res < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        converged = true;
                        fprintf('      Iter %d: Converged (Res=%.4e, Rel=%.2e)\n', iter, norm_res, rel_res);
                        break;
                    end
                    
                    delta_x = obj.LinearSolver.solve(J_bc, -Res_bc);
                    
                    if isempty(delta_x) || any(isnan(delta_x)), warning('Solver failed.'); break; end
                    
                    alpha = obj.Relaxation;
                    if obj.UseLineSearch && ~isLinearSystem
                        res_prev = norm_res;
                        for k = 0:obj.MaxLineSearchIters
                            alpha_try = alpha * (0.5^k);
                            x_try = x_curr + alpha_try * delta_x;
                            [~, R_try] = obj.evaluateSystem(x_try, ctx, false, ctx.eps_val);
                            R_try(all_fixed) = 0; 
                            if norm(R_try) < res_prev || (iter==1 && k<3)
                                alpha = alpha_try;
                                if k > 0, fprintf('      [LS] Step %.4f\n', alpha); end
                                break;
                            end
                        end
                    end
                    
                    x_curr = x_curr + alpha * delta_x;
                    
                    [J_sys, R_sys] = obj.evaluateSystem(x_curr, ctx, true, ctx.eps_val);
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                    norm_res = norm(Res_bc);
                    
                    I_val = x_curr(end);
                    fprintf('    Iter %d: Res = %.4e, RelRes = %.4e, I = %.6e\n', ...
                        iter, norm_res, norm_res/initial_res, I_val);
                    
                    if isLinearSystem && norm_res < 1e-9, converged = true; break; end
                end
                
                if ~converged, fprintf('   [Warn] Not converged.\n'); end
                x_prev = x_curr;
                solutionResults = x_curr;
                currentHistory(t_idx) = x_curr(end);
            end
            
            info.FinalTime = currentTime; info.CurrentHistory = currentHistory;
        end
    end
    
    methods (Access = private)
        function [J_out, R_out, eps_scale] = evaluateSystem(obj, x, ctx, calc_J, eps_reg)
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            % 1. 物理刚度 (Curl-Curl)
            [J_fem_local, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            if size(R_nu, 1) < ctx.numTotalDofs, R_nu(ctx.numTotalDofs, 1) = 0; end
            
            dx_dt = (x - ctx.x_prev) / ctx.dt;
            
            % 2. 涡流项 + 正则化项
            % R_total = R_nu + M_sigma * dx/dt + eps * M_reg * A
            % 注意: 正则化是对 A 本身进行惩罚 (类似 A term)，不是对 dx/dt
            term_eddy = ctx.M_sigma * dx_dt; 
            term_reg  = eps_reg * ctx.M_reg * x; % 惩罚项: eps * M * A
            
            R_out = R_nu + term_eddy + term_reg;
            
            % 3. 电路耦合 (RHS)
            R_out = R_out - ctx.W_vec * I_val;
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt; 
            Res_I = ctx.V_source_val - term_R - term_L - term_EMF;
            R_out(end) = -Res_I * ctx.CircuitRowScale;
            
            J_out = [];
            eps_scale = 0;
            if calc_J
                [ji, jj, jv] = find(J_fem_local);
                J_fem_ext = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
                
                % J = K(nu) + M_sigma/dt + eps * M_reg
                J_out = J_fem_ext + ctx.M_sigma / ctx.dt + eps_reg * ctx.M_reg;
                
                % 电路耦合项
                [w_idx, ~, w_val] = find(ctx.W_vec); num_nz = length(w_idx);
                J_col_I = sparse(w_idx, repmat(ctx.numTotalDofs, num_nz, 1), w_val, ctx.numTotalDofs, ctx.numTotalDofs);
                J_out = J_out - J_col_I;
                
                J_row_I = sparse(repmat(ctx.numTotalDofs, num_nz, 1), w_idx, w_val / ctx.dt, ctx.numTotalDofs, ctx.numTotalDofs);
                J_out = J_out + J_row_I * ctx.CircuitRowScale;
                
                val_II = ctx.circuit.R + ctx.circuit.L / ctx.dt;
                J_out(end, end) = J_out(end, end) + val_II * ctx.CircuitRowScale;
            end
        end
        
        function isLin = checkLinearity(~, matLib)
            isLin = true;
            if isa(matLib, 'containers.Map')
               k = matLib.keys;
               for i=1:length(k), m = matLib(k{i}); if strcmpi(m.Type, 'Nonlinear'), isLin=false; return; end; end
            end
        end
    end
end