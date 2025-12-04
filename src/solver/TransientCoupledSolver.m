classdef TransientCoupledSolver < handle
    % TRANSIENTCOUPLEDSOLVER 场路耦合非线性瞬态求解器 (v3.5 - Stiffness Matching Fix)
    %
    % 更新日志:
    %   v3.5: [Physics Fix] 彻底修复 Locking 问题。
    %         将 LagrangeScale 的计算基准从"空气刚度"改为"铁芯刚度"。
    %         大幅降低约束项的数值强度，防止其掩盖铁芯内的物理磁化过程。
    %   v3.3: [Fix] 修复矩阵转置块缺失。
    
    properties
        Assembler
        LinearSolver
        
        Dt = 0.01
        
        % --- 非线性迭代配置 ---
        MaxIterations = 50
        Tolerance = 1e-4
        RelTolerance = 1e-4    
        UseLineSearch = true
        MaxLineSearchIters = 10
        Relaxation = 1.0
        
        % [Scaling]
        LagrangeScale = 1.0; 
        CircuitScale = 1e5;
    end
    
    methods
        function obj = TransientCoupledSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 400; 
            obj.LinearSolver.MumpsSymmetry = 0; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, fixedDofs_A, fixedDofs_P, init_State)
            
            fprintf('==================================================\n');
            fprintf('   Transient Coupled Solver (v3.5 Physics Unlocked)\n');
            fprintf('==================================================\n');
            
            if isempty(space_P), error('Need space_P for Lagrange method.'); end
            dofHandler = obj.Assembler.DofHandler;
            
            % 1. 自由度分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            [~, offset_A] = dofHandler.distributeDofs(space_A);
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            if ~dofHandler.DofMaps.isKey(space_P.toString())
                dofHandler.distributeDofs(space_P, unique(obj.Assembler.Mesh.RegionTags));
            end
            offset_P = dofHandler.SpaceOffsets(space_P.toString());
            num_P = dofHandler.SpaceLocalSizes(space_P.toString());
            
            numFemDofs = dofHandler.NumGlobalDofs;
            numTotalDofs = numFemDofs + 1;
            
            % 2. 预计算与 Scaling
            ctx.num_A = num_A; ctx.num_P = num_P; ctx.offset_P = offset_P;
            ctx.numTotalDofs = numTotalDofs;
            ctx.space_A = space_A; ctx.matLib = matLib; ctx.circuit = circuitProps;
            
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, numTotalDofs, numTotalDofs); 
            
            % G 矩阵
            unitMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
            tags = unique(obj.Assembler.Mesh.RegionTags);
            for i=1:length(tags), unitMap(tags(i))=1.0; end
            G_local = obj.Assembler.assembleCoupling(space_A, space_P, unitMap);
            
            % [Auto-Scaling Correction]
            % 关键修改：使用 1e3 (典型铁芯磁阻率) 而不是 1e7 (真空) 作为刚度基准
            % 这样约束项的量级将与铁芯物理项匹配，避免"锁死"铁芯。
            target_nu = 1e3;
            temp_K = obj.Assembler.assembleStiffness(space_A, target_nu); 
            
            avg_K = mean(abs(nonzeros(temp_K)));
            avg_G = mean(abs(nonzeros(G_local)));
            
            % 稍微给一点安全余量 (*10)，保证约束依然有效，但不至于压倒物理
            obj.LagrangeScale = (avg_K / avg_G) * 10.0;
            
            fprintf('   [Scaling] Target_Nu=%.0f -> Scale=%.1e (Iron-Matched)\n', ...
                target_nu, obj.LagrangeScale);
            
            [gi, gj, gv] = find(G_local);
            ctx.G_mat = sparse(gi, gj, gv * obj.LagrangeScale, numTotalDofs, numTotalDofs);
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            ctx.CircuitRowScale = obj.CircuitScale;
            
            % 3. 初始化与BC
            x_prev = zeros(numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n_c = min(length(init_State), numTotalDofs);
                x_prev(1:n_c) = init_State(1:n_c);
            end
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            
            idx_P_global = idx_P + offset_P;
            all_fixed = [idx_A; idx_P_global];
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
                
                [J_sys, R_sys] = obj.evaluateSystem(x_curr, ctx, true);
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                norm_res = norm(Res_bc);
                
                initial_res = norm_res;
                if initial_res < 1e-20, initial_res = 1.0; end 
                
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
                            [~, R_try] = obj.evaluateSystem(x_try, ctx, false);
                            R_try(all_fixed) = 0; 
                            if norm(R_try) < res_prev || (iter==1 && k<3)
                                alpha = alpha_try;
                                if k > 0, fprintf('      [LS] Step %.4f\n', alpha); end
                                break;
                            end
                        end
                    end
                    
                    x_curr = x_curr + alpha * delta_x;
                    
                    [J_sys, R_sys] = obj.evaluateSystem(x_curr, ctx, true);
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                    norm_res = norm(Res_bc);
                    
                    I_val = x_curr(end);
                    fprintf('      Iter %d: Res=%.4e, I=%.4f A\n', iter, norm_res, I_val);
                    
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
        function [J_out, R_out] = evaluateSystem(obj, x, ctx, calc_J)
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            [J_fem_local, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            if size(R_nu, 1) < ctx.numTotalDofs, R_nu(ctx.numTotalDofs, 1) = 0; end
            
            dx_dt = (x - ctx.x_prev) / ctx.dt;
            
            term_eddy = ctx.M_sigma * dx_dt; 
            term_constraint = ctx.G_mat * x; 
            term_constraint_T = ctx.G_mat.' * x;
            
            R_out = R_nu + term_eddy + term_constraint + term_constraint_T;
            R_out = R_out - ctx.W_vec * I_val;
            
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt; 
            Res_I = ctx.V_source_val - term_R - term_L - term_EMF;
            R_out(end) = -Res_I * ctx.CircuitRowScale;
            
            J_out = [];
            if calc_J
                [ji, jj, jv] = find(J_fem_local);
                J_fem_ext = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
                
                J_out = J_fem_ext + ctx.M_sigma / ctx.dt;
                J_out = J_out + ctx.G_mat + ctx.G_mat.';
                
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