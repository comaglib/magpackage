classdef TransientCoupledSolver < handle
    % TRANSIENTCOUPLEDSOLVER 场路耦合非线性瞬态求解器 (v1.2 - Dimension Fix)
    %
    % 修复:
    %   1. [Crash Fix] 修复 evaluateCoupledSystem 中 J_fem (N x N) 与 
    %      M_sigma (N+1 x N+1) 直接相加导致的维度不匹配错误。
    %      现在先将 J_fem 扩展为全系统尺寸再进行加法。
    
    properties
        Assembler
        LinearSolver
        
        Dt = 0.01
        
        % --- 非线性迭代配置 ---
        MaxIterations = 50
        Tolerance = 1e-5
        UseLineSearch = true
        MaxLineSearchIters = 10
        Relaxation = 1.0
        
        % --- 正则化配置 ---
        AutoRegularization = true 
        RegScaleRel = 1e-5 
    end
    
    methods
        function obj = TransientCoupledSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_V, matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, fixedDofs_A, fixedDofs_V, init_State)
            
            fprintf('==================================================\n');
            fprintf('   Transient Coupled Solver (v1.2 Fixed)          \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % --- 1. 自由度分配 ---
            [~, offset_A] = dofHandler.distributeDofs(space_A);
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            conductingTags = obj.getConductingRegions(sigmaMap);
            enable_V = ~isempty(conductingTags) && ~isempty(space_V);
            offset_V = -1;
            auto_ground_idx = [];
            
            if enable_V
                if ~dofHandler.DofMaps.isKey(space_V.toString())
                    [~, offset_V] = dofHandler.distributeDofs(space_V, conductingTags);
                else
                    offset_V = dofHandler.SpaceOffsets(space_V.toString());
                end
                
                hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                if ~hasFixedV
                    fprintf('   [Auto-Fix] V-field floating. Grounding 1st node to 0V.\n');
                    auto_ground_idx = offset_V + 1;
                end
            else
                fprintf('   [Info] No eddy current regions (Massive). Pure Coil mode.\n');
            end
            
            numFemDofs = dofHandler.NumGlobalDofs;
            offset_I = numFemDofs;
            numTotalDofs = numFemDofs + 1;
            
            % --- 2. 预计算常数矩阵 ---
            ctx.num_A = num_A;
            ctx.numFemDofs = numFemDofs;
            ctx.numTotalDofs = numTotalDofs;
            ctx.enable_V = enable_V;
            ctx.space_A = space_A;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % M_sigma (扩展到 N_total)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, numTotalDofs, numTotalDofs); 
            
            ctx.M_reg = sparse(numTotalDofs, numTotalDofs);
            if obj.AutoRegularization
                M_reg_local = obj.Assembler.assembleMass(space_A);
                [mr_i, mr_j, mr_v] = find(M_reg_local);
                ctx.M_reg = sparse(mr_i, mr_j, mr_v, numTotalDofs, numTotalDofs);
            end
            
            ctx.C_AV = sparse(numTotalDofs, numTotalDofs);
            ctx.L_VV = sparse(numTotalDofs, numTotalDofs);
            if enable_V
                C_AV_local = obj.Assembler.assembleCoupling(space_A, space_V, sigmaMap);
                L_VV_local = obj.Assembler.assembleScalarLaplacian(space_V, sigmaMap);
                [ci, cj, cv] = find(C_AV_local); ctx.C_AV = sparse(ci, cj, cv, numTotalDofs, numTotalDofs);
                [li, lj, lv] = find(L_VV_local); ctx.L_VV = sparse(li, lj, lv, numTotalDofs, numTotalDofs);
            end
            
            fprintf('   [Init] Assembling winding coupling...\n');
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            % --- 3. 初始化 ---
            x_prev = zeros(numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n_c = min(length(init_State), numTotalDofs);
                x_prev(1:n_c) = init_State(1:n_c);
            end
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed = idx_A;
            if enable_V
                if isempty(fixedDofs_V), idx_V = [];
                elseif islogical(fixedDofs_V), idx_V = find(fixedDofs_V);
                else, idx_V = fixedDofs_V(:); 
                end
                if ~isempty(auto_ground_idx), idx_V = unique([idx_V; auto_ground_idx]); end
                idx_V = idx_V(idx_V <= numFemDofs);
                all_fixed = [all_fixed; idx_V];
            end
            
            isLinearSystem = obj.checkLinearity(matLib);
            ctx.RegScaleAbs = 0; 
            
            numSteps = length(timeSteps);
            solutionResults = cell(numSteps, 1);
            currentTime = 0;
            
            % --- 4. 时间步进 ---
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f) ---\n', t_idx, numSteps, dt, currentTime);
                
                ctx.dt = dt;
                ctx.x_prev = x_prev;
                ctx.V_source_val = circuitProps.V_source_func(currentTime);
                
                x_curr = x_prev;
                converged = false;
                
                % 初次评估
                [J_sys, R_sys] = obj.evaluateCoupledSystem(x_curr, ctx, true);
                
                % [Auto-Reg]
                if ctx.RegScaleAbs == 0 && obj.AutoRegularization
                    % 只统计 FEM 部分
                    J_fem_diag = abs(diag(J_sys(1:numFemDofs, 1:numFemDofs)));
                    avg_diag = mean(J_fem_diag(J_fem_diag > 0));
                    if isnan(avg_diag) || avg_diag==0, avg_diag = 1.0; end
                    ctx.RegScaleAbs = full(avg_diag * obj.RegScaleRel);
                    
                    fprintf('   [Auto-Reg] Scale: %.2e\n', ctx.RegScaleAbs);
                    J_sys = J_sys + ctx.RegScaleAbs * ctx.M_reg;
                    R_sys = R_sys + (ctx.RegScaleAbs * ctx.M_reg) * x_curr;
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                norm_res = norm(Res_bc);
                
                for iter = 1:obj.MaxIterations
                    if norm_res < obj.Tolerance
                        converged = true;
                        fprintf('      Iter %d: Converged (Res=%.4e)\n', iter, norm_res);
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
                            
                            [~, R_try] = obj.evaluateCoupledSystem(x_try, ctx, false);
                            R_try(all_fixed) = 0; 
                            if obj.AutoRegularization
                                R_try = R_try + (ctx.RegScaleAbs * ctx.M_reg) * x_try;
                            end
                            
                            if norm(R_try) < res_prev || (iter==1 && k<3)
                                alpha = alpha_try;
                                if k > 0, fprintf('      [LS] Step %.4f\n', alpha); end
                                break;
                            end
                        end
                    end
                    
                    x_curr = x_curr + alpha * delta_x;
                    
                    [J_sys, R_sys] = obj.evaluateCoupledSystem(x_curr, ctx, true);
                    if obj.AutoRegularization
                        J_sys = J_sys + ctx.RegScaleAbs * ctx.M_reg;
                        R_sys = R_sys + (ctx.RegScaleAbs * ctx.M_reg) * x_curr;
                    end
                    
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed);
                    norm_res = norm(Res_bc);
                    
                    I_val = x_curr(end);
                    fprintf('      Iter %d: Res=%.4e, I=%.4f A\n', iter, norm_res, I_val);
                    
                    if isLinearSystem && norm_res < 1e-9, converged = true; break; end
                end
                
                if ~converged, fprintf('   [Warn] Not converged.\n'); end
                x_prev = x_curr;
                solutionResults{t_idx} = x_curr;
            end
            
            info.FinalTime = currentTime;
        end
        
        function tags = getConductingRegions(~, sigmaMap)
            tags = [];
            if isa(sigmaMap, 'containers.Map')
                keys = sigmaMap.keys;
                for i = 1:length(keys)
                    k = keys{i}; val = sigmaMap(k);
                    if max(abs(val)) > 1e-12
                        if ischar(k), k = str2double(k); end
                        tags = [tags, k]; 
                    end
                end
            elseif isnumeric(sigmaMap)
                tags = find(sigmaMap > 1e-12);
            end
        end
    end
    
    methods (Access = private)
        function [J_out, R_out] = evaluateCoupledSystem(obj, x, ctx, calc_J)
            A_vec = x(1:ctx.num_A);
            I_val = x(end);
            
            % 1. 物理 Jacobian (Curl-Curl)
            [J_fem_local, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            
            % [CRITICAL FIX] 尺寸对齐
            % J_fem_local 是 (N_fem x N_fem)，必须扩展为 (N_total x N_total)
            % R_nu 是 (N_fem x 1)，必须扩展为 (N_total x 1)
            
            if size(R_nu, 1) < ctx.numTotalDofs
                R_nu(ctx.numTotalDofs, 1) = 0; 
            end
            
            % 2. 瞬态项 (dA/dt)
            dx_dt = (x - ctx.x_prev) / ctx.dt;
            dA_dt_term = ctx.M_sigma * dx_dt;
            
            R_out = R_nu + dA_dt_term;
            
            % Row A: - W * I
            R_out = R_out - ctx.W_vec * I_val;
            
            if ctx.enable_V
                R_out = R_out + ctx.C_AV * x;
                R_V_rows = (ctx.C_AV' * dx_dt) + ctx.L_VV * x;
                R_out = R_out + R_V_rows; 
            end
            
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt; 
            Res_I = ctx.V_source_val - term_R - term_L - term_EMF;
            
            R_out(end) = -Res_I; 
            
            J_out = [];
            if calc_J
                % [FIX] 先将 J_fem_local 扩展为 N_total x N_total
                [ji, jj, jv] = find(J_fem_local);
                J_fem_ext = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
                
                % Block FEM
                J_out = J_fem_ext + ctx.M_sigma / ctx.dt;
                
                if ctx.enable_V
                    J_out = J_out + ctx.C_AV + (ctx.C_AV') / ctx.dt + ctx.L_VV;
                end
                
                % Coupling Blocks
                [w_idx, ~, w_val] = find(ctx.W_vec);
                num_nz = length(w_idx);
                
                % Row A, Col I: -W
                J_col_I = sparse(w_idx, repmat(ctx.numTotalDofs, num_nz, 1), w_val, ctx.numTotalDofs, ctx.numTotalDofs);
                J_out = J_out - J_col_I;
                
                % Row I, Col A: +W'/dt
                J_row_I = sparse(repmat(ctx.numTotalDofs, num_nz, 1), w_idx, w_val / ctx.dt, ctx.numTotalDofs, ctx.numTotalDofs);
                J_out = J_out + J_row_I;
                
                % Block Circuit
                val_II = ctx.circuit.R + ctx.circuit.L / ctx.dt;
                J_out(end, end) = J_out(end, end) + val_II;
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