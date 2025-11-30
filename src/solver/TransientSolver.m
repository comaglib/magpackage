classdef TransientSolver < handle
    % TRANSIENTSOLVER 瞬态非线性磁场求解器 (v3.1 - Fix Sparse Print)
    %
    % 修复:
    %   1. [Fix] fprintf 报错: 强制将 RegScaleAbs 转换为 full 类型。
    %   2. 保持 v3.0 的所有特性 (Auto-Grounding, Adaptive Reg, NaN Check)。
    
    properties
        Assembler
        LinearSolver
        
        Dt = 0.01
        
        % 非线性迭代配置
        MaxIterations = 50
        Tolerance = 1e-5 
        UseLineSearch = true
        MaxLineSearchIters = 10
        Relaxation = 1.0
        
        % 正则化配置
        AutoRegularization = true 
        RegScaleRel = 1e-5 % 相对比例系数
    end
    
    methods
        function obj = TransientSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_V, matLib, sigmaMap, sourceFunc, timeSteps, fixedDofs_A, fixedDofs_V, init_A, init_V)
            
            fprintf('==============================================\n');
            fprintf('   Transient Solver (Robust A-V Formulation)  \n');
            fprintf('==============================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % --- 1. 自由度分配 ---
            [~, offset_A] = dofHandler.distributeDofs(space_A); 
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            conductingTags = obj.getConductingRegions(sigmaMap);
            enable_V = ~isempty(conductingTags) && ~isempty(space_V);
            offset_V = -1;
            
            if enable_V
                if ~dofHandler.DofMaps.isKey(space_V.toString())
                    [~, offset_V] = dofHandler.distributeDofs(space_V, conductingTags);
                else
                    offset_V = dofHandler.SpaceOffsets(space_V.toString());
                end
                
                % [Auto-Grounding] 检查悬浮电位
                hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                if ~hasFixedV
                    fprintf('   [Auto-Fix] V-field is floating (no BCs). Grounding 1st node to 0V.\n');
                    v_dof_start = offset_V + 1;
                    if isempty(fixedDofs_V)
                        fixedDofs_V = []; 
                    end
                    auto_ground_idx = v_dof_start;
                else
                    auto_ground_idx = [];
                end
            else
                fprintf('   [Info] Pure A-field solve.\n');
                auto_ground_idx = [];
            end
            
            numTotalDofs = dofHandler.NumGlobalDofs;
            
            % --- 2. 预计算常数矩阵 ---
            sysCtx.num_A = num_A;
            sysCtx.numTotalDofs = numTotalDofs;
            sysCtx.enable_V = enable_V;
            sysCtx.space_A = space_A;
            sysCtx.matLib = matLib;
            
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            sysCtx.M_sigma = sparse(mi, mj, mv, numTotalDofs, numTotalDofs);
            
            sysCtx.M_reg = sparse(numTotalDofs, numTotalDofs);
            if obj.AutoRegularization
                M_reg_local = obj.Assembler.assembleMass(space_A);
                [mr_i, mr_j, mr_v] = find(M_reg_local);
                sysCtx.M_reg = sparse(mr_i, mr_j, mr_v, numTotalDofs, numTotalDofs);
            end
            
            sysCtx.C_AV = sparse(numTotalDofs, numTotalDofs);
            sysCtx.L_VV = sparse(numTotalDofs, numTotalDofs);
            if enable_V
                sysCtx.C_AV = obj.Assembler.assembleCoupling(space_A, space_V, sigmaMap);
                sysCtx.L_VV = obj.Assembler.assembleScalarLaplacian(space_V, sigmaMap);
            end
            
            % --- 3. 初始化 ---
            x_prev = zeros(numTotalDofs, 1);
            if nargin >= 10 && ~isempty(init_A), x_prev(1:min(length(init_A), num_A)) = init_A; end
            if nargin >= 11 && ~isempty(init_V) && enable_V, x_prev(offset_V+1:end) = init_V; end
            
            % 合并边界条件
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed_dofs = idx_A;
            
            if enable_V
                if isempty(fixedDofs_V)
                    idx_V = [];
                elseif islogical(fixedDofs_V)
                    idx_V = find(fixedDofs_V);
                else
                    idx_V = fixedDofs_V(:);
                end
                
                if ~isempty(auto_ground_idx)
                    idx_V = unique([idx_V; auto_ground_idx]);
                end
                idx_V = idx_V(idx_V <= numTotalDofs); 
                all_fixed_dofs = [all_fixed_dofs; idx_V];
            end
            
            isLinearSystem = obj.checkLinearity(matLib);
            
            % 初始化正则化系数
            sysCtx.RegScaleAbs = 0; 
            
            numSteps = length(timeSteps);
            solutionResults = cell(numSteps, 1);
            currentTime = 0;
            
            % --- 4. 时间步进 ---
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f) ---\n', t_idx, numSteps, dt, currentTime);
                
                sysCtx.dt = dt;
                sysCtx.x_prev = x_prev;
                
                currentSourceMap = sourceFunc(currentTime);
                F_A = obj.Assembler.assembleSource(space_A, currentSourceMap);
                sysCtx.F_sys = sparse(numTotalDofs, 1);
                if ~isempty(F_A), sysCtx.F_sys(1:length(F_A)) = F_A; end
                
                x_curr = x_prev;
                converged = false;
                
                [J_sys, R_sys] = obj.evaluateSystemState(x_curr, sysCtx, true);
                
                % [自适应正则化]
                if sysCtx.RegScaleAbs == 0 && obj.AutoRegularization
                    diag_J = abs(diag(J_sys));
                    avg_diag = mean(diag_J(diag_J > 0));
                    % [CRITICAL FIX] 转换为 full 类型，防止 fprintf 报错
                    if isnan(avg_diag) || avg_diag == 0, avg_diag = 1.0; end
                    sysCtx.RegScaleAbs = full(avg_diag * obj.RegScaleRel);
                    
                    fprintf('   [Auto-Reg] Absolute Scale set to %.2e (Rel: %.1e)\n', sysCtx.RegScaleAbs, obj.RegScaleRel);
                    
                    % 重新添加正则化
                    J_sys = J_sys + sysCtx.RegScaleAbs * sysCtx.M_reg;
                    R_sys = R_sys + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x_curr;
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed_dofs);
                norm_res = norm(Res_bc);
                
                for iter = 1:obj.MaxIterations
                    if norm_res < obj.Tolerance
                        converged = true;
                        fprintf('      Iter %d: Converged (Res=%.4e)\n', iter, norm_res);
                        break;
                    end
                    
                    delta_x = obj.LinearSolver.solve(J_bc, -Res_bc);
                    
                    if isempty(delta_x) || any(isnan(delta_x)), 
                        warning('Linear solver returned NaN/Empty.'); 
                        break; 
                    end
                    
                    % 线搜索
                    alpha = obj.Relaxation;
                    if obj.UseLineSearch && ~isLinearSystem
                        res_prev = norm_res;
                        for k = 0:obj.MaxLineSearchIters
                            alpha_try = alpha * (0.5^k);
                            x_try = x_curr + alpha_try * delta_x;
                            
                            [~, R_try] = obj.evaluateSystemState(x_try, sysCtx, false);
                            R_try(all_fixed_dofs) = 0; 
                            res_try = norm(R_try);
                            
                            if res_try < res_prev || (iter==1 && k<3)
                                alpha = alpha_try;
                                if k > 0, fprintf('      [LS] Step %.4f, Res %.2e\n', alpha, res_try); end
                                break;
                            end
                            
                            if res_try > res_prev * 10
                                alpha = alpha * 0.1; 
                            end
                        end
                    end
                    
                    x_curr = x_curr + alpha * delta_x;
                    
                    [J_sys, R_sys] = obj.evaluateSystemState(x_curr, sysCtx, true);
                    
                    if obj.AutoRegularization
                        J_sys = J_sys + sysCtx.RegScaleAbs * sysCtx.M_reg;
                        R_sys = R_sys + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x_curr;
                    end
                    
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed_dofs);
                    norm_res = norm(Res_bc);
                    
                    fprintf('      Iter %d: Res=%.4e\n', iter, norm_res);
                    
                    if norm_res > 1e20
                        error('Divergence detected (Res > 1e20). Check BCs or Time Step.');
                    end
                    
                    if isLinearSystem && norm_res < 1e-9, converged = true; break; end
                end
                
                if ~converged, fprintf('   [Warn] Not converged at step %d.\n', t_idx); end
                x_prev = x_curr;
                solutionResults{t_idx} = x_curr;
            end
            
            info.FinalTime = currentTime;
        end
    end
    
    methods (Access = private)
        function [J_out, R_out] = evaluateSystemState(obj, x_eval, ctx, calc_J)
            A_vec = x_eval(1:ctx.num_A);
            
            % 1. 物理 Jacobian
            [J_nu, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            if size(R_nu, 1) < ctx.numTotalDofs, R_nu(ctx.numTotalDofs, 1) = 0; end
            
            % 2. 瞬态项
            dA_dt_term = ctx.M_sigma * (x_eval - ctx.x_prev) / ctx.dt;
            
            % 3. 残差
            R_out = R_nu + dA_dt_term - ctx.F_sys;
            
            % 4. 耦合
            if ctx.enable_V
                R_out = R_out + ctx.C_AV * x_eval;
                R_V_rows = (ctx.C_AV' * (x_eval - ctx.x_prev)) / ctx.dt + ctx.L_VV * x_eval;
                R_out = R_out + R_V_rows;
            end
            
            J_out = [];
            if calc_J
                J_out = J_nu + ctx.M_sigma / ctx.dt;
                if ctx.enable_V
                    J_out = J_out + ctx.C_AV + (ctx.C_AV') / ctx.dt + ctx.L_VV;
                end
            end
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
        
        function isLin = checkLinearity(~, matLib)
            isLin = true;
            if isa(matLib, 'containers.Map')
               k = matLib.keys;
               for i=1:length(k)
                   m = matLib(k{i});
                   if strcmpi(m.Type, 'Nonlinear'), isLin = false; return; end
               end
            else
                for i=1:length(matLib)
                    if strcmpi(matLib(i).Type, 'Nonlinear'), isLin = false; return; end
                end
            end
        end
    end
end