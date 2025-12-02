classdef TransientSolver < handle
    % TRANSIENTSOLVER 瞬态非线性磁场求解器 (v3.2 - Robust & Optimized)
    %
    % 更新日志:
    %   1. [Robust] 引入相对容差 (RelTolerance)，解决大数系统(如高电导率)的收敛停滞问题。
    %   2. [Perf] 优化 Jacobian 组装，将线性恒定项移出牛顿循环。
    %   3. [Log] 增强收敛信息的输出。
    
    properties
        Assembler
        LinearSolver
        
        Dt = 0.01
        
        % 非线性迭代配置
        MaxIterations = 50
        Tolerance = 1e-5       % 绝对容差 (适用于场值较小的情况)
        RelTolerance = 1e-6    % [新增] 相对容差 (关键! 适用于刚度矩阵数值很大的情况)
        UseLineSearch = true
        MaxLineSearchIters = 10
        Relaxation = 1.0
        
        % 正则化配置
        AutoRegularization = true 
        RegScaleRel = 1e-5 
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
                
                % [Auto-Grounding]
                hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                if ~hasFixedV
                    fprintf('   [Auto-Fix] V-field is floating (no BCs). Grounding 1st node to 0V.\n');
                    v_dof_start = offset_V + 1;
                    if isempty(fixedDofs_V), fixedDofs_V = []; end
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
            
            % M_sigma (Eddy term)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            sysCtx.M_sigma = sparse(mi, mj, mv, numTotalDofs, numTotalDofs);
            
            % Regularization Matrix
            sysCtx.M_reg = sparse(numTotalDofs, numTotalDofs);
            if obj.AutoRegularization
                M_reg_local = obj.Assembler.assembleMass(space_A);
                [mr_i, mr_j, mr_v] = find(M_reg_local);
                sysCtx.M_reg = sparse(mr_i, mr_j, mr_v, numTotalDofs, numTotalDofs);
            end
            
            % Coupling Matrices
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
            
            % 合并边界条件索引
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
                
                % 组装源向量
                currentSourceMap = sourceFunc(currentTime);
                F_A = obj.Assembler.assembleSource(space_A, currentSourceMap);
                sysCtx.F_sys = sparse(numTotalDofs, 1);
                if ~isempty(F_A), sysCtx.F_sys(1:length(F_A)) = F_A; end
                
                % [Optimization] 预计算当前时间步的恒定 Jacobian 部分 (M/dt + C + C'/dt + L)
                J_const = sysCtx.M_sigma / dt;
                if enable_V
                    J_const = J_const + sysCtx.C_AV + (sysCtx.C_AV') / dt + sysCtx.L_VV;
                end
                
                x_curr = x_prev;
                converged = false;
                initial_res = 1.0; % 用于相对收敛判断
                
                % --- 牛顿迭代 ---
                for iter = 1:obj.MaxIterations
                    
                    % 1. 计算状态 (J_nu: 仅刚度部分, R_sys: 完整残差)
                    [J_nu, R_sys] = obj.evaluateSystemState(x_curr, sysCtx, true);
                    
                    % 2. 组装完整 Jacobian
                    J_sys = J_nu + J_const;
                    
                    % 3. 自适应正则化 (基于第一次迭代的对角线)
                    if sysCtx.RegScaleAbs == 0 && obj.AutoRegularization
                        diag_J = abs(diag(J_sys));
                        avg_diag = mean(diag_J(diag_J > 0));
                        if isnan(avg_diag) || avg_diag == 0, avg_diag = 1.0; end
                        sysCtx.RegScaleAbs = full(avg_diag * obj.RegScaleRel);
                        fprintf('   [Auto-Reg] Absolute Scale set to %.2e (Rel: %.1e)\n', sysCtx.RegScaleAbs, obj.RegScaleRel);
                    end
                    
                    % 应用正则化
                    if obj.AutoRegularization && sysCtx.RegScaleAbs > 0
                        reg_term = sysCtx.RegScaleAbs * sysCtx.M_reg;
                        J_sys = J_sys + reg_term;
                        R_sys = R_sys + reg_term * x_curr;
                    end
                    
                    % 4. 边界条件处理
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed_dofs);
                    norm_res = norm(Res_bc);
                    
                    % 5. 收敛检查 (绝对 + 相对)
                    if iter == 1
                        initial_res = norm_res;
                        if initial_res < 1e-20, initial_res = 1.0; end % 防止除零
                    end
                    
                    rel_res = norm_res / initial_res;
                    
                    % 打印状态
                    if iter == 1
                        fprintf('      Iter 1: Res=%.4e (Init)\n', norm_res);
                    else
                        fprintf('      Iter %d: Res=%.4e (Rel=%.2e)\n', iter, norm_res, rel_res);
                    end

                    if norm_res < obj.Tolerance || rel_res < obj.RelTolerance
                        converged = true;
                        fprintf('      -> Converged.\n');
                        break;
                    end
                    
                    % 6. 线性求解
                    delta_x = obj.LinearSolver.solve(J_bc, -Res_bc);
                    
                    if isempty(delta_x) || any(isnan(delta_x))
                        warning('Linear solver returned NaN/Empty.'); 
                        break; 
                    end
                    
                    % 7. 线搜索 (Line Search)
                    alpha = obj.Relaxation;
                    if obj.UseLineSearch && ~isLinearSystem
                        res_prev_val = norm_res;
                        for k = 0:obj.MaxLineSearchIters
                            alpha_try = alpha * (0.5^k);
                            x_try = x_curr + alpha_try * delta_x;
                            
                            [~, R_try] = obj.evaluateSystemState(x_try, sysCtx, false);
                            
                            % 应用正则化残差 (仅用于评估)
                            if obj.AutoRegularization && sysCtx.RegScaleAbs > 0
                                R_try = R_try + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x_try;
                            end
                            
                            % Mask BCs
                            R_try(all_fixed_dofs) = 0; 
                            res_try = norm(R_try);
                            
                            if res_try < res_prev_val || (iter==1 && k<3)
                                alpha = alpha_try;
                                if k > 0, fprintf('      [LS] Step %.4f, Res %.2e\n', alpha, res_try); end
                                break;
                            end
                            
                            if res_try > res_prev_val * 10
                                alpha = alpha * 0.1; 
                            end
                        end
                    end
                    
                    x_curr = x_curr + alpha * delta_x;
                    
                    if norm_res > 1e20
                        error('Divergence detected (Res > 1e20). Check BCs or Time Step.');
                    end
                    
                    if isLinearSystem && norm_res < 1e-9, converged = true; break; end
                end
                
                if ~converged
                    fprintf('   [Warn] Not converged at step %d (Final Res=%.4e).\n', t_idx, norm_res); 
                end
                
                x_prev = x_curr;
                solutionResults{t_idx} = x_curr;
            end
            
            info.FinalTime = currentTime;
        end
    end
    
    methods (Access = private)
        function [J_nu, R_out] = evaluateSystemState(obj, x_eval, ctx, calc_J)
            % EVALUATESYSTEMSTATE 计算系统的非线性残差和切线刚度
            % 
            % Output:
            %   J_nu: 仅包含非线性磁阻率部分的 Jacobian (Stiffness Matrix)
            %   R_out: 完整的系统残差 (Stiffness + Eddy + Coupling - Source)
            
            A_vec = x_eval(1:ctx.num_A);
            
            % 1. 物理 Stiffness (Curl-Curl)
            [J_nu, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            
            % 确保维度匹配 (Assembler 可能返回压缩的矩阵)
            if size(R_nu, 1) < ctx.numTotalDofs
                R_nu(ctx.numTotalDofs, 1) = 0; 
            end
            if calc_J && size(J_nu, 1) < ctx.numTotalDofs
                [ji, jj, jv] = find(J_nu);
                J_nu = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            end
            if ~calc_J, J_nu = []; end
            
            % 2. 瞬态项 (Eddy Current): M * (x - x_prev) / dt
            dA_dt_term = ctx.M_sigma * (x_eval - ctx.x_prev) / ctx.dt;
            
            % 3. 基础残差
            R_out = R_nu + dA_dt_term - ctx.F_sys;
            
            % 4. 耦合项
            if ctx.enable_V
                % A-Equation: + C_AV * V
                R_out = R_out + ctx.C_AV * x_eval;
                
                % V-Equation: C_AV' * (A - A_prev)/dt + L_VV * V
                R_V_rows = (ctx.C_AV' * (x_eval - ctx.x_prev)) / ctx.dt + ctx.L_VV * x_eval;
                R_out = R_out + R_V_rows;
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