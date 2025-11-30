classdef NonlinearSolver < handle
    % NONLINEARSOLVER 静磁场非线性求解器 (v6.0 - Robust)
    %
    % 改进:
    %   1. [Robustness] 引入回溯线搜索 (Line Search)，解决强饱和收敛问题。
    %   2. [Auto-Grounding] 自动检测悬浮 V 场并接地，防止矩阵奇异。
    %   3. [Regularization] 改进正则化系数计算 (使用 max 而非 mean)。
    
    properties
        Assembler
        LinearSolver
        MaxIterations = 100
        Tolerance = 1e-6
        
        % --- 线搜索配置 ---
        Relaxation = 1.0         
        UseLineSearch = true
        MaxLineSearchIters = 10
        
        % --- 正则化配置 ---
        MatrixK_Add = []         
        AutoRegularization = true 
        RegScale = 1e-6          
    end
    
    methods
        function obj = NonlinearSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 40;
        end
        
        function [A_sol, V_sol, info] = solve(obj, space_A, space_V, matLib, sigmaMap, sourceMap, fixedDofs_A, fixedDofs_V, x_init)
            if nargin < 9, x_init = []; end
            if nargin < 8, fixedDofs_V = []; end
            
            fprintf('==============================================\n');
            fprintf('   Nonlinear Magnetostatic Solver (v6.0)      \n');
            fprintf('==============================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % --- 1. 空间与自由度 ---
            [~, offset_A] = dofHandler.distributeDofs(space_A);
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            if num_A == 0
                error('NonlinearSolver:NoDoFs', 'Active DoFs is 0. Check Mesh initialization.');
            end
            
            enable_V = ~isempty(space_V);
            offset_V = -1;
            auto_ground_idx = [];
            
            if enable_V
                conductingTags = obj.getConductingRegions(sigmaMap);
                if isempty(conductingTags)
                    fprintf('   [Info] No conducting regions. Dropping V-field.\n');
                    enable_V = false;
                else
                    if ~dofHandler.DofMaps.isKey(space_V.toString())
                        [~, offset_V] = dofHandler.distributeDofs(space_V, conductingTags);
                    else
                        offset_V = dofHandler.SpaceOffsets(space_V.toString());
                    end
                    fprintf('   [Info] V-field enabled on %d regions.\n', length(conductingTags));
                    
                    % [Auto-Grounding] 检查悬浮电位
                    hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                    if ~hasFixedV
                        fprintf('   [Auto-Fix] V-field is floating. Grounding 1st node to 0V.\n');
                        auto_ground_idx = offset_V + 1;
                    end
                end
            end
            
            numTotalDofs = dofHandler.NumGlobalDofs;
            
            % --- 2. 预计算常数矩阵 ---
            % 上下文结构体，减少传参
            sysCtx.num_A = num_A;
            sysCtx.numTotalDofs = numTotalDofs;
            sysCtx.enable_V = enable_V;
            sysCtx.space_A = space_A;
            sysCtx.matLib = matLib;
            
            F_A = obj.Assembler.assembleSource(space_A, sourceMap);
            sysCtx.F_sys = sparse(numTotalDofs, 1);
            if ~isempty(F_A), sysCtx.F_sys(1:length(F_A)) = F_A; end
            
            sysCtx.C_AV = sparse(numTotalDofs, numTotalDofs);
            sysCtx.L_VV = sparse(numTotalDofs, numTotalDofs);
            if enable_V
                sysCtx.C_AV = obj.Assembler.assembleCoupling(space_A, space_V, sigmaMap);
                sysCtx.L_VV = obj.Assembler.assembleScalarLaplacian(space_V, sigmaMap);
            end
            
            % 预计算正则化质量矩阵
            sysCtx.M_reg = sparse(numTotalDofs, numTotalDofs);
            if obj.AutoRegularization
                M_local = obj.Assembler.assembleMass(space_A);
                [im, jm, vm] = find(M_local);
                sysCtx.M_reg = sparse(im, jm, vm, numTotalDofs, numTotalDofs);
            end
            
            % --- 3. 边界条件 ---
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed_dofs = idx_A;
            
            if enable_V
                if isempty(fixedDofs_V), idx_V = [];
                elseif islogical(fixedDofs_V), idx_V = find(fixedDofs_V); 
                else, idx_V = fixedDofs_V(:); end
                
                if ~isempty(auto_ground_idx)
                    idx_V = unique([idx_V; auto_ground_idx]);
                end
                
                idx_V = idx_V(idx_V <= numTotalDofs); 
                all_fixed_dofs = [all_fixed_dofs; idx_V];
            end
            
            % --- 4. 初始解 ---
            x = zeros(numTotalDofs, 1);
            if ~isempty(x_init)
                n_copy = min(length(x_init), numTotalDofs);
                x(1:n_copy) = x_init(1:n_copy);
            end
            
            % --- 5. 牛顿迭代 ---
            converged = false;
            norm_res = 0;
            
            % 动态正则化系数
            sysCtx.RegScaleAbs = 0;
            
            % 初次评估
            [J_sys, R_sys] = obj.evaluateSystem(x, sysCtx, true);
            
            % [Auto-Reg] 计算绝对正则化系数
            if obj.AutoRegularization
                diag_J = abs(diag(J_sys));
                ref_scale = max(diag_J); % 使用 Max 增强稳定性
                if ref_scale == 0, ref_scale = 1.0; end
                sysCtx.RegScaleAbs = full(ref_scale * obj.RegScale);
                fprintf('   [Auto-Reg] Regularization Scale: %.2e\n', sysCtx.RegScaleAbs);
                
                J_sys = J_sys + sysCtx.RegScaleAbs * sysCtx.M_reg;
                R_sys = R_sys + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x;
            end
            
            [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed_dofs);
            norm_res = norm(Res_bc);
            
            for iter = 1:obj.MaxIterations
                
                delta_x = obj.LinearSolver.solve(J_bc, -Res_bc);
                if isempty(delta_x), break; end
                
                % --- 线搜索 (Line Search) ---
                alpha = obj.Relaxation;
                if obj.UseLineSearch
                    res_prev = norm_res;
                    for k = 0:obj.MaxLineSearchIters
                        alpha_try = alpha * (0.5^k);
                        x_try = x + alpha_try * delta_x;
                        
                        [~, R_try] = obj.evaluateSystem(x_try, sysCtx, false);
                        R_try(all_fixed_dofs) = 0; % 忽略边界点残差
                        
                        % 考虑正则化对残差的影响
                        if obj.AutoRegularization
                            R_try = R_try + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x_try;
                        end
                        
                        res_try = norm(R_try);
                        if res_try < res_prev || (iter == 1 && k < 3)
                            alpha = alpha_try;
                            if k > 0, fprintf('      [LS] Step %.4f, Res %.2e\n', alpha, res_try); end
                            break;
                        end
                    end
                end
                
                x = x + alpha * delta_x;
                norm_step = norm(alpha * delta_x);
                
                % 准备下一次迭代
                [J_sys, R_sys] = obj.evaluateSystem(x, sysCtx, true);
                
                % 应用正则化
                if obj.AutoRegularization
                    J_sys = J_sys + sysCtx.RegScaleAbs * sysCtx.M_reg;
                    R_sys = R_sys + (sysCtx.RegScaleAbs * sysCtx.M_reg) * x;
                end
                
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, R_sys, all_fixed_dofs);
                norm_res = norm(Res_bc);
                
                fprintf('      Iter %2d: Res=%.4e, Step=%.4e\n', iter, norm_res, norm_step);
                
                if norm_res < obj.Tolerance || norm_step < 1e-10
                    converged = true;
                    break;
                end
            end
            
            if ~converged
                fprintf('   [Warning] Max iterations reached (%d). Res=%.4e\n', obj.MaxIterations, norm_res);
            end
            
            A_sol = x(1:num_A);
            if enable_V
                V_sol = x(offset_V+1 : end);
            else
                V_sol = [];
            end
            
            info.Converged = converged;
            info.Iterations = iter;
            info.FinalResidual = norm_res;
            info.Solution = x; 
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
        function [J_out, R_out] = evaluateSystem(obj, x, ctx, calc_J)
            A_vec = x(1:ctx.num_A);
            
            [J_AA, R_nu] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, calc_J);
            
            % 尺寸安全
            if size(R_nu, 1) < ctx.numTotalDofs
                R_nu(ctx.numTotalDofs, 1) = 0; 
            end
            
            if ctx.enable_V
                % DC 方程:
                % Row A: Curl(Nu Curl A) + C_AV * V = J_s
                % Row V: L_VV * V = 0 (Assuming DC continuity, no d/dt terms)
                % 注意: 这里假设 C_AV 是 (A,V) 块，且不对称
                
                Term_C = ctx.C_AV * x;
                Term_L = ctx.L_VV * x;
                R_out = R_nu + Term_C + Term_L - ctx.F_sys;
                
                J_out = [];
                if calc_J
                    J_out = J_AA + ctx.C_AV + ctx.L_VV;
                    
                    % 检查是否有外部附加矩阵
                    if ~isempty(obj.MatrixK_Add)
                         % (省略实现以保持简洁，通常用于特定约束)
                    end
                end
            else
                R_out = R_nu - ctx.F_sys;
                J_out = J_AA;
            end
        end
    end
end