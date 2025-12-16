classdef TransientBDF2Solver < handle
    % TRANSIENTBDF2SOLVER 场路耦合非线性瞬态求解器 (BDF2 二阶精度版)
    % 
    % 描述:
    %   本类实现了二维/三维瞬态磁场与外电路的强耦合求解。
    %   与基类不同，本类采用二阶向后差分公式 (BDF2) 进行时间步进，
    %   相比 BDF1 (Backward Euler) 具有更高的代数精度和更小的数值耗散。
    %
    % 核心算法特性:
    %   1. 时间积分: 
    %       - 第1步: Backward Euler (自启动)
    %       - 后续步: BDF2 (二阶精度: 1.5*x_n - 2*x_{n-1} + 0.5*x_{n-2}) / dt
    %   2. 其他特性(Newton-Raphson, GeoBalance, LineSearch) 与原版保持一致。
    
    properties
        Assembler       % 有限元组装器对象
        LinearSolver    % 线性方程组求解器对象
        
        % --- 迭代控制参数 ---
        MaxIterations = 30        % 最大非线性迭代次数
        Tolerance = 1e-3          % 绝对残差容差
        RelTolerance = 1e-3       % 相对残差容差
        
        % --- 线搜索控制参数 ---
        UseLineSearch = true      % 开关：是否启用自适应线搜索
        MaxLineSearchIters = 10   % 线搜索最大回溯次数
        MinRelaxation = 1e-4      % 最小允许松弛因子
    end
    
    methods
        function obj = TransientBDF2Solver(assembler)
            % 构造函数: 初始化求解器实例
            obj.Assembler = assembler;
            
            % 初始化线性求解器 (MUMPS)
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 200; 
            obj.LinearSolver.MumpsSymmetry = 0;    
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            % SOLVE 执行瞬态求解主循环 (BDF2 版本)
            
            fprintf('==================================================\n');
            fprintf('   Transient Solver (BDF2 High-Order Stepping)    \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % 1. 自由度检查与分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            
            ctx.space_A = space_A;
            ctx.space_P = space_P;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            if nargin < 12, monitorHandle = []; end
            if nargin < 13, probePoint = []; end
            
            % 初始化磁密探测器
            postProc = [];
            probeB_History = [];
            if ~isempty(probePoint)
                fprintf('   [Init] Initializing B-Field Probe at [%.3f, %.3f, %.3f]...\n', ...
                    probePoint(1), probePoint(2), probePoint(3));
                postProc = PostProcessor(obj.Assembler);
                probeB_History = zeros(length(timeSteps), 1);
            end
            
            % 2. 组装时不变矩阵
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            fprintf('   [Init] Assembling Coulomb gauge constraint...\n');
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            % 3. 计算自适应电路缩放因子
            dt_init = timeSteps(1);
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            
            % 4. 合并边界条件索引
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs); 
            
            % 5. 初始化状态向量
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_curr(1:n) = init_State(1:n);
            end
            
            % --- BDF2 状态初始化 ---
            % x_prev : n-1 时刻解
            % x_prev2: n-2 时刻解
            ctx.x_prev = x_curr; 
            ctx.x_prev2 = x_curr; 
            
            numSteps = length(timeSteps);
            solutionResults = []; 
            currentHistory = zeros(numSteps, 1); 
            currentTime = 0;
            
            % ==================================================
            % 6. 时间步循环
            % ==================================================
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                % --- 计算 BDF 系数 ---
                % dx/dt = alpha * x_curr + beta
                % 其中 beta 包含了历史项 (x_prev, x_prev2)
                if t_idx == 1
                    % 第一步：使用 BDF1 (Backward Euler) 启动
                    % dx/dt = (x_n - x_{n-1}) / dt
                    %       = (1/dt)*x_n + (-1/dt)*x_{n-1}
                    alpha_t = 1.0 / dt;
                    beta_t_vec = -ctx.x_prev / dt;
                    fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f, BDF1 Start) ---\n', t_idx, numSteps, dt, currentTime);
                else
                    % 后续步：使用 BDF2
                    % dx/dt = (1.5*x_n - 2*x_{n-1} + 0.5*x_{n-2}) / dt
                    %       = (1.5/dt)*x_n + (-2*x_{n-1} + 0.5*x_{n-2})/dt
                    alpha_t = 1.5 / dt;
                    beta_t_vec = (-2.0 * ctx.x_prev + 0.5 * ctx.x_prev2) / dt;
                    fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f, BDF2) ---\n', t_idx, numSteps, dt, currentTime);
                end
                
                % 将 BDF 系数存入上下文
                ctx.bdf_alpha = alpha_t;
                ctx.bdf_beta_vec = beta_t_vec;
                ctx.V_source_val = circuitProps.V_source_func(currentTime);
                
                res_norm_init = 0;
                
                % ==============================================
                % 7. 牛顿迭代循环
                % ==============================================
                for iter = 1:obj.MaxIterations
                    
                    % [Step A] 组装 Newton 系统
                    [J_sys, Res_vec] = obj.assembleNewtonSystem(x_curr, ctx);
                    
                    Res_check = Res_vec;
                    Res_check(ctx.fixedDofs) = 0;
                    current_res_norm = norm(Res_check);
                    
                    if iter == 1
                        res_norm_init = current_res_norm;
                        if res_norm_init < 1e-12, res_norm_init = 1.0; end
                    end
                    rel_res = current_res_norm / res_norm_init;
                    
                    I_val = x_curr(end);
                    fprintf('    Iter %d: Res = %.4e, RelRes = %.4e, I = %.6e', ...
                            iter, current_res_norm, rel_res, full(I_val));
                    
                    if current_res_norm < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        fprintf(' -> Converged.\n');
                        break;
                    end
                    fprintf('\n');
                    
                    % [Step B] 求解线性方程组
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res_vec, ctx.fixedDofs);
                    
                    lastwarn(''); 
                    dx = obj.LinearSolver.solve(J_bc, -Res_bc); 
                    
                    [warnMsg, ~] = lastwarn;
                    if contains(warnMsg, 'singular') || contains(warnMsg, 'badly scaled') || any(isnan(dx))
                        fprintf('      [Warn] Matrix ill-conditioned. Adding regularization shift.\n');
                        J_fix = J_bc + 1e-8 * speye(size(J_bc));
                        dx = obj.LinearSolver.solve(J_fix, -Res_bc);
                    end
                    
                    % [Step C] 自适应线搜索
                    alpha = 1.0; 
                    if obj.UseLineSearch
                        accepted = false;
                        for k = 0:obj.MaxLineSearchIters
                            x_try = x_curr + alpha * dx;
                            try_res_norm = obj.computeResidualNorm(x_try, ctx);
                            
                            if try_res_norm < current_res_norm * 1.0001
                                if k > 0
                                    fprintf('      [LS] Backtrack: Alpha=%.4f (Res: %.2e -> %.2e)\n', ...
                                        alpha, current_res_norm, try_res_norm);
                                end
                                accepted = true; 
                                break;
                            else
                                alpha = alpha * 0.5;
                            end
                            
                            if alpha < obj.MinRelaxation
                                accepted = true; 
                                break; 
                            end
                        end
                        if ~accepted; alpha = obj.MinRelaxation; end
                    end
                    
                    % [Step D] 更新解向量
                    x_curr = x_curr + alpha * dx;
                end
                
                % --- [BDF2 更新逻辑] ---
                % 推进时间步历史：x_{n-2} = x_{n-1}, x_{n-1} = x_n
                ctx.x_prev2 = ctx.x_prev; 
                ctx.x_prev = x_curr; 
                
                % 记录结果
                solutionResults = x_curr;
                currentHistory(t_idx) = x_curr(end);
                
                % 磁密探针计算
                B_mag_current = 0;
                if ~isempty(postProc)
                    A_sol_step = x_curr(1:ctx.num_A);
                    B_vec = postProc.probeB(A_sol_step, probePoint);
                    B_mag_current = norm(B_vec);
                    probeB_History(t_idx) = B_mag_current;
                    fprintf('    [Probe] |B| at [%.2f,%.2f,%.2f] = %.4f T\n', ...
                        probePoint(1), probePoint(2), probePoint(3), B_mag_current);
                end
                
                % 实时监控与绘图
                if ~isempty(monitorHandle)
                    try
                        t_vec = cumsum(timeSteps(1:t_idx));
                        I_vec = currentHistory(1:t_idx);
                        
                        if ~isempty(postProc) && nargin(monitorHandle) >= 6
                            B_vec_hist = probeB_History(1:t_idx);
                            monitorHandle(currentTime, x_curr(end), t_vec, I_vec, B_mag_current, B_vec_hist);
                        else
                            monitorHandle(currentTime, x_curr(end), t_vec, I_vec);
                        end
                        grid on; drawnow; 
                    catch ME
                        fprintf('Monitor handle execution failed: %s', ME.message);
                    end
                end
            end
            
            info.FinalTime = currentTime; 
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc)
                info.ProbeB_History = probeB_History;
            end
        end
    end
    
    % ==================================================
    % 私有辅助方法 (Private Methods)
    % ==================================================
    methods (Access = private)
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            % CALCULATECIRCUITSCALE 计算自适应电路方程缩放因子
            fprintf('   [Init] Calculating smart circuit scaling (GeoBalance)...\n');
            
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            
            % 注意: 在 BDF2 中，有效刚度项是 alpha * M。
            % 粗略估计 scaling 时，仍使用 1/dt 量级即可，或者使用 1.5/dt
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            
            nu_vec_soft = obj.getBoundNuVec(ctx, 'soft');
            K_soft = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_soft);
            soft_vals = full(abs(diag(K_soft)) + abs(diag_M(1:num_fem)) / dt);
            
            meshTags = obj.Assembler.Mesh.RegionTags;
            matLib = ctx.matLib;
            dofMap = obj.Assembler.DofHandler.DofMaps(ctx.space_A.toString());
            
            target_tags = [];
            k = matLib.keys;
            for i = 1:length(k)
                tag = k{i};
                if ischar(tag), tag_val = str2double(tag); else, tag_val = tag; end
                mat = matLib(tag);
                if strcmpi(mat.Type, 'Nonlinear'), target_tags = [target_tags, tag_val]; end
            end
            
            if isempty(target_tags)
                k_soft_rep = k_hard_median; 
            else
                mask_elem = ismember(meshTags, target_tags);
                target_dofs = unique(dofMap(:, mask_elem));
                target_dofs = target_dofs(target_dofs > 0);
                k_soft_rep = median(soft_vals(target_dofs));
            end
            
            k_target = sqrt(k_hard_median * k_soft_rep);
            
            % 电路估算
            Z_cir = abs(ctx.circuit.R + ctx.circuit.L/dt); 
            max_W = max(abs(nonzeros(ctx.W_vec)));         
            max_coupling_val = max_W / dt;                 
            max_circuit_val = max(Z_cir, max_coupling_val);
            if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
            
            scaleFactor = k_target / max_circuit_val;
            
            fprintf('     -> Smart Scale Factor: %.2e\n', scaleFactor);
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            % 辅助函数: 获取极限状态下的磁阻率向量 nu
            meshTags = obj.Assembler.Mesh.RegionTags;
            numElems = length(meshTags);
            nu_vec = zeros(numElems, 1);
            matLib = ctx.matLib;
            nu0 = 1 / (4 * pi * 1e-7);
            
            uniqueTags = unique(meshTags);
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k);
                if ~matLib.isKey(tag), continue; end
                mat = matLib(tag);
                val = 0;
                if strcmp(mat.Type, 'Linear')
                    val = mat.Nu_Linear;
                else
                    if strcmp(mode, 'hard')
                        val = nu0;
                    else
                        mu_r_max = obj.findMaxMuR(mat);
                        val = nu0 / mu_r_max;
                    end
                end
                nu_vec(meshTags == tag) = val;
            end
        end
        
        function max_mu_r = findMaxMuR(~, mat)
            max_mu_r = 1000.0;
            try
                b_samples = linspace(0.01, 2.5, 50);
                b_sq_samples = b_samples.^2;
                nu_vals = zeros(size(b_samples));
                for i = 1:length(b_samples)
                    [nu, ~] = MaterialLib.evaluate(b_sq_samples(i), mat);
                    nu_vals(i) = nu;
                end
                mu0 = 4*pi*1e-7;
                mu_r_vals = 1 ./ (nu_vals * mu0);
                max_mu_r = max(mu_r_vals);
                if max_mu_r < 1, max_mu_r = 1; end
            catch
                warning('Failed to find Max Mu_r. Using default 1000.');
            end
        end
        
        function [J_sys, Res_vec] = assembleNewtonSystem(obj, x, ctx)
            % 组装 Newton 系统 (BDF2 版本)
            % dx/dt = alpha * x + beta_vec
            
            % 1. FEM 部分 (A场)
            A_vec = x(1:ctx.num_A);
            [J_fem, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            
            [ji, jj, jv] = find(J_fem);
            J_sys = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            % 2. 瞬态项 (M * dx/dt)
            % dx/dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec
            dx_dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec;
            
            % 涡流残差: + M * (alpha*x + beta)
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            % 涡流雅可比: + alpha * M
            J_sys = J_sys + ctx.bdf_alpha * ctx.M_sigma;
            
            % 规范残差: + (C_AP + C_AP') * x
            J_sys = J_sys + ctx.C_AP + ctx.C_AP';
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            
            % 3. 电路耦合
            I_val = x(end);
            Scale = ctx.CircuitRowScale; 
            [w_idx, ~, w_val] = find(ctx.W_vec);
            
            % 列耦合: -W * I
            J_col = sparse(w_idx, repmat(ctx.numTotalDofs, length(w_idx), 1), -w_val, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = Res_vec - ctx.W_vec * I_val;
            
            % 行耦合: (W' * dx/dt) -> (W' * (alpha*A + beta_A)) * Scale
            % 对应 Jacobian 部分: alpha * W' * Scale
            J_row_entries = w_val * ctx.bdf_alpha * Scale;
            J_row = sparse(repmat(ctx.numTotalDofs, length(w_idx), 1), w_idx, J_row_entries, ctx.numTotalDofs, ctx.numTotalDofs);
            J_sys = J_sys + J_col + J_row;
            
            % 4. 电路对角项与残差
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt; % d/dt (W'*A)
            
            % 电路残差
            res_I_val = (term_R + term_L + term_EMF - ctx.V_source_val) * Scale;
            Res_vec(end) = res_I_val;
            
            % 电路雅可比对角: (R + alpha*L) * Scale
            Z_circuit = (ctx.circuit.R + ctx.circuit.L * ctx.bdf_alpha) * Scale;
            J_sys(end, end) = J_sys(end, end) + Z_circuit;
        end
        
        function res_norm = computeResidualNorm(obj, x, ctx)
            % 快速计算残差范数 (BDF2 版本)
            
            A_vec = x(1:ctx.num_A);
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            % dx/dt = alpha * x + beta_vec
            dx_dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec;
            
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            Res_vec = Res_vec - ctx.W_vec * x(end);
            
            Scale = ctx.CircuitRowScale;
            term_R = ctx.circuit.R * x(end);
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt;
            res_I_val = (term_R + term_L + term_EMF - ctx.V_source_val) * Scale;
            Res_vec(end) = res_I_val;
            
            Res_vec(ctx.fixedDofs) = 0;
            res_norm = norm(Res_vec);
        end
    end
end