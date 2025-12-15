classdef SDCSolver < handle
    % SDCSOLVER 基于谱延迟校正 (SDC) 的高阶场路耦合瞬态求解器 (v7.16 - MUMPS Fix)
    % 
    % 更新日志:
    %   v7.16: 适配 LinearSolver v3.9，在 factorize/solve 分离模式下
    %          正确存储和传递系统矩阵 A_bc，解决了 dmumps 参数不足的错误。
    
    properties
        Assembler       % 有限元组装器对象
        LinearSolver    % 线性方程组求解器 (需支持复数求解)
        PolyOrder = 3        
        MaxSDCIters = 5      
        SDCTolerance = 1e-6  
        MaxNewtonIters = 15  
        NewtonTolerance = 1e-4 
        UseLineSearch = true
    end
    
    methods
        function obj = SDCSolver(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 300; 
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle, ...
                                                 probePoint)
            % SOLVE 执行 SDC 瞬态求解主循环
            
            fprintf('==================================================\n');
            fprintf('   SDC High-Order Solver (P=%d, SDC_Iters=%d)     \n', obj.PolyOrder, obj.MaxSDCIters);
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
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
            
            postProc = [];
            probeB_History = []; 
            if ~isempty(probePoint)
                fprintf('   [Init] Initializing B-Field Probe at [%.3f, %.3f, %.3f]...\n', ...
                    probePoint(1), probePoint(2), probePoint(3));
                postProc = PostProcessor(obj.Assembler);
                probeB_History = zeros(length(timeSteps), 1);
            end
            
            % 预计算谱时间基础数据
            fprintf('   [Init] Pre-computing Spectral Basis (GLL Order %d)...\n', obj.PolyOrder);
            [gll_nodes, gll_weights] = SpectralTimeUtils.gll(obj.PolyOrder);
            Q_mat = SpectralTimeUtils.integration_matrix(gll_nodes);
            [S, Lambda, S_inv, pair_map] = SpectralTimeUtils.decompose_for_sdc(Q_mat);
            
            ctx.Spectral.Nodes = gll_nodes;     
            ctx.Spectral.Weights = gll_weights;
            ctx.Spectral.Q = Q_mat;
            ctx.Spectral.S = S;
            ctx.Spectral.Lambda = diag(Lambda); 
            ctx.Spectral.S_inv = S_inv;
            ctx.Spectral.PairMap = pair_map;    
            ctx.Spectral.NumNodes = length(gll_nodes); 
            
            % 组装时不变矩阵
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            dt_init = timeSteps(1); 
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs);
            
            x_start = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_start(1:n) = init_State(1:n);
            end
            
            currentHistory = zeros(length(timeSteps), 1);
            globalTime = 0;
            
            for slab_idx = 1:length(timeSteps)
                dt_slab = timeSteps(slab_idx);
                t_start = globalTime;
                t_end = t_start + dt_slab;
                
                fprintf('\n=== Time Slab %d / %d (dt=%.1e, T=%.4f -> %.4f) ===\n', ...
                    slab_idx, length(timeSteps), dt_slab, t_start, t_end);
                
                ctx.dt_slab = dt_slab;
                ctx.t_start = t_start;
                
                fprintf('   [Step 1] Prediction (Backward Euler Sweep)\n');
                X_slab = obj.predictionStep(x_start, ctx);
                
                fprintf('   [Step 2] Linearization & Decoupling\n');
                DecoupledSystems = obj.linearizationStep(X_slab, ctx);
                
                fprintf('   [Step 3] SDC Correction Loop\n');
                for sdc_iter = 1:obj.MaxSDCIters
                    [X_new, error_norm] = obj.sdcCorrectionStep(X_slab, x_start, DecoupledSystems, ctx);
                    
                    fprintf('      SDC Iter %d: Defect Norm = %.4e', sdc_iter, error_norm);
                    
                    X_slab = X_new;
                    if error_norm < obj.SDCTolerance
                        fprintf(' -> Converged.\n');
                        break;
                    else
                        fprintf('\n');
                    end
                end
                
                % [Phase 4] 清理资源
                obj.cleanUpSystems(DecoupledSystems);
                
                x_start = X_slab(:, end); 
                globalTime = t_end;
                
                currentHistory(slab_idx) = x_start(end);
                
                B_mag_current = 0;
                if ~isempty(postProc)
                    A_sol = x_start(1:ctx.num_A);
                    B_val = postProc.probeB(A_sol, probePoint);
                    B_mag_current = norm(B_val);
                    probeB_History(slab_idx) = B_mag_current;
                    fprintf('      [Probe] |B| at point = %.4f T, Current = %.4f A\n', ...
                        B_mag_current, full(x_start(end)));
                end
                
                if ~isempty(monitorHandle)
                    try
                        t_vec = cumsum(timeSteps(1:slab_idx));
                        I_vec = currentHistory(1:slab_idx);
                        
                        if ~isempty(postProc) && nargin(monitorHandle) >= 6
                            B_vec_hist = probeB_History(1:slab_idx);
                            monitorHandle(globalTime, x_start(end), t_vec, I_vec, B_mag_current, B_vec_hist);
                        else
                            monitorHandle(globalTime, x_start(end), t_vec, I_vec);
                        end
                        drawnow limitrate;
                    catch
                    end
                end
            end
            
            solutionResults = x_start;
            info.FinalTime = globalTime;
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc), info.ProbeB_History = probeB_History; end
        end
    end
    
    methods (Access = private)
        function X_pred = predictionStep(obj, x0, ctx)
            % 预测步 (Backward Euler)
            num_nodes = ctx.Spectral.NumNodes; 
            num_dofs = ctx.numTotalDofs;
            
            X_pred = zeros(num_dofs, num_nodes);
            X_pred(:, 1) = x0;
            
            x_prev = x0;
            t_start_slab = ctx.t_start;
            dt_slab = ctx.dt_slab;
            nodes = ctx.Spectral.Nodes;
            
            for m = 2:num_nodes
                tau_curr = nodes(m);
                tau_prev = nodes(m-1);
                
                t_curr = t_start_slab + (tau_curr + 1)/2 * dt_slab;
                dt_sub = (tau_curr - tau_prev)/2 * dt_slab;
                
                fprintf('      - Node %d/%d (dt=%.1e): ', m, num_nodes, dt_sub);
                
                V_source = ctx.circuit.V_source_func(t_curr);
                [x_next, ~] = obj.solveBackwardEulerStep(x_prev, dt_sub, V_source, ctx);
                
                X_pred(:, m) = x_next;
                x_prev = x_next;
            end
        end
        
        function Systems = linearizationStep(obj, X_ref, ctx)
            % 线性化与解耦
            Systems = struct();
            fprintf('      - Computing Frozen Jacobian (from Slab End-Point)...\n');
            
            % 使用终值线性化
            x_ref_state = X_ref(:, end);
            A_ref = x_ref_state(1:ctx.num_A);
            
            [K_fem_small, ~] = obj.Assembler.assembleJacobian(ctx.space_A, A_ref, ctx.matLib, true);
            [ki, kj, kv] = find(K_fem_small);
            K_fem = sparse(ki, kj, kv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            Scale = ctx.CircuitRowScale;
            [w_idx, ~, w_val] = find(ctx.W_vec);
            num_total = ctx.numTotalDofs;
            
            M_cir_row = sparse(repmat(num_total, length(w_idx), 1), w_idx, w_val * Scale, num_total, num_total);
            M_cir_diag = sparse(num_total, num_total, ctx.circuit.L * Scale, num_total, num_total);
            M_sys = ctx.M_sigma + M_cir_row + M_cir_diag;
            
            K_base = K_fem + ctx.C_AP + ctx.C_AP';
            K_col_coup = sparse(w_idx, repmat(num_total, length(w_idx), 1), -w_val, num_total, num_total);
            K_cir_diag = sparse(num_total, num_total, ctx.circuit.R * Scale, num_total, num_total);
            K_sys = K_base + K_col_coup + K_cir_diag;
            
            pair_map = ctx.Spectral.PairMap;
            lambda = ctx.Spectral.Lambda;
            dt_factor = ctx.dt_slab / 2.0; 
            
            Systems.Factors = cell(length(lambda), 1);
            count_factorized = 0;
            
            for m = 1:length(lambda)
                type = pair_map(m).Type;
                lam = lambda(m);
                
                if abs(lam) < 1e-9
                    Systems.Factors{m} = []; 
                    continue;
                end
                
                if strcmp(type, 'Real') || strcmp(type, 'ComplexA')
                    alpha = lam * dt_factor;
                    A_complex = M_sys + alpha * K_sys;
                    
                    [A_bc, ~] = BoundaryCondition.applyDirichlet(A_complex, zeros(num_total,1), ctx.fixedDofs);
                    
                    try
                        % [Fix] 存储分解因子 ID 和 矩阵 A_bc
                        mumpsID = obj.LinearSolver.factorize(A_bc);
                        Systems.Factors{m} = struct('id', mumpsID, 'mat', A_bc);
                        count_factorized = count_factorized + 1;
                    catch ME
                        warning('SDC Subsystem %d factorization failed: %s', m, ME.message);
                        Systems.Factors{m} = [];
                    end
                end
            end
            fprintf('      - Factorized %d unique subsystems (MUMPS Optimized).\n', count_factorized);
            
            Systems.M_sys = M_sys;
            Systems.K_sys = K_sys; 
        end
        
        function [X_new, error_norm] = sdcCorrectionStep(obj, X_curr, x0, Systems, ctx)
            % SDCCORRECTIONSTEP 执行一次完整的 SDC 修正扫描 (带自适应阻尼)
            % 
            % [v7.17 Fix] 引入阻尼因子 (Damping Factor)
            % 原因: 在极强饱和区 (Slab 4)，冻结 Jacobian 与实际刚度偏差过大，
            %       全量更新 (alpha=1.0) 会导致震荡发散。
            
            num_nodes = ctx.Spectral.NumNodes;
            num_dofs = ctx.numTotalDofs;
            dt_half = ctx.dt_slab / 2.0;
            
            % 1. 计算物理残差 F(t, x)
            F_all = zeros(num_dofs, num_nodes);
            Scale = ctx.CircuitRowScale;
            C_AP_tot = ctx.C_AP + ctx.C_AP';
            W_vec = ctx.W_vec;
            R_cir = ctx.circuit.R;
            
            nodes = ctx.Spectral.Nodes;
            t_physical = ctx.t_start + (nodes + 1) * dt_half;
            
            for m = 1:num_nodes
                x_m = X_curr(:, m);
                t_m = t_physical(m);
                
                A_vec = x_m(1:ctx.num_A);
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                
                R_fem_ext = sparse(num_dofs, 1);
                R_fem_ext(1:length(R_fem)) = R_fem;
                F_val = -R_fem_ext - C_AP_tot * x_m + W_vec * x_m(end);
                
                V_source = ctx.circuit.V_source_func(t_m);
                I_val = x_m(end);
                F_cir = Scale * (V_source - R_cir * I_val);
                F_val(end) = F_cir;
                
                F_val(ctx.fixedDofs) = 0;
                F_all(:, m) = F_val;
            end
            
            % 2. 谱积分
            Q = ctx.Spectral.Q;
            Integral_F = (F_all * Q.') * dt_half;
            
            % 3. Defect (Picard 误差)
            M_sys = Systems.M_sys;
            M_x0 = M_sys * x0; 
            M_X = M_sys * X_curr;
            Defect = M_x0 + Integral_F - M_X;
            
            % 4. 变换到特征空间
            S_inv = ctx.Spectral.S_inv;
            G = Defect * S_inv.';
            
            % 5. 并行求解 (修正量 Delta_Y)
            Delta_Y = zeros(size(G));
            pair_map = ctx.Spectral.PairMap;
            Factors = Systems.Factors;
            
            for m = 1:num_nodes
                ptype = pair_map(m).Type;
                
                if isempty(Factors{m})
                    Delta_Y(:, m) = 0;
                    continue; 
                end
                
                if strcmp(ptype, 'Real') || strcmp(ptype, 'ComplexA')
                    rhs = G(:, m);
                    
                    sys_data = Factors{m};
                    dy = obj.LinearSolver.solveFromFactors(sys_data.id, rhs, sys_data.mat);
                    
                    Delta_Y(:, m) = dy;
                    
                    if strcmp(ptype, 'ComplexA')
                        idx_B = pair_map(m).ConjIdx;
                        Delta_Y(:, idx_B) = conj(dy);
                    end
                end
            end
            
            % 6. 反变换回物理空间
            S = ctx.Spectral.S;
            Delta_X = Delta_Y * S.';
            Delta_X = real(Delta_X);
            
            % 7. 更新 (带阻尼)
            % [STRATEGY] 使用简单的阻尼策略。
            % 初始几次迭代使用较小阻尼以稳定，后期全量更新。
            % 或者监测 Defect Norm 是否上升（这里简化为固定阻尼）
            
            % 对于极强非线性问题，0.5 ~ 0.8 是比较安全的范围
            damping = 0.8; 
            
            X_new = X_curr + damping * Delta_X;
            
            % 强制边界和初始条件
            X_new(ctx.fixedDofs, :) = repmat(x0(ctx.fixedDofs), 1, num_nodes);
            X_new(:, 1) = x0; 
            
            % 误差范数基于全量 Delta_X 计算
            error_norm = norm(Delta_X, 'fro') / sqrt(num_nodes * num_dofs);
        end
        
        function cleanUpSystems(obj, Systems)
            % [Fix] 释放 MUMPS 内存 (适配 struct 结构)
            if ~isfield(Systems, 'Factors'), return; end
            Factors = Systems.Factors;
            for k = 1:length(Factors)
                if ~isempty(Factors{k})
                    % struct 中包含 id 字段
                    obj.LinearSolver.clearFactors(Factors{k}.id);
                end
            end
        end
        
        function [x_new, converged] = solveBackwardEulerStep(obj, x_prev, dt, V_source, ctx)
            % 辅助函数: 求解单步 Backward Euler (带回溯线搜索)
            x_curr = x_prev;
            converged = false;
            tol = obj.NewtonTolerance; 
            
            % 初始残差计算
            [Res, J_sys] = obj.computeResidual(x_curr, x_prev, dt, V_source, ctx, true);
            res_norm = norm(Res);
            
            for iter = 1:obj.MaxNewtonIters
                % 1. 求解牛顿步
                [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res, ctx.fixedDofs);
                dx = obj.LinearSolver.solve(J_bc, -Res_bc);
                
                % 2. 回溯线搜索 (Backtracking Line Search)
                alpha = 1.0;
                dec_factor = 0.5;
                min_alpha = 1e-2; % 防止死循环
                accepted = false;
                
                res_norm_prev = res_norm;
                
                while alpha > min_alpha
                    x_try = x_curr + alpha * dx;
                    
                    % 仅计算残差 (calcJ=false) 以节省时间
                    [Res_try, ~] = obj.computeResidual(x_try, x_prev, dt, V_source, ctx, false);
                    
                    % 强制边界残差为0 (不计入范数)
                    Res_try(ctx.fixedDofs) = 0; 
                    res_norm_new = norm(Res_try);
                    
                    % Armijo 准则简化版: 只要残差下降即可接受
                    if res_norm_new < res_norm_prev
                        x_curr = x_try;
                        Res = Res_try;
                        res_norm = res_norm_new;
                        accepted = true;
                        
                        % 如果是第一步且 alpha < 1，说明非线性极强
                        % if alpha < 1.0, fprintf('      [LS] Backtrack alpha=%.2f\n', alpha); end
                        break;
                    end
                    alpha = alpha * dec_factor;
                end
                
                % 如果线搜索失败，强制接受最小步长（避免死锁），但很可能不收敛
                if ~accepted
                    x_curr = x_curr + min_alpha * dx;
                    % 重新计算完整 Jacobian 供下一步使用
                    [Res, J_sys] = obj.computeResidual(x_curr, x_prev, dt, V_source, ctx, true);
                    res_norm = norm(Res);
                else
                    % 如果接受了，且需要进行下一次迭代，需更新 Jacobian
                    % 注意: 如果已经收敛，则不需要这一步，但在循环头判断更清晰
                    if res_norm > tol
                         [Res, J_sys] = obj.computeResidual(x_curr, x_prev, dt, V_source, ctx, true);
                    end
                end
                
                % 3. 收敛判断 (使用更新量 dx 和 残差 res_norm 双重判断)
                % 注意: dx 是全步长，如果 alpha 小，实际移动量是 alpha*dx
                if res_norm < tol || (norm(dx)*alpha < tol && iter > 1)
                    converged = true; 
                    fprintf('- Newton Converged (Iter %d, Res=%.2e)\n', iter, res_norm);
                    break; 
                end
            end
            
            if ~converged
                fprintf('      [Warn] Prediction Newton Max Iter Reached (Res=%.2e)\n', res_norm);
            end
            x_new = x_curr;
        end
        
        function [Res, J_sys] = computeResidual(obj, x_curr, x_prev, dt, V_source, ctx, calcJ)
            % 统一计算残差和雅可比矩阵
            % calcJ = false 时，J_sys 返回空，加速线搜索
            
            A_vec = x_curr(1:ctx.num_A);
            
            if calcJ
                [J_fem, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
                [ji, jj, jv] = find(J_fem);
                J_sys = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            else
                [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
                J_sys = [];
            end
            
            Res = sparse(ctx.numTotalDofs, 1);
            Res(1:length(R_fem)) = R_fem;
            
            % 瞬态项 (M/dt)
            dx_dt = (x_curr - x_prev) / dt;
            Res = Res + ctx.M_sigma * dx_dt;
            if calcJ
                J_sys = J_sys + ctx.M_sigma / dt; 
            end
            
            % 规范项 (C_AP)
            % 注意: C_AP 是对称的，K_sys += C + C'
            C_tot = ctx.C_AP + ctx.C_AP';
            Res = Res + C_tot * x_curr;
            if calcJ
                J_sys = J_sys + C_tot;
            end
            
            % 电路耦合
            I_val = x_curr(end);
            [w_idx, ~, w_val] = find(ctx.W_vec);
            Scale = ctx.CircuitRowScale;
            
            % FEM 方程中的电路项 (-W * I)
            Res = Res - ctx.W_vec * I_val;
            if calcJ
                % 列耦合
                J_col = sparse(w_idx, repmat(ctx.numTotalDofs, length(w_idx), 1), -w_val, ctx.numTotalDofs, ctx.numTotalDofs);
                % 行耦合 (W'/dt * A) * Scale
                J_row = sparse(repmat(ctx.numTotalDofs, length(w_idx), 1), w_idx, (w_val/dt)*Scale, ctx.numTotalDofs, ctx.numTotalDofs);
                J_sys = J_sys + J_col + J_row;
            end
            
            % 电路方程 (V - RI - L dI/dt - dPsi/dt = 0)
            % 乘以 Scale
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt; % dPsi/dt approx
            
            Res(end) = (term_R + term_L + term_EMF - V_source) * Scale;
            
            if calcJ
                 J_sys(end, end) = J_sys(end, end) + (ctx.circuit.R + ctx.circuit.L/dt) * Scale;
            end
            
            % 边界条件残差置零 (不影响 LinearSolver 求解，但影响 Norm 判断)
            Res(ctx.fixedDofs) = 0;
        end
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            % 复用 TransientCoupledSolver 的 GeoBalance 逻辑
            fprintf('   [Init] Calculating smart circuit scaling (GeoBalance)...\n');
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
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
            Z_cir = abs(ctx.circuit.R + ctx.circuit.L/dt); 
            max_W = max(abs(nonzeros(ctx.W_vec)));         
            max_coupling_val = max_W / dt;               
            max_circuit_val = max(Z_cir, max_coupling_val);
            if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
            scaleFactor = k_target / max_circuit_val;
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
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
            % 辅助函数: 采样非线性材料的 B-H 曲线，查找最大相对导磁率
            max_mu_r = 1000.0; % 默认值
            try
                b_samples = linspace(0.01, 2.5, 50); % 采样范围 B = 0.01 ~ 2.5 T
                b_sq_samples = b_samples.^2;
                nu_vals = zeros(size(b_samples));
                
                % 计算每个点的 nu(B^2)
                for i = 1:length(b_samples)
                    [nu, ~] = MaterialLib.evaluate(b_sq_samples(i), mat);
                    nu_vals(i) = nu;
                end
                
                % 转换回相对导磁率 mu_r = 1 / (nu * mu0)
                mu0 = 4*pi*1e-7;
                mu_r_vals = 1 ./ (nu_vals * mu0);
                max_mu_r = max(mu_r_vals);
                
                if max_mu_r < 1, max_mu_r = 1; end
            catch
                warning('Failed to find Max Mu_r. Using default 1000.');
            end
        end
    end
end