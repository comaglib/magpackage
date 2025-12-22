classdef TransientBDF2Solver < handle
    % TRANSIENTBDF2SOLVER 场路耦合非线性瞬态求解器 (二阶 BDF2 多绕组增强版)
    % 
    % 算法特性:
    %   - 时间离散: 二阶向后差分公式 (BDF2)。第一步使用 BDF1 启动，之后切换至 BDF2。
    %   - 特点: 相比 BDF1，具有二阶时间精度，数值阻尼更小，能更精确地捕捉瞬态波形转折。
    %   - 非线性迭代: 带有自适应线搜索的 Newton-Raphson 法，处理变压器铁芯深度饱和。
    %   - 场路耦合: 全耦合 (Monolithic) 策略，保证磁场能与电路能的严格守恒。
    %   - 数值调理: 严格 GeoBalance 逻辑，解决 FEM 方程与电路方程数量级相差悬殊导致的病态问题。
    
    properties
        Assembler       % 有限元组装器对象
        LinearSolver    % 线性直接求解器 (MUMPS)
        
        % --- 迭代控制参数 ---
        MaxIterations = 30        % 最大 Newton 迭代次数
        Tolerance = 1e-3          % 绝对残差容差 (||R||)
        RelTolerance = 1e-3       % 相对残差容差 (||R||/||R0||)
        
        % --- 线搜索控制参数 ---
        UseLineSearch = true      % 是否启用线搜索 (推荐开启)
        MaxLineSearchIters = 10   % 最大回溯次数
        MinRelaxation = 1e-4      % 最小松弛因子 (步长下限)
    end
    
    methods
        function obj = TransientBDF2Solver(assembler)
            % 构造函数: 初始化求解器
            obj.Assembler = assembler;
            
            % 初始化 MUMPS 求解器: 场路耦合是非对称系统，必须 Symmetry=0
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
            
            fprintf('==================================================\n');
            fprintf('   Transient Solver (BDF2) - Multi-Winding Strict \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % 1. 确保自由度已分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            % 2. 系统维度初始化
            numWindings = length(windingObj);
            numFemDofs = dofHandler.NumGlobalDofs; 
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numFemDofs = numFemDofs;
            ctx.numWindings = numWindings;
            ctx.numTotalDofs = numFemDofs + numWindings; 
            
            % 封装计算上下文
            ctx.space_A = space_A;
            ctx.space_P = space_P;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            % 3. 组装时不变物理矩阵 (Invariant Matrices)
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % [M_sigma] 涡流项质量矩阵: 描述 sigma * dA/dt
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % [C_AP] 库伦规范矩阵: 施加 div A = 0 约束
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % [W_mat] 场路耦合系数矩阵: 描述电流密度分布 J
            W_fem_mat = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_mat = sparse(ctx.numTotalDofs, numWindings);
            ctx.W_mat(1:numFemDofs, :) = W_fem_mat;
            
            % 4. 计算 GeoBalance 缩放因子 (用于数值调理)
            dt_init = timeSteps(1);
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            
            % 5. 处理边界条件 (仅固定 FEMDoF，电路 DoF 为待求量)
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs <= numFemDofs); 
            
            % 6. 初始化状态向量
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                x_curr(1:length(init_State)) = init_State;
            end
            
            % BDF2 状态序列初始化 (x_prev: t-1, x_prev2: t-2)
            ctx.x_prev = x_curr; 
            ctx.x_prev2 = x_curr; 
            
            % 7. [关键] 预处理 t=0 初值点用于绘图
            I_initial = x_curr(numFemDofs+1 : end); % 初始电流向量
            B_mag_initial = 0;
            postProc = [];
            if nargin >= 13 && ~isempty(probePoint)
                postProc = PostProcessor(obj.Assembler);
                B0_vec = postProc.probeB(x_curr(1:ctx.num_A), probePoint);
                B_mag_initial = norm(B0_vec);
            end
            
            % 在循环前执行第一次回调，绘制 t=0 点
            if ~isempty(monitorHandle)
                if nargin(monitorHandle) >= 6
                    monitorHandle(0, I_initial, 0, I_initial', B_mag_initial, B_mag_initial);
                else
                    monitorHandle(0, I_initial, 0, I_initial');
                end
            end
            
            % 8. 初始化历史容器
            numSteps = length(timeSteps);
            currentHistory = zeros(numSteps, numWindings); 
            probeB_History = zeros(numSteps, 1);
            times = zeros(numSteps, 1);
            currentTime = 0;
            
            % ==================================================
            % 时间步主循环 (Time Stepping Loop)
            % ==================================================
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                times(t_idx) = currentTime;
                
                % --- [BDF 系数选择] ---
                if t_idx == 1
                    % 第一步使用 BDF1 (Backward Euler)
                    alpha_t = 1.0 / dt;
                    beta_t_vec = -ctx.x_prev / dt;
                    fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f, BDF1) ---\n', t_idx, numSteps, dt, currentTime);
                else
                    % 之后使用二阶 BDF2
                    alpha_t = 1.5 / dt;
                    beta_t_vec = (-2.0 * ctx.x_prev + 0.5 * ctx.x_prev2) / dt;
                    fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f, BDF2) ---\n', t_idx, numSteps, dt, currentTime);
                end
                
                ctx.bdf_alpha = alpha_t;
                ctx.bdf_beta_vec = beta_t_vec;
                
                % 获取各绕组当前时刻电压源值
                V_vals = zeros(numWindings, 1);
                for k = 1:numWindings
                    V_vals(k) = ctx.circuit(k).V_source_func(currentTime);
                end
                ctx.V_source_vals = V_vals;
                
                res_norm_init = 0;
                
                % ==============================================
                % Newton-Raphson 非线性迭代
                % ==============================================
                for iter = 1:obj.MaxIterations
                    
                    % [Step A] 组装 Newton 系统 (J, Res)
                    [J_sys, Res_vec] = obj.assembleNewtonSystem(x_curr, ctx);
                    
                    % 计算残差范数 (判断收敛)
                    Res_check = Res_vec;
                    Res_check(ctx.fixedDofs) = 0;
                    current_res_norm = norm(Res_check);
                    
                    if iter == 1
                        res_norm_init = current_res_norm;
                        if res_norm_init < 1e-12, res_norm_init = 1.0; end
                    end
                    rel_res = current_res_norm / res_norm_init;
                    
                    % 打印日志
                    I_vec_now = x_curr(numFemDofs+1 : end);
                    fprintf('    Iter %d: Res=%.4e, RelRes=%.4e, I1=%.3e', iter, current_res_norm, rel_res, I_vec_now(1));
                    
                    if current_res_norm < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        fprintf(' -> Converged.\n');
                        break;
                    end
                    fprintf('\n');
                    
                    % [Step B] 求解线性增量 dx
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res_vec, ctx.fixedDofs);
                    dx = obj.LinearSolver.solve(J_bc, -Res_bc); 
                    
                    % [Step C] 自适应线搜索 (Backtracking Line Search)
                    alpha = 1.0; 
                    if obj.UseLineSearch
                        for k = 0:obj.MaxLineSearchIters
                            if obj.computeResidualNorm(x_curr + alpha * dx, ctx) < current_res_norm * 1.0001
                                break;
                            end
                            alpha = alpha * 0.5;
                        end
                    end
                    
                    % [Step D] 更新解向量
                    x_curr = x_curr + alpha * dx;
                end
                
                % --- [BDF 状态滑动] ---
                ctx.x_prev2 = ctx.x_prev; 
                ctx.x_prev = x_curr; 
                
                % 记录结果
                solutionResults = x_curr;
                I_step = x_curr(numFemDofs+1 : end);
                currentHistory(t_idx, :) = I_step(:)';
                
                % 计算磁密
                B_mag_step = 0;
                if ~isempty(postProc)
                    A_sol = x_curr(1:ctx.num_A);
                    B_vec = postProc.probeB(A_sol, probePoint);
                    B_mag_step = norm(B_vec);
                    probeB_History(t_idx) = B_mag_step;
                end
                
                % 实时绘图 (拼接 t=0 点)
                if ~isempty(monitorHandle)
                    try
                        t_plot = [0; cumsum(timeSteps(1:t_idx))];
                        I_plot = [I_initial'; currentHistory(1:t_idx, :)];
                        
                        if nargin(monitorHandle) >= 6
                            B_plot = [B_mag_initial; probeB_History(1:t_idx)];
                            monitorHandle(currentTime, I_step, t_plot, I_plot, B_mag_step, B_plot);
                        else
                            monitorHandle(currentTime, I_step, t_plot, I_plot);
                        end
                    catch
                    end
                end
            end
            
            info.times = times;
            info.FinalTime = currentTime; 
            info.CurrentHistory = currentHistory;
            if ~isempty(probePoint), info.ProbeB_History = probeB_History; end
        end
    end
    
    methods (Access = protected)
        
        function scaleFactors = calculateCircuitScale(obj, ctx, dt)
            % [GeoBalance Restore] 严格遵循原版采样逻辑
            % 平衡空气(硬刚度)与铁芯(软刚度)的矩阵权重
            fprintf('   [Init] Calculating smart circuit scaling (GeoBalance)...\n');
            
            % 1. 估算硬刚度中位数 (空气/饱和极限)
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            
            % 2. 估算软刚度中位数 (铁芯最大导磁率)
            nu_vec_soft = obj.getBoundNuVec(ctx, 'soft');
            K_soft = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_soft);
            soft_vals = full(abs(diag(K_soft)) + abs(diag_M(1:num_fem)) / dt);
            
            meshTags = obj.Assembler.Mesh.RegionTags;
            matLib = ctx.matLib;
            dofMap = obj.Assembler.DofHandler.DofMaps(ctx.space_A.toString());
            
            target_tags = [];
            keys = matLib.keys;
            for i = 1:length(keys)
                tag = keys{i};
                if ischar(tag), tag_val = str2double(tag); else, tag_val = tag; end
                if strcmpi(matLib(tag).Type, 'Nonlinear'), target_tags = [target_tags, tag_val]; end
            end
            
            if isempty(target_tags)
                k_soft_rep = k_hard_median; 
            else
                mask_elem = ismember(meshTags, target_tags);
                target_dofs = unique(dofMap(:, mask_elem));
                target_dofs = target_dofs(target_dofs > 0);
                k_soft_rep = median(soft_vals(target_dofs));
            end
            
            % 3. 目标权重: 几何平均值
            k_target = sqrt(k_hard_median * k_soft_rep);
            
            % 4. 为每个绕组生成独立 Scale
            numW = ctx.numWindings;
            scaleFactors = zeros(numW, 1);
            for k = 1:numW
                Z_cir = abs(ctx.circuit(k).R + ctx.circuit(k).L/dt); 
                w_col = ctx.W_mat(:, k);
                max_W = max(abs(nonzeros(w_col)));
                max_circuit_val = max(Z_cir, max_W / dt);
                if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
                scaleFactors(k) = k_target / max_circuit_val;
            end
            fprintf('     -> Smart Scale Factors: [%s]\n', sprintf('%.2e ', scaleFactors));
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            tags = obj.Assembler.Mesh.RegionTags;
            nu_vec = zeros(length(tags), 1);
            nu0 = 1 / (4*pi*1e-7);
            uTags = unique(tags);
            for i = 1:length(uTags)
                tag = uTags(i); if ~ctx.matLib.isKey(tag), continue; end
                mat = ctx.matLib(tag);
                if strcmp(mat.Type, 'Linear'), val = mat.Nu_Linear;
                else, if strcmp(mode, 'hard'), val = nu0; else, val = nu0 / 1000.0; end
                end
                nu_vec(tags == tag) = val;
            end
        end
        
        function [J_sys, Res_vec] = assembleNewtonSystem(obj, x, ctx)
            % 组装 Newton 系统 (BDF2 矩阵化版本)
            numFem = ctx.numFemDofs;
            numW = ctx.numWindings;
            
            % [A] FEM 雅可比与残差
            A_vec = x(1:ctx.num_A);
            [J_fem, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            [ji, jj, jv] = find(J_fem);
            J_sys = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            % 时间导数项 dx/dt = alpha*x + beta
            dx_dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec;
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            J_sys = J_sys + ctx.bdf_alpha * ctx.M_sigma;
            
            % 库伦规范约束
            J_sys = J_sys + ctx.C_AP + ctx.C_AP';
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            
            % [B] 场路耦合逻辑
            I_vec = x(numFem+1 : end);
            Scale = ctx.CircuitRowScale; 
            
            % 1. 列耦合 (Field <- Circuit): -W * I (不带缩放)
            Res_vec = Res_vec - ctx.W_mat * I_vec;
            [wi, wj, wv] = find(ctx.W_mat);
            J_sys = J_sys + sparse(wi, numFem + wj, -wv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 2. 行耦合 (Circuit <- Field): Scale * W' * dA/dt
            % Jacobian 块组装
            W_scaled = ctx.W_mat * spdiags(Scale, 0, numW, numW);
            [wsi, wsj, wsv] = find(W_scaled);
            J_sys = J_sys + sparse(numFem + wsj, wsi, ctx.bdf_alpha * wsv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 感应电动势 (全维度矩阵乘法，保证数值稳定性)
            EMF_vec = ctx.W_mat' * dx_dt; 
            
            % [C] 电路对角项逻辑
            R_vals = [ctx.circuit.R]'; L_vals = [ctx.circuit.L]';
            
            % 电路方程残差: Scale * (R*I + L*di/dt + EMF - V)
            term_circuit = R_vals .* I_vec + L_vals .* dx_dt(numFem+1:end) + EMF_vec - ctx.V_source_vals;
            Res_vec(numFem+1 : end) = term_circuit .* Scale;
            
            % 电路 Jacobian 对角块: Scale * (R + L*alpha)
            diag_val = (R_vals + L_vals * ctx.bdf_alpha) .* Scale;
            idx_circuit = (numFem + 1 : ctx.numTotalDofs)';
            J_sys = J_sys + sparse(idx_circuit, idx_circuit, diag_val, ctx.numTotalDofs, ctx.numTotalDofs);
        end
        
        function res_norm = computeResidualNorm(obj, x, ctx)
            % 快速计算残差范数 (用于线搜索)
            numFem = ctx.numFemDofs;
            numW = ctx.numWindings;
            A_vec = x(1:ctx.num_A);
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            dx_dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec;
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            Res_vec = Res_vec - ctx.W_mat * x(numFem+1 : end);
            
            R_vals = [ctx.circuit.R]'; L_vals = [ctx.circuit.L]';
            EMF_vec = ctx.W_mat' * dx_dt;
            term_circuit = R_vals .* x(numFem+1 : end) + L_vals .* dx_dt(numFem+1:end) + EMF_vec - ctx.V_source_vals;
            Res_vec(numFem+1 : end) = term_circuit .* ctx.CircuitRowScale;
            
            Res_vec(ctx.fixedDofs) = 0;
            res_norm = norm(Res_vec);
        end
    end
end