classdef TransientBDF1Solver < handle
    % TRANSIENTBDF1SOLVER 场路耦合非线性瞬态求解器 (多绕组增强精度版)
    % 
    % 算法特性:
    %   - 时间离散: 一阶后向欧拉法 (Backward Euler / BDF1)，L-稳定，抗数值振荡。
    %   - 非线性迭代: Newton-Raphson 法配合自适应线搜索，处理铁芯饱和。
    %   - 场路耦合: 全耦合 (Monolithic) 矩阵求解，严格保证场-路能量平衡。
    %   - 数值调理: 恢复严格 GeoBalance 逻辑，平衡空气与铁芯区域的矩阵条件数。
    %   - 绘图增强: 支持从 t=0 初值点开始实时绘制波形。
    
    properties
        Assembler       % 有限元组装器
        LinearSolver    % 线性直接求解器 (MUMPS)
        
        % --- 迭代控制参数 ---
        MaxIterations = 30        
        Tolerance = 1e-3          
        RelTolerance = 1e-3       
        
        % --- 线搜索控制参数 ---
        UseLineSearch = true      
        MaxLineSearchIters = 10   
        MinRelaxation = 1e-4      
    end
    
    methods
        function obj = TransientBDF1Solver(assembler)
            obj.Assembler = assembler;
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
            fprintf('   Transient Solver (BDF1) - Multi-Winding Strict \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % 1. 自由度分配
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            % 2. 维度与上下文初始化
            numWindings = length(windingObj);
            numFemDofs = dofHandler.NumGlobalDofs; 
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numFemDofs = numFemDofs;
            ctx.numWindings = numWindings;
            ctx.numTotalDofs = numFemDofs + numWindings; 
            
            ctx.space_A = space_A;
            ctx.space_P = space_P;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            % 3. 组装时不变矩阵
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            W_fem_mat = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_mat = sparse(ctx.numTotalDofs, numWindings);
            ctx.W_mat(1:numFemDofs, :) = W_fem_mat;
            
            % 4. GeoBalance 缩放因子计算
            dt_init = timeSteps(1);
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            
            % 5. 边界条件
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs <= numFemDofs); 
            
            % 6. 状态向量初值设置
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                x_curr(1:length(init_State)) = init_State;
            end
            ctx.x_prev = x_curr; 
            
            % 7. 预处理初值点用于绘图 (t=0)
            I_initial = x_curr(numFemDofs+1 : end); % 初始电流向量
            B_mag_initial = 0;
            postProc = [];
            if nargin >= 13 && ~isempty(probePoint)
                postProc = PostProcessor(obj.Assembler);
                B0_vec = postProc.probeB(x_curr(1:ctx.num_A), probePoint);
                B_mag_initial = norm(B0_vec);
            end
            
            % 如果提供了 monitorHandle，先绘制 t=0 的点
            if ~isempty(monitorHandle)
                if nargin(monitorHandle) >= 6
                    monitorHandle(0, I_initial, 0, I_initial', B_mag_initial, B_mag_initial);
                else
                    monitorHandle(0, I_initial, 0, I_initial');
                end
            end
            
            % 8. 循环容器初始化
            numSteps = length(timeSteps);
            currentHistory = zeros(numSteps, numWindings); 
            probeB_History = zeros(numSteps, 1);
            currentTime = 0;
            
            % ==================================================
            % 时间步主循环
            % ==================================================
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                % 更新步进上下文
                ctx.x_prev = x_curr; 
                ctx.dt = dt;
                ctx.bdf_alpha = 1.0 / dt;
                ctx.bdf_beta_vec = -ctx.x_prev / dt;
                
                % 电压源赋值
                V_vals = zeros(numWindings, 1);
                for k = 1:numWindings
                    V_vals(k) = ctx.circuit(k).V_source_func(currentTime);
                end
                ctx.V_source_vals = V_vals;
                
                fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f) ---\n', t_idx, numSteps, dt, currentTime);
                
                res_norm_init = 0;
                % Newton-Raphson 迭代循环
                for iter = 1:obj.MaxIterations
                    [J_sys, Res_vec] = obj.assembleNewtonSystem(x_curr, ctx);
                    
                    Res_check = Res_vec;
                    Res_check(ctx.fixedDofs) = 0;
                    current_res_norm = norm(Res_check);
                    
                    if iter == 1
                        res_norm_init = current_res_norm;
                        if res_norm_init < 1e-12, res_norm_init = 1.0; end
                    end
                    rel_res = current_res_norm / res_norm_init;
                    
                    I_vec_now = x_curr(numFemDofs+1 : end);
                    fprintf('    Iter %d: Res=%.4e, RelRes=%.4e, I1=%.3e', iter, current_res_norm, rel_res, I_vec_now(1));
                    
                    if current_res_norm < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        fprintf(' -> Converged.\n');
                        break;
                    end
                    fprintf('\n');
                    
                    % 线性求解
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res_vec, ctx.fixedDofs);
                    dx = obj.LinearSolver.solve(J_bc, -Res_bc); 
                    
                    % 自适应线搜索
                    alpha = 1.0;
                    if obj.UseLineSearch
                        for k = 0:obj.MaxLineSearchIters
                            if obj.computeResidualNorm(x_curr + alpha * dx, ctx) < current_res_norm * 1.0001
                                break;
                            end
                            alpha = alpha * 0.5;
                        end
                    end
                    x_curr = x_curr + alpha * dx;
                end
                
                % 保存当前步结果
                solutionResults = x_curr;
                I_step = x_curr(numFemDofs+1 : end);
                currentHistory(t_idx, :) = I_step(:)';
                
                % 计算当前步磁密
                B_mag_step = 0;
                if ~isempty(postProc)
                    B_vec = postProc.probeB(x_curr(1:ctx.num_A), probePoint);
                    B_mag_step = norm(B_vec);
                    probeB_History(t_idx) = B_mag_step;
                end
                
                % 实时监控绘制 (包含 t=0 点)
                if ~isempty(monitorHandle)
                    % 构造带初值的时间轴和数据矩阵
                    t_plot = [0; cumsum(timeSteps(1:t_idx))];
                    I_plot = [I_initial'; currentHistory(1:t_idx, :)];
                    
                    if nargin(monitorHandle) >= 6
                        B_plot = [B_mag_initial; probeB_History(1:t_idx)];
                        monitorHandle(currentTime, I_step, t_plot, I_plot, B_mag_step, B_plot);
                    else
                        monitorHandle(currentTime, I_step, t_plot, I_plot);
                    end
                end
            end
            
            info.FinalTime = currentTime; 
            info.CurrentHistory = currentHistory;
            if ~isempty(probePoint), info.ProbeB_History = probeB_History; end
        end
    end
    
    methods (Access = private)
        
        function scaleFactors = calculateCircuitScale(obj, ctx, dt)
            % [GeoBalance Restore] 严格遵循原版采样逻辑
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
            
            k_target = sqrt(k_hard_median * k_soft_rep);
            
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
            numFem = ctx.numFemDofs;
            numW = ctx.numWindings;
            
            % FEM 项
            A_vec = x(1:ctx.num_A);
            [J_fem, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            [ji, jj, jv] = find(J_fem);
            J_sys = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            % 时间导数算子应用
            dx_dt = ctx.bdf_alpha * x + ctx.bdf_beta_vec;
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            J_sys = J_sys + ctx.bdf_alpha * ctx.M_sigma;
            
            % 规范约束
            J_sys = J_sys + ctx.C_AP + ctx.C_AP';
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            
            % 场路耦合逻辑
            I_vec = x(numFem+1 : end);
            Scale = ctx.CircuitRowScale; 
            
            % 磁场方程中的电流源 (不缩放)
            Res_vec = Res_vec - ctx.W_mat * I_vec;
            [wi, wj, wv] = find(ctx.W_mat);
            J_sys = J_sys + sparse(wi, numFem + wj, -wv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 电路方程中的感应电动势 Jacobian (带缩放)
            W_scaled = ctx.W_mat * spdiags(Scale, 0, numW, numW);
            [wsi, wsj, wsv] = find(W_scaled);
            J_sys = J_sys + sparse(numFem + wsj, wsi, ctx.bdf_alpha * wsv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % 物理 EMF 计算 (全维度矩阵乘法，保证精度)
            EMF_vec = ctx.W_mat' * dx_dt; 
            
            % 电路参数准备
            R_vals = [ctx.circuit.R]'; L_vals = [ctx.circuit.L]';
            
            % 电路残差: Scale * (R*I + L*di/dt + EMF - V)
            term_circuit = R_vals .* I_vec + L_vals .* dx_dt(numFem+1:end) + EMF_vec - ctx.V_source_vals;
            Res_vec(numFem+1 : end) = term_circuit .* Scale;
            
            % 电路雅可比对角块
            diag_val = (R_vals + L_vals * ctx.bdf_alpha) .* Scale;
            J_sys = J_sys + sparse(numFem+1:ctx.numTotalDofs, numFem+1:ctx.numTotalDofs, diag_val, ctx.numTotalDofs, ctx.numTotalDofs);
        end
        
        function res_norm = computeResidualNorm(obj, x, ctx)
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