classdef BDF2BreakerSolver < TransientBDF2Solver
    % BDF2BREAKERSOLVER 带分闸断路器逻辑 + 能量审计 + 完整监测功能的求解器
    % 
    % 功能描述:
    %   1. Breaker Logic: t > BreakerTime 后，强行置 I=0，模拟理想分闸。
    %   2. Energy Audit: 实时计算数值耗散功率 P_num，评估算法保真度。
    %   3. Full Monitoring: 支持实时绘制电流、电压、磁密波形 (同父类功能)。
    
    properties
        BreakerTime = 0.01       % [s] 分闸时刻
        R_winding_Audit = 1.5    % [Ohm] 用于能量审计的绕组电阻
    end
    
    properties (Access = protected)
        % --- 能量审计历史变量 ---
        Hist_W_mag = 0           % 当前磁场储能 [J]
        Hist_E_ohm_cum = 0       % 累计电阻损耗 [J]
        Hist_E_num_cum = 0       % 累计数值耗散 [J]
        Hist_E_eddy_cum = 0      % 累计涡流损耗 [J]
        Lambda_prev = []         % 上一步的磁链向量
        I_prev_audit = []        % 上一步的电流向量
    end
    
    methods
        function obj = BDF2BreakerSolver(assembler)
            obj@TransientBDF2Solver(assembler);
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
            fprintf('   Breaker Solver (BDF2) - PostProc Integrated    \n');
            fprintf('==================================================\n');
            
            % -----------------------------------------------------------
            % 1. 初始化与矩阵组装
            % -----------------------------------------------------------
            dofHandler = obj.Assembler.DofHandler;
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            numWindings = length(windingObj);
            numFemDofs = dofHandler.NumGlobalDofs; 
            
            % Context 构建
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            ctx.numFemDofs = numFemDofs;
            ctx.numWindings = numWindings;
            ctx.numTotalDofs = numFemDofs + numWindings; 
            ctx.space_A = space_A;
            ctx.space_P = space_P;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            % 矩阵组装
            fprintf('   [Init] Assembling invariant matrices...\n');
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            W_fem_mat = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_mat = sparse(ctx.numTotalDofs, numWindings);
            n_rows_W = size(W_fem_mat, 1);
            ctx.W_mat(1:n_rows_W, :) = W_fem_mat;
            
            % Scaling
            dt_init = timeSteps(1);
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            fprintf('   [Init] Circuit Scale Factor: %.4e\n', ctx.CircuitRowScale(1));
            
            % Boundary Conditions
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            ctx.fixedDofs = [idx_A; idx_P];
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs <= numFemDofs); 
            
            % Initial State
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State), x_curr(1:length(init_State)) = init_State; end
            
            ctx.x_prev = x_curr; 
            ctx.x_prev2 = x_curr; 
            
            % -----------------------------------------------------------
            % 2. 监测器与后处理对象初始化
            % -----------------------------------------------------------
            I_initial = x_curr(numFemDofs+1 : end); 
            
            % 初始化 PostProcessor (复用 Assembler)
            postProc = PostProcessor(obj.Assembler); 
            
            B_mag_initial = 0;
            if nargin >= 13 && ~isempty(probePoint)
                B0_vec = postProc.probeB(x_curr(1:ctx.num_A), probePoint);
                B_mag_initial = norm(B0_vec);
            end
            
            if ~isempty(monitorHandle)
                try
                    if nargin(monitorHandle) >= 6
                        monitorHandle(0, I_initial, 0, I_initial', B_mag_initial, B_mag_initial);
                    else
                        monitorHandle(0, I_initial, 0, I_initial');
                    end
                catch
                end
            end
            
            % 计算 t=0 时刻的精确磁能 (利用 PostProcessor)
            % 注意: PostProcessor 需要纯 A 场向量
            A_init = x_curr(1:ctx.num_A);
            W_mag_prev_exact = postProc.computeMagneticEnergy(A_init, matLib);
            
            % -----------------------------------------------------------
            % 3. 历史记录与能量审计初始化
            % -----------------------------------------------------------
            numSteps = length(timeSteps);
            currentHistory = zeros(numSteps, numWindings); 
            probeB_History = zeros(numSteps, 1); 
            
            energyLog.time = zeros(numSteps, 1);
            energyLog.W_mag = zeros(numSteps, 1);
            energyLog.P_num = zeros(numSteps, 1);
            energyLog.P_eddy = zeros(numSteps, 1); 
            energyLog.LossRatio = zeros(numSteps, 1);
            
            % 复位累积状态
            obj.I_prev_audit = zeros(numWindings, 1);
            obj.Hist_W_mag = W_mag_prev_exact; 
            obj.Hist_E_ohm_cum = 0; 
            obj.Hist_E_num_cum = 0;
            obj.Hist_E_eddy_cum = 0; 
            
            currentTime = 0;
            
            % ==================================================
            % 4. 时间步主循环
            % ==================================================
            fprintf('\nStarting Time Integration (%d steps)...\n', numSteps);
            
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                % Breaker Logic: > (Time + 1e-9)
                isBreakerOpen = currentTime > (obj.BreakerTime + 1e-9);
                breakerStatus = "CLOSED"; if isBreakerOpen, breakerStatus = "OPEN"; end
                
                fprintf('\n--- Step %d/%d | t=%.5fs | dt=%.1e | Breaker: %s ---\n', ...
                        t_idx, numSteps, currentTime, dt, breakerStatus);
                
                % BDF coefficients
                if t_idx == 1
                    alpha_t = 1.0 / dt;
                    beta_t_vec = -ctx.x_prev / dt;
                else
                    alpha_t = 1.5 / dt;
                    beta_t_vec = (-2.0 * ctx.x_prev + 0.5 * ctx.x_prev2) / dt;
                end
                ctx.bdf_alpha = alpha_t;
                ctx.bdf_beta_vec = beta_t_vec;
                
                % Voltage Source
                V_vals = zeros(numWindings, 1);
                if ~isBreakerOpen
                    for k = 1:numWindings
                        V_vals(k) = ctx.circuit(k).V_source_func(currentTime);
                    end
                end
                ctx.V_source_vals = V_vals;
                
                % --- Newton Iteration ---
                fprintf('   %-5s | %-12s | %-12s | %-12s\n', 'Iter', 'Abs Res', 'Rel Res', 'I_coil (A)');
                fprintf('   --------------------------------------------------------\n');
                res_norm_init = 1.0;
                
                for iter = 1:obj.MaxIterations
                    [J_sys, Res_vec] = obj.assembleNewtonSystem(x_curr, ctx);
                    
                    % Apply Breaker Constraints
                    if isBreakerOpen
                        idx_circuit = (numFemDofs + 1 : ctx.numTotalDofs)';
                        J_sys(idx_circuit, :) = 0;
                        Scale = ctx.CircuitRowScale;
                        ind_diag = sub2ind(size(J_sys), idx_circuit, idx_circuit);
                        J_sys(ind_diag) = Scale; 
                        Res_vec(idx_circuit) = x_curr(idx_circuit) .* Scale;
                    end
                    
                    Res_check = Res_vec; Res_check(ctx.fixedDofs) = 0;
                    current_res_norm = norm(Res_check);
                    if iter == 1, res_norm_init = max(current_res_norm, 1e-12); end
                    rel_res = current_res_norm / res_norm_init;
                    
                    fprintf('   %5d | %12.4e | %12.4e | %12.4e\n', iter, current_res_norm, rel_res, x_curr(numFemDofs+1));
                    
                    if current_res_norm < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        fprintf('   >> Converged.\n');
                        break;
                    end
                    
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res_vec, ctx.fixedDofs);
                    dx = obj.LinearSolver.solve(J_bc, -Res_bc);
                    
                    % Line Search
                    alpha = 1.0; 
                    if obj.UseLineSearch
                        for k = 0:obj.MaxLineSearchIters
                            x_try = x_curr + alpha * dx;
                            [~, Res_try] = obj.assembleNewtonSystem(x_try, ctx);
                            if isBreakerOpen
                                idx_c = (numFemDofs + 1 : ctx.numTotalDofs)';
                                Res_try(idx_c) = x_try(idx_c) .* ctx.CircuitRowScale;
                            end
                            Res_try(ctx.fixedDofs) = 0;
                            if norm(Res_try) < current_res_norm * 1.001, break; end
                            alpha = alpha * 0.5;
                        end
                        if alpha < 1.0, fprintf('      [LineSearch] Backtrack: alpha = %.4f\n', alpha); end
                    end
                    x_curr = x_curr + alpha * dx;
                end
                
                % Update State
                ctx.x_prev2 = ctx.x_prev; 
                ctx.x_prev = x_curr;
                
                % Post-Process Data
                A_sol = x_curr(1:ctx.num_A);
                I_step = x_curr(numFemDofs+1 : end);
                currentHistory(t_idx, :) = I_step(:)';
                
                if ~isempty(postProc)
                    B_vec = postProc.probeB(A_sol, probePoint);
                    probeB_History(t_idx) = norm(B_vec);
                end
                
                % =======================================================
                % 精确能量审计 (直接调用 PostProcessor)
                % =======================================================
                
                % 1. 调用 PostProcessor 计算全场磁能
                %    (内部会自动处理非线性积分)
                W_mag_current_exact = postProc.computeMagneticEnergy(A_sol, matLib);
                
                if t_idx > 1
                    % 1.1 磁场储能变化
                    d_W_mag = W_mag_current_exact - W_mag_prev_exact;
                    obj.Hist_W_mag = W_mag_current_exact;
                    
                    % 2. 线圈电阻损耗
                    d_E_ohm = sum((I_step.^2) * obj.R_winding_Audit * dt);
                    obj.Hist_E_ohm_cum = obj.Hist_E_ohm_cum + d_E_ohm;
                    
                    % 3. 涡流损耗 (直接利用质量矩阵计算)
                    %    P_eddy = (dA/dt)' * M_sigma * (dA/dt)
                    dx_dt = ctx.bdf_alpha * x_curr + ctx.bdf_beta_vec;
                    P_eddy_inst = dx_dt' * ctx.M_sigma * dx_dt;
                    d_E_eddy = P_eddy_inst * dt;
                    obj.Hist_E_eddy_cum = obj.Hist_E_eddy_cum + d_E_eddy;
                    energyLog.P_eddy(t_idx) = P_eddy_inst;
                    
                    % 4. 源输入
                    d_E_source = 0;
                    if ~isBreakerOpen
                        d_E_source = sum(V_vals .* I_step * dt); 
                    end
                    
                    % 5. 数值耗散
                    d_E_num = d_E_source - (d_W_mag + d_E_ohm + d_E_eddy);
                    obj.Hist_E_num_cum = obj.Hist_E_num_cum + d_E_num;
                    
                    energyLog.P_num(t_idx) = d_E_num / dt;
                    
                    total_phys = obj.Hist_E_ohm_cum + obj.Hist_E_eddy_cum;
                    if abs(total_phys) > 1e-9
                        energyLog.LossRatio(t_idx) = (obj.Hist_E_num_cum / total_phys) * 100;
                    end
                end
                
                W_mag_prev_exact = W_mag_current_exact;
                obj.I_prev_audit = I_step;
                
                energyLog.time(t_idx) = currentTime;
                energyLog.W_mag(t_idx) = obj.Hist_W_mag;
                
                if ~isempty(monitorHandle)
                    try
                         t_plt = [0; energyLog.time(1:t_idx)];
                         I_plt = [I_initial'; currentHistory(1:t_idx, :)];
                         if nargin(monitorHandle)>=6
                             monitorHandle(currentTime, I_step, t_plt, I_plt, 0, [B_mag_initial; probeB_History(1:t_idx)]);
                         else
                             monitorHandle(currentTime, I_step, t_plt, I_plt);
                         end
                    catch
                    end
                end
            end
            
            solutionResults = x_curr;
            info.times = energyLog.time;
            info.CurrentHistory = currentHistory;
            if ~isempty(postProc), info.ProbeB_History = probeB_History; end
            info.EnergyLog = energyLog; 
            
            obj.plotEnergyAudit(energyLog);
        end
        
        function plotEnergyAudit(~, log)
            figure('Name', 'BDF2 Breaker Energy Audit (Full Physics)', 'NumberTitle', 'off');
            
            ax1 = subplot(3,1,1);
            plot(log.time, log.W_mag, 'LineWidth', 1.5);
            xline(0.01, 'k--', 'Breaker');
            ylabel('Energy (J)'); title('Magnetic Energy W_{mag}'); grid on;
            
            ax2 = subplot(3,1,2);
            % 对比 物理损耗 vs 数值误差
            plot(log.time, log.P_eddy, 'g', 'LineWidth', 1.2, 'DisplayName', 'Physical Eddy Loss');
            hold on;
            plot(log.time, log.P_num, 'r', 'LineWidth', 1.0, 'DisplayName', 'Numerical Error');
            yline(0, 'k-');
            xline(0.01, 'k--', 'Breaker');
            ylabel('Power (W)'); title('Power Balance: Physical vs Numerical'); 
            legend('show'); grid on;
            
            ax3 = subplot(3,1,3);
            plot(log.time, log.LossRatio, 'b', 'LineWidth', 1.5);
            xline(0.01, 'k--', 'Breaker');
            ylabel('Ratio (%)'); title('Cumulative Numerical Error Ratio'); 
            xlabel('Time (s)'); grid on;
            
            linkaxes([ax1, ax2, ax3], 'x');
        end
    end
end