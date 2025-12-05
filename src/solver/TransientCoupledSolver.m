classdef TransientCoupledSolver < handle
    % TRANSIENTCOUPLEDSOLVER 场路耦合非线性瞬态求解器 (v7.12 - Full Commented)
    % 
    % 描述:
    %   本类实现了二维/三维瞬态磁场与外电路的强耦合求解。
    %   采用全耦合 (Monolithic) 策略，将有限元方程 (A, P) 和电路方程 (I) 
    %   组装成一个统一的大型稀疏线性系统进行联立求解。
    %
    % 核心算法特性:
    %   1. 时间积分: Backward Euler (隐式欧拉法)，一阶精度，无条件稳定。
    %   2. 非线性迭代: Newton-Raphson (牛顿-拉夫逊法)，利用切线刚度矩阵实现二次收敛。
    %   3. 规范条件: Coulomb Gauge (库伦规范, div A = 0)，通过引入 Lagrange 乘子 P 
    %               来施加，有效消除 A 场方程在低频/静态下的零空间奇异性。
    %   4. 数值调理: 采用 "GeoBalance" (几何平衡) 策略自动计算电路方程的缩放因子，
    %               平衡铁芯(软)与空气(硬)的刚度差异，防止矩阵病态。
    %   5. 全局收敛: 内置自适应回溯线搜索 (Backtracking Line Search)，防止迭代发散。
    %   6. 实时监控: 支持通过 monitorHandle 回调函数在迭代过程中实时绘制曲线。
    
    properties
        Assembler       % 有限元组装器对象 (负责计算单元刚度、质量矩阵、耦合项等)
        LinearSolver    % 线性方程组求解器对象 (通常封装了 MUMPS 或 MATLAB 的 mldivide)
        
        % --- 迭代控制参数 ---
        MaxIterations = 30        % 最大非线性迭代次数 (通常 5-10 次即可收敛)
        Tolerance = 1e-3          % 绝对残差容差 (收敛判据: ||R|| < Tol)
        RelTolerance = 1e-3       % 相对残差容差 (收敛判据: ||R||/||R_init|| < RelTol)
        
        % --- 线搜索控制参数 ---
        UseLineSearch = true      % 开关：是否启用自适应线搜索 (建议开启以保证稳定性)
        MaxLineSearchIters = 10   % 线搜索最大回溯次数 (例如 1 -> 0.5 -> 0.25 ...)
        MinRelaxation = 1e-4      % 最小允许松弛因子 (防止步长过小导致停滞)
    end
    
    methods
        function obj = TransientCoupledSolver(assembler)
            % 构造函数: 初始化求解器实例
            obj.Assembler = assembler;
            
            % 初始化线性求解器 (MUMPS)
            % 注意: 场路耦合矩阵通常是非对称的 (含有电路耦合项)，且由于
            % 引入了 Lagrange 乘子，矩阵是鞍点型的 (Indefinite，即特征值有正有负)，
            % 因此必须配置求解器为非对称模式 (Symmetry=0) 以确保稳定性。
            obj.LinearSolver = LinearSolver('Auto');
            obj.LinearSolver.MumpsICNTL.i14 = 200; % 增加工作空间百分比 (防止内存不足错误)
            obj.LinearSolver.MumpsSymmetry = 0;    % 强制使用非对称求解策略
        end
        
        function [solutionResults, info] = solve(obj, space_A, space_P, ...
                                                 matLib, sigmaMap, ...
                                                 circuitProps, windingObj, ...
                                                 timeSteps, ...
                                                 fixedDofs_A, fixedDofs_P, ...
                                                 init_State, ...
                                                 monitorHandle)
            % SOLVE 执行瞬态求解主循环
            %
            % 输入参数:
            %   space_A:      磁矢位 A 的有限元空间对象 (通常是 Nedelec 棱边元)
            %   space_P:      Lagrange 乘子 P 的有限元空间对象 (通常是 Nodal 节点元)
            %   matLib:       材料库对象 (包含 B-H 曲线数据)
            %   sigmaMap:     导电率分布映射 (Map<RegionTag, Sigma>)
            %   circuitProps: 电路参数结构体 (包含 R, L, V_source_func)
            %   windingObj:   绕组几何对象 (用于计算耦合系数 W)
            %   timeSteps:    时间步长数组 [dt1, dt2, ...]
            %   fixedDofs_*:  边界条件 (Dirichlet) 的自由度索引
            %   init_State:   初始状态向量 (可选，用于续算)
            %   monitorHandle:(可选) 回调函数句柄 @(t, I, t_hist, I_hist)，用于实时绘图
            
            fprintf('==================================================\n');
            fprintf('   Transient Solver (v7.12 Newton + Monitor)      \n');
            fprintf('==================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % --------------------------------------------------
            % 1. 自由度检查与分配
            % --------------------------------------------------
            % 确保 A 场和 P 场空间都已正确分配全局自由度编号
            if ~dofHandler.DofMaps.isKey(space_A.toString()), dofHandler.distributeDofs(space_A); end
            if ~dofHandler.DofMaps.isKey(space_P.toString()), dofHandler.distributeDofs(space_P); end
            
            % 获取局部自由度数
            ctx.num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            % 计算全局总自由度数: FEM自由度 (A+P) + 1 个电路自由度 (电流 I)
            ctx.numTotalDofs = dofHandler.NumGlobalDofs + 1; 
            
            % 构建上下文结构体 (Context)，用于在各私有函数间传递参数，避免参数列表过长
            ctx.space_A = space_A;
            ctx.space_P = space_P;
            ctx.matLib = matLib;
            ctx.circuit = circuitProps;
            
            % 处理可选参数 monitorHandle
            if nargin < 12, monitorHandle = []; end
            
            % --------------------------------------------------
            % 2. 组装时不变矩阵 (Invariant Matrices)
            % --------------------------------------------------
            % 这些矩阵只与几何和线性材料属性(如电导率)有关，无需在每步非线性迭代中重算，
            % 预先组装可以显著提高效率。
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % [M_sigma] 涡流质量矩阵: int( sigma * N_i * N_j ) dV
            % 对应物理项: sigma * dA/dt (涡流效应)
            M_sigma_local = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
            [mi, mj, mv] = find(M_sigma_local);
            ctx.M_sigma = sparse(mi, mj, mv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % [C_AP] 库伦规范约束矩阵: int( grad(P) * A ) dV
            % 对应方程: div A = 0 的弱形式，通过 Lagrange 乘子 P 施加
            fprintf('   [Init] Assembling Coulomb gauge constraint...\n');
            C_AP_local = obj.Assembler.assembleCoupling(space_A, space_P, 1.0);
            [ci, cj, cv] = find(C_AP_local);
            ctx.C_AP = sparse(ci, cj, cv, ctx.numTotalDofs, ctx.numTotalDofs);
            
            % [W_vec] 绕组耦合向量
            % 描述绕组电流密度分布 J_coil = (N/S) * I * dir
            W_fem = obj.Assembler.assembleWinding(space_A, windingObj);
            ctx.W_vec = sparse(ctx.numTotalDofs, 1);
            ctx.W_vec(1:length(W_fem)) = W_fem;
            
            % --------------------------------------------------
            % 3. 计算自适应电路缩放因子 (Scale Factor)
            % --------------------------------------------------
            % 电路方程 (V=IR) 的系数量级通常远小于 FEM 方程 (K*A)，这会导致矩阵条件数极差。
            % 这里采用 "GeoBalance" 策略，计算一个缩放因子，使电路行与 FEM 行的量级匹配。
            dt_init = timeSteps(1);
            ctx.CircuitRowScale = obj.calculateCircuitScale(ctx, dt_init);
            
            % --------------------------------------------------
            % 4. 合并边界条件索引
            % --------------------------------------------------
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            if islogical(fixedDofs_P), idx_P = find(fixedDofs_P); else, idx_P = fixedDofs_P(:); end
            
            ctx.fixedDofs = [idx_A; idx_P];
            % 确保不误伤电路自由度 (最后一项是电流 I，它是未知量，不应被固定)
            ctx.fixedDofs = ctx.fixedDofs(ctx.fixedDofs < ctx.numTotalDofs); 
            
            % 5. 初始化状态向量
            x_curr = zeros(ctx.numTotalDofs, 1);
            if nargin >= 11 && ~isempty(init_State)
                n = min(length(init_State), ctx.numTotalDofs);
                x_curr(1:n) = init_State(1:n);
            end
            
            numSteps = length(timeSteps);
            solutionResults = []; 
            currentHistory = zeros(numSteps, 1); 
            currentTime = 0;
            
            % ==================================================
            % 6. 时间步循环 (Time Stepping Loop)
            % ==================================================
            for t_idx = 1:numSteps
                dt = timeSteps(t_idx);
                currentTime = currentTime + dt;
                
                fprintf('\n--- Time Step %d / %d (dt=%.1e, t=%.4f) ---\n', t_idx, numSteps, dt, currentTime);
                
                % 更新上下文中的时间相关变量
                ctx.dt = dt;
                ctx.x_prev = x_curr; % 保存上一时间步的解 (用于计算 dx/dt)
                ctx.V_source_val = circuitProps.V_source_func(currentTime); % 计算当前电压源值
                
                res_norm_init = 0;
                
                % ==============================================
                % 7. 牛顿迭代循环 (Newton-Raphson Loop)
                % ==============================================
                for iter = 1:obj.MaxIterations
                    
                    % [Step A] 组装 Newton 系统 (雅可比矩阵 J 和 残差向量 Res)
                    % J_sys: 全局切线刚度矩阵 (Tangent Stiffness Matrix)
                    % Res_vec: 物理残差向量 (Residual Vector)
                    [J_sys, Res_vec] = obj.assembleNewtonSystem(x_curr, ctx);
                    
                    % 计算当前残差范数 (用于收敛判断)
                    Res_check = Res_vec;
                    Res_check(ctx.fixedDofs) = 0; % 排除边界点 (Dirichlet 边界残差恒为0)
                    current_res_norm = norm(Res_check);
                    
                    % 记录初始残差，用于计算相对收敛率
                    if iter == 1
                        res_norm_init = current_res_norm;
                        if res_norm_init < 1e-12, res_norm_init = 1.0; end
                    end
                    rel_res = current_res_norm / res_norm_init;
                    
                    % 输出迭代日志
                    I_val = x_curr(end);
                    fprintf('    Iter %d: Res = %.4e, RelRes = %.4e, I = %.6e', ...
                            iter, current_res_norm, rel_res, full(I_val));
                    
                    % [收敛检查] 满足绝对容差 或 相对容差 即认为收敛
                    if current_res_norm < obj.Tolerance || (iter > 1 && rel_res < obj.RelTolerance)
                        fprintf(' -> Converged.\n');
                        break;
                    end
                    fprintf('\n');
                    
                    % [Step B] 求解线性方程组 J * dx = -Res
                    % 应用 Dirichlet 边界条件 (将对应行/列置零，对角线置1)
                    [J_bc, Res_bc] = BoundaryCondition.applyDirichlet(J_sys, Res_vec, ctx.fixedDofs);
                    
                    lastwarn(''); % 清空警告
                    dx = obj.LinearSolver.solve(J_bc, -Res_bc); % 注意右端项为负残差
                    
                    % 检查矩阵病态 (Singularity Check)
                    [warnMsg, ~] = lastwarn;
                    if contains(warnMsg, 'singular') || contains(warnMsg, 'badly scaled') || any(isnan(dx))
                        fprintf('      [Warn] Matrix ill-conditioned. Adding regularization shift.\n');
                        % 挽救措施：在对角线添加微小位移 (Shift) 以保证矩阵可逆
                        % 这通常发生在材料进入深度饱和或 Lagrange 乘子块零对角线导致数值不稳定时
                        J_fix = J_bc + 1e-8 * speye(size(J_bc));
                        dx = obj.LinearSolver.solve(J_fix, -Res_bc);
                    end
                    
                    % [Step C] 自适应线搜索 (Adaptive Line Search)
                    % 目标: 寻找步长 alpha (0 < alpha <= 1)，使得新解 x_new = x + alpha*dx
                    %       的残差范数显著小于当前残差范数。
                    alpha = 1.0; 
                    if obj.UseLineSearch
                        accepted = false;
                        for k = 0:obj.MaxLineSearchIters
                            x_try = x_curr + alpha * dx;
                            
                            % 快速计算试探步的残差范数 (仅计算向量，不组装 Jacobian 矩阵，速度快)
                            try_res_norm = obj.computeResidualNorm(x_try, ctx);
                            
                            % 判定条件: 残差下降 (允许 0.01% 的数值波动容差)
                            if try_res_norm < current_res_norm * 1.0001
                                if k > 0
                                    fprintf('      [LS] Backtrack: Alpha=%.4f (Res: %.2e -> %.2e)\n', ...
                                        alpha, current_res_norm, try_res_norm);
                                end
                                accepted = true; 
                                break;
                            else
                                % 残差未下降，步长减半回溯
                                alpha = alpha * 0.5;
                            end
                            
                            % 达到最小步长限制，强制更新 (避免迭代死锁)
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
                
                % 记录本时间步的最终解
                solutionResults = x_curr;
                currentHistory(t_idx) = x_curr(end);
                
                % [实时监控] 调用用户传入的句柄进行绘图
                if ~isempty(monitorHandle)
                    try
                        % 构造累积时间轴和电流历史
                        t_vec = cumsum(timeSteps(1:t_idx));
                        I_vec = currentHistory(1:t_idx);
                        
                        % 调用回调函数
                        monitorHandle(currentTime, x_curr(end), t_vec, I_vec);
                        grid on;
                        drawnow limitrate; % 强制刷新图形窗口 (limitrate 限制刷新率以保证性能)
                    catch ME
                        fprintf('Monitor handle execution failed: %s', ME.message);
                    end
                end
            end
            
            info.FinalTime = currentTime; 
            info.CurrentHistory = currentHistory;
        end
    end
    
    % ==================================================
    % 私有辅助方法 (Private Methods)
    % ==================================================
    methods (Access = private)
        
        function scaleFactor = calculateCircuitScale(obj, ctx, dt)
            % CALCULATECIRCUITSCALE 计算自适应电路方程缩放因子
            % 
            % 策略 (Geometric Balance Strategy):
            %   Scale = Target_Stiffness / Max_Circuit_Coeff
            %   Target_Stiffness = sqrt( K_hard * K_soft )
            %
            % 原理:
            %   FEM 矩阵中包含极高刚度(空气, K_hard)和较低刚度(铁芯, K_soft)的区域。
            %   如果电路方程仅匹配空气刚度，会掩盖铁芯的非线性变化；如果仅匹配铁芯，
            %   又可能被空气刚度数值淹没。几何平均值能在这两者之间取得最佳平衡。
            
            fprintf('   [Init] Calculating smart circuit scaling (GeoBalance)...\n');
            
            % 1. 估算 Hard Stiffness (空气/饱和极限)
            % 使用真空磁阻率 nu0 计算全场刚度
            nu_vec_hard = obj.getBoundNuVec(ctx, 'hard'); 
            K_hard = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_hard);
            
            diag_M = diag(ctx.M_sigma);
            num_fem = size(K_hard, 1);
            
            % 计算有效对角元 (包含刚度项 K 和 惯性项 M/dt)
            hard_vals = abs(diag(K_hard)) + abs(diag_M(1:num_fem)) / dt;
            % 使用中位数代表"环境/空气"的刚度基准
            k_hard_median = full(median(hard_vals(hard_vals > 1e-12)));
            
            % 2. 估算 Soft Stiffness (铁芯最大导磁率极限)
            % 使用 nu0 / mu_r_max 计算全场刚度
            nu_vec_soft = obj.getBoundNuVec(ctx, 'soft');
            K_soft = obj.Assembler.assembleStiffness(ctx.space_A, nu_vec_soft);
            
            soft_vals = full(abs(diag(K_soft)) + abs(diag_M(1:num_fem)) / dt);
            
            % 筛选非线性材料区域 (铁芯)
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
            
            % 如果存在铁芯，提取铁芯区域的中位数作为软基准；否则回退到全场基准
            if isempty(target_tags)
                k_soft_rep = k_hard_median; 
            else
                mask_elem = ismember(meshTags, target_tags);
                target_dofs = unique(dofMap(:, mask_elem));
                target_dofs = target_dofs(target_dofs > 0);
                k_soft_rep = median(soft_vals(target_dofs));
            end
            
            % 3. 计算几何平均目标刚度
            k_target = sqrt(k_hard_median * k_soft_rep);
            
            fprintf('     -> Stiffness Bounds: Hard(Air)=%.2e, Soft(Iron)=%.2e\n', ...
                    k_hard_median, k_soft_rep);
            fprintf('     -> Balanced Target (GeoMean): %.2e\n', k_target);
            
            % 4. 估算电路行的最大系数
            Z_cir = abs(ctx.circuit.R + ctx.circuit.L/dt); % 电路阻抗
            max_W = max(abs(nonzeros(ctx.W_vec)));         % 绕组耦合系数
            max_coupling_val = max_W / dt;                 % 耦合项的时间导数系数
            
            % 取两者较大者作为分母
            max_circuit_val = max(Z_cir, max_coupling_val);
            if max_circuit_val < 1e-20, max_circuit_val = 1.0; end
            
            % 5. 计算最终缩放因子
            scaleFactor = k_target / max_circuit_val;
            
            fprintf('     -> Max Circuit Coeff : %.2e (Z=%.1e, W/dt=%.1e)\n', ...
                    max_circuit_val, Z_cir, max_coupling_val);
            fprintf('     -> Smart Scale Factor: %.2e\n', scaleFactor);
        end
        
        function nu_vec = getBoundNuVec(obj, ctx, mode)
            % 辅助函数: 获取极限状态下的磁阻率向量 nu
            meshTags = obj.Assembler.Mesh.RegionTags;
            numElems = length(meshTags);
            nu_vec = zeros(numElems, 1);
            matLib = ctx.matLib;
            nu0 = 1 / (4 * pi * 1e-7); % 真空磁阻率
            
            uniqueTags = unique(meshTags);
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k);
                if ~matLib.isKey(tag), continue; end
                mat = matLib(tag);
                
                val = 0;
                if strcmp(mat.Type, 'Linear')
                    val = mat.Nu_Linear; % 线性材料属性不变
                else
                    if strcmp(mode, 'hard')
                        val = nu0; % 硬极限: 完全饱和 (空气磁阻率)
                    else
                        % 软极限: 极高导磁率状态 (nu0 / mu_r_max)
                        % 通过遍历 B-H 曲线找到最大导磁率点
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
        
        function [J_sys, Res_vec] = assembleNewtonSystem(obj, x, ctx)
            % 组装完整的 Newton-Raphson 线性系统
            % 目标方程: J_sys * dx = -Res_vec
            
            % 1. FEM 部分 (A场)
            A_vec = x(1:ctx.num_A);
            % assembleJacobian 调用内核计算切线刚度 J_fem 和物理残差 R_fem
            [J_fem, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, true);
            
            % 扩展矩阵大小到总自由度
            [ji, jj, jv] = find(J_fem);
            J_sys = sparse(ji, jj, jv, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            % 2. 瞬态项 (M/dt) 和 规范项 (Lagrange)
            dx_dt = (x - ctx.x_prev) / ctx.dt;
            
            % 涡流残差: + M * (x - x_prev)/dt
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            % 涡流雅可比: + M / dt
            J_sys = J_sys + ctx.M_sigma / ctx.dt;
            
            % 规范残差: + (C_AP + C_AP') * x  (即: div A = 0, grad P = 0)
            J_sys = J_sys + ctx.C_AP + ctx.C_AP';
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            
            % 3. 电路耦合 (应用 Scale 缩放)
            I_val = x(end);
            Scale = ctx.CircuitRowScale; 
            [w_idx, ~, w_val] = find(ctx.W_vec);
            
            % 列耦合: -W * I (物理方程中的电流源项，通常不缩放)
            J_col = sparse(w_idx, repmat(ctx.numTotalDofs, length(w_idx), 1), -w_val, ctx.numTotalDofs, ctx.numTotalDofs);
            Res_vec = Res_vec - ctx.W_vec * I_val;
            
            % 行耦合: (W'/dt * A) * Scale (电路方程中的感应电动势项，需缩放)
            J_row_entries = (w_val / ctx.dt) * Scale;
            J_row = sparse(repmat(ctx.numTotalDofs, length(w_idx), 1), w_idx, J_row_entries, ctx.numTotalDofs, ctx.numTotalDofs);
            J_sys = J_sys + J_col + J_row;
            
            % 4. 电路对角项与残差
            term_R = ctx.circuit.R * I_val;
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt;
            
            % 电路残差: (R*I + L*dI/dt + EMF - V_source) * Scale
            res_I_val = (term_R + term_L + term_EMF - ctx.V_source_val) * Scale;
            Res_vec(end) = res_I_val;
            
            % 电路雅可比对角: (R + L/dt) * Scale
            Z_circuit = (ctx.circuit.R + ctx.circuit.L / ctx.dt) * Scale;
            J_sys(end, end) = J_sys(end, end) + Z_circuit;
        end
        
        function res_norm = computeResidualNorm(obj, x, ctx)
            % 快速计算残差范数 (用于线搜索，跳过 Jacobian 组装)
            
            A_vec = x(1:ctx.num_A);
            % assembleJacobian(..., false) 仅返回残差向量 R_fem
            [~, R_fem] = obj.Assembler.assembleJacobian(ctx.space_A, A_vec, ctx.matLib, false);
            
            Res_vec = sparse(ctx.numTotalDofs, 1);
            Res_vec(1:length(R_fem)) = R_fem;
            
            dx_dt = (x - ctx.x_prev) / ctx.dt;
            Res_vec = Res_vec + ctx.M_sigma * dx_dt;
            Res_vec = Res_vec + (ctx.C_AP + ctx.C_AP') * x;
            Res_vec = Res_vec - ctx.W_vec * x(end);
            
            Scale = ctx.CircuitRowScale;
            term_R = ctx.circuit.R * x(end);
            term_L = ctx.circuit.L * dx_dt(end);
            term_EMF = ctx.W_vec' * dx_dt;
            res_I_val = (term_R + term_L + term_EMF - ctx.V_source_val) * Scale;
            Res_vec(end) = res_I_val;
            
            % 消除 Dirichlet 边界上的伪残差 (因为该处解被固定，不应计入误差)
            Res_vec(ctx.fixedDofs) = 0;
            res_norm = norm(Res_vec);
        end
    end
end