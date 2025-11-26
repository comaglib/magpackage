function [x, info] = solve_nr(Model, b_ext, x0)
% SOLVE_NR 牛顿-拉夫逊非线性求解器 (修正版)
% 
% 求解方程: K(A) * A = b
% 迭代格式: J * dx = -Res
%
% 输入:
%   Model  - 模型结构体 (需包含 Materials, Solver 配置)
%   b_ext  - 外部载荷向量 (N_edges x 1)
%   x0     - 初值 (可选). 如果未提供，默认为 0.
%
% 输出:
%   x      - 收敛解
%   info   - 收敛历史信息
% 
% 改进:
% 1. 收敛检查时屏蔽 Dirichlet 边界的约束反力，只检查自由自由度
% 2. 引入自适应阻尼 (Damping) 防止饱和区震荡

    % fprintf('==============================================\n');
    % fprintf('   Nonlinear Solver (Newton-Raphson)          \n');
    % fprintf('==============================================\n');

    % 1. 配置参数
    SolverConf = Model.Solver.Nonlinear;
    MaxIter = SolverConf.MaxIter;
    Tol = SolverConf.Tolerance;
    
    GaugeMethod = 'Penalty'; 
    if isfield(Model.Solver, 'Gauge'), GaugeMethod = Model.Solver.Gauge.Method; end
    
    numEdges = size(Model.Mesh.Edges, 2);
    numNodes = size(Model.Mesh.P, 2);
    
    use_lagrange = strcmpi(GaugeMethod, 'Lagrange');
    sys_size = numEdges + (use_lagrange * numNodes);
    
    % 初始化解
    if nargin < 3 || isempty(x0)
        x = zeros(sys_size, 1);
    else
        x = x0;
    end
    
    % 2. 准备扩增矩阵 (G_mat) 和 边界条件
    if use_lagrange
        G_mat = assemble_coupling_matrix(Model); 
        if isfield(Model.Runtime, 'FixedEdges')
            fixed_edges = Model.Runtime.FixedEdges;
            fixed_vals = Model.Runtime.FixedVals;
        else
            fixed_edges = []; fixed_vals = [];
        end
        % Lambda BC
        edge_defs = Model.Mesh.Edges(:, fixed_edges);
        bnd_nodes = unique(edge_defs(:));
        lambda_dofs = numEdges + bnd_nodes;
        
        all_fixed_dofs = [fixed_edges(:); lambda_dofs(:)];
        all_fixed_vals = [fixed_vals(:); zeros(length(bnd_nodes), 1)];
        
    else % Penalty
        G_mat = assemble_mass_matrix(Model, 'All');
        if isfield(Model.Runtime, 'FixedEdges')
            all_fixed_dofs = Model.Runtime.FixedEdges;
            all_fixed_vals = Model.Runtime.FixedVals;
        else
            all_fixed_dofs = []; all_fixed_vals = [];
        end
    end
    
    % [关键] 构建自由度掩码 (用于收敛检查)
    % 我们只关心非边界点的残差是否趋于0
    is_free_dof = true(sys_size, 1);
    if ~isempty(all_fixed_dofs)
        is_free_dof(all_fixed_dofs) = false;
    end
    
    % 3. 迭代主循环
    fprintf('  [NR] Starting iterations (Max: %d, Tol: %.1e)\n', MaxIter, Tol);
    
    conv_history = [];
    damping = 1.0;          % 初始阻尼
    prev_res_norm = inf;
    
    for iter = 1:MaxIter
        
        % 3.1 拆分变量
        A_curr = x(1:numEdges);
        lambda_curr = [];
        if use_lagrange, lambda_curr = x(numEdges+1:end); end
        
        % 3.2 组装非线性部分
        [K_tan_mag, Res_mag] = assemble_nonlinear_jacobian(Model, A_curr, b_ext);
        
        % 3.3 组装扩增系统
        if use_lagrange
            C = G_mat;
            Res_global = [Res_mag + C * lambda_curr; C' * A_curr];
            K_global = [K_tan_mag, C; C', sparse(numNodes, numNodes)];
        else
            M = G_mat;
            ref_scale = mean(abs(nonzeros(K_tan_mag)));
            alpha = ref_scale * 1e-6;
            K_global = K_tan_mag + alpha * M;
            Res_global = Res_mag + alpha * M * A_curr;
        end
        
        % 3.4 检查收敛性 (仅检查自由 DOFs)
        % 边界上的残差是反力，通常很大，必须忽略
        res_free = Res_global(is_free_dof);
        res_norm = norm(res_free);
        
        % 相对残差 (归一化)
        norm_b = norm(b_ext);
        if norm_b < 1e-10, norm_b = 1.0; end
        rel_res = res_norm / norm_b;
        
        fprintf('    Iter %2d: Res=%.4e (Rel=%.4e) Damp=%.2f', iter, res_norm, rel_res, damping);
        conv_history(end+1) = rel_res; %#ok<AGROW>
        
        if rel_res < Tol
            fprintf(' -> Converged!\n');
            break;
        end
        
        % 3.5 阻尼调整策略 (Adaptive Damping)
        % 如果残差反弹，说明步长太大，减小阻尼
        % 如果残差下降，尝试恢复阻尼
        if iter > 1
            if res_norm > prev_res_norm
                damping = max(damping * 0.5, 0.1); % 减速
                fprintf(' [Oscillation detected, reducing step]');
            else
                damping = min(damping * 1.2, 1.0); % 加速
            end
        end
        prev_res_norm = res_norm;
        fprintf('\n');
        
        % 3.6 处理边界条件 (Newton Update)
        RHS_step = -Res_global;
        
        if ~isempty(all_fixed_dofs)
            curr_vals = x(all_fixed_dofs);
            target_dx = all_fixed_vals - curr_vals;
            [K_step, RHS_step] = apply_dirichlet_bc(K_global, RHS_step, all_fixed_dofs, target_dx);
        else
            K_step = K_global;
        end
        
        % 3.7 线性求解
        Model.Solver.Linear.Symmetric = true; 
        dx = linear_solve(K_step, RHS_step, Model);
        
        % 3.8 更新解
        x = x + damping * dx;
    end
    
    if iter == MaxIter
        warning('非线性迭代达到最大步数，可能未完全收敛。');
    end
    
    info.Iterations = iter;
    info.Residuals = conv_history;
end