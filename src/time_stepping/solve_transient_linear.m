function [Solution, Info] = solve_transient_linear(Model, Coil, TimeParams)
% SOLVE_TRANSIENT_LINEAR (压缩 DoF 优化版)
% 
% 输入:
%   Model      - 模型结构体
%   Coil       - 线圈对象 (包含几何信息，电流 I 将被时间函数覆盖)
%   TimeParams - 时间控制结构体
%       .dt        - 时间步长 (s)
%       .N_steps   - 总步数
%       .I_func    - 电流时间函数句柄 @(t) -> current_val
%
% 输出:
%   Solution   - 包含每一步的 A 和 V (可选: 只存最后一步或关键点)
%   Info       - 诊断信息
% 
% 修正:
% 1. 增加了感应电动势源项 (-sigma * dA/dt)，解决了非磁性导体涡流为0的问题
% 2. 修复了逻辑冗余和不可达代码
% 3. 增加了绝缘区 V 场的正则化 (V-Regularization)，防止矩阵奇异导致的数值爆炸。
% 4. 清理了 SymmetricForm 的冗余分支和死代码。
% 5. 修正了 MumpsICNTL 的冗余设置。
% 6. 仅在导电区域求解 V，彻底消除空气区奇异性

    fprintf('==============================================\n');
    fprintf('   Transient Eddy Current Solver (Linear A-V) \n');
    fprintf('==============================================\n');

    dt = TimeParams.dt;
    N_steps = TimeParams.N_steps;
    I_func = TimeParams.I_func;
    
    % 1. 生成自由度映射
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    % 2. 组装物理矩阵 (使用 DofData)
    fprintf('  [Transient] Assembling physical matrices...\n');
    K = assemble_magnetic_stiffness(Model);
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    
    % A 场正则化 (空气中 Curl-Curl 零空间仍需处理)
    M_reg_A = assemble_mass_matrix(Model, 'All');
    scale_K = mean(abs(nonzeros(K)));
    eps_A = 1e-8 * scale_K;
    
    % V 场正则化? 不需要了！因为 G_sigma 现在是正定的 (仅含导体)
    
    % 3. 构建系统矩阵
    fprintf('  [Transient] Building System Matrix (Symmetric)...\n');
    
    K_eff = K + (1/dt) * M_sigma + eps_A * M_reg_A;
    C_eff = C_sigma;
    Ct_eff = C_sigma';     
    G_eff  = G_sigma * dt; 
    
    % 维度检查:
    % K_eff: Ne x Ne
    % C_eff: Ne x Nv
    % Ct_eff: Nv x Ne
    % G_eff: Nv x Nv
    
    SysK = [K_eff, C_eff; Ct_eff, G_eff];
    
    % 4. 预计算源场投影
    CoilUnit = Coil;
    CoilUnit.I(:) = 1.0;
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    
    % 5. 边界条件
    % 仅需处理 A 的 BC。V 在导电区内部自然连续，边界自然 Neumann
    if isfield(Model, 'Runtime') && isfield(Model.Runtime, 'FixedEdges')
        fixed_edges = Model.Runtime.FixedEdges;
    else
        fixed_edges = [];
    end
    fixed_vals_A = zeros(length(fixed_edges), 1);
    
    % 6. 求解器配置
    Model.Solver.Linear.Interface = 'MUMPS'; 
    Model.Solver.Linear.Symmetric = true;
    if ~isfield(Model.Solver.Linear, 'MumpsICNTL')
        Model.Solver.Linear.MumpsICNTL = struct();
    end
    Model.Solver.Linear.MumpsICNTL.i8 = 77;
    
    % 7. 时间循环
    fprintf('  [Transient] Starting Time Stepping (dt=%.1e, Steps=%d)...\n', dt, N_steps);
    
    A_prev = zeros(numEdges, 1);
    I_prev = 0;
    if nargin > 2 && isfield(TimeParams, 'I_init'), I_prev = TimeParams.I_init; end
    
    for n = 1:N_steps
        t_now = n * dt;
        I_now = I_func(t_now);
        
        % RHS Part A
        Coil.I(:) = I_now;
        b_mag = assemble_rhs_reduced(Model, Coil);
        
        dI_dt = (I_now - I_prev) / dt;
        dAs_dt_vec = As_unit_vec * dI_dt;
        b_ind = - M_sigma * dAs_dt_vec;
        b_ext = b_mag + b_ind;
        
        % RHS Part B
        rhs_A = b_ext + (1/dt) * M_sigma * A_prev;
        rhs_V = C_sigma' * A_prev; % V 维度为 Nv
        
        SysRHS = [rhs_A; rhs_V];
        
        [SysK_bc, SysRHS_bc] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, fixed_vals_A);
        
        x_curr = linear_solve(SysK_bc, SysRHS_bc, Model);
        
        A_curr = x_curr(1:numEdges);
        % V_curr = x_curr(numEdges+1:end); % 仅含 Active V
        
        A_prev = A_curr;
        I_prev = I_now;
        
        if mod(n, 1) == 0
            fprintf('    Step %d/%d (t=%.2e, I=%.1f) |A|_max=%.2e\n', ...
                n, N_steps, t_now, I_now, max(abs(A_curr)));
        end
        
        Solution.A = A_curr;
        % Solution.V ... 需要映射回全局节点才方便绘图，暂时存 raw
    end
    
    Info.FinalTime = t_now;
end