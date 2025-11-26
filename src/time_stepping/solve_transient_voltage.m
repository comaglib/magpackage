function [Solution, Info] = solve_transient_voltage(Model, Coil, CircuitProps, TimeParams)
% SOLVE_TRANSIENT_VOLTAGE 电压驱动瞬态求解器 (强耦合)
%
% 求解: A_r 和 I
%
% 输入:
%   CircuitProps: .R (电阻), .L_leak (漏感/外感), .V_func (电压函数 handle)

    fprintf('==============================================\n');
    fprintf('   Transient Solver (Voltage Driven)          \n');
    fprintf('==============================================\n');

    dt = TimeParams.dt;
    N_steps = TimeParams.N_steps;
    V_func = CircuitProps.V_func;
    R_circ = CircuitProps.R;
    L_circ = CircuitProps.L_leak;
    
    numEdges = size(Model.Mesh.Edges, 2);
    
    % 1. 组装基础矩阵
    K = assemble_magnetic_stiffness(Model);
    M = assemble_mass_matrix(Model, 'ConductorOnly'); % 涡流区
    % 稳定化
    M_reg = assemble_mass_matrix(Model, 'All');
    eps_A = 1e-8 * mean(abs(nonzeros(K)));
    
    % A 方程左端 (不含 I)
    K_sys_A = K + (1/dt)*M + eps_A*M_reg;
    
    % 2. 组装耦合向量
    % Q: 电流对 A 的源项贡献 (Unit Current RHS)
    % 只有磁性材料产生 Q (nu0-nu)。如果是空心线圈，Q=0。
    CoilUnit = Coil; CoilUnit.I(:)=1.0;
    Q_vec = assemble_rhs_reduced(Model, CoilUnit);
    
    % W: A 对电路的磁链贡献 (Winding Vector)
    W_vec = assemble_winding_vector(Model, Coil);
    
    % 3. 构建耦合系统矩阵
    % Row 1 (Field):   K_sys_A * A - Q * I = M/dt * A_prev
    % Row 2 (Circuit): W' * (A - A_prev)/dt + R*I + L*(I - I_prev)/dt = V
    %
    % 整理 LHS:
    % [ K_sys_A    -Q           ] [ A_n+1 ]
    % [ W'/dt      R + L/dt     ] [ I_n+1 ]
    
    % 注意: W'/dt 这一项表示 d(Psi_r)/dt。
    % 这里的符号定义：Circuit Eq: V = RI + L dI/dt + d(Psi)/dt
    % 所以 dPsi/dt 项移到左边是 + W'/dt
    
    Block_AA = K_sys_A;
    Block_AI = -Q_vec;
    Block_IA = (1/dt) * W_vec';
    Block_II = sparse([R_circ + L_circ/dt]);
    
    SysK = [Block_AA, Block_AI; Block_IA, Block_II];
    
    % 4. 边界条件 (仅针对 A)
    if isfield(Model.Runtime, 'FixedEdges')
        fixed_edges = Model.Runtime.FixedEdges;
    else
        fixed_edges = [];
    end
    fixed_vals = zeros(length(fixed_edges), 1);
    
    % 扩展 BC 到 I (I 是自由的，不受 Dirichlet 约束)
    % 但 apply_dirichlet_bc 需要知道总维度。
    % 我们只约束前 numEdges 个自由度。
    
    % 5. 时间步进
    fprintf('  [Transient] Starting Steps (dt=%.1e)...\n', dt);
    
    A_prev = zeros(numEdges, 1);
    I_prev = 0;
    
    TimeHist.I = zeros(N_steps, 1);
    
    for n = 1:N_steps
        t_now = n * dt;
        V_now = V_func(t_now);
        
        % 构建 RHS
        % RHS_A = (M/dt) * A_prev
        rhs_A = (1/dt) * M * A_prev;
        
        % RHS_I = V_now + (W'/dt)*A_prev + (L/dt)*I_prev
        % (注意移项: dPsi/dt = (Psi_new - Psi_old)/dt. Old 项移到右边是 +)
        rhs_I = V_now + (1/dt)*(W_vec' * A_prev) + (L_circ/dt)*I_prev;
        
        SysRHS = [rhs_A; rhs_I];
        
        % 施加 BC
        [K_step, R_step] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, fixed_vals);
        
        % 求解
        % I 是最后一个未知数
        % 矩阵可能不对称 (右上 -Q, 左下 W'/dt)。除非 Q ~ -W (Reciprocity?)
        % 在 A-formulation 中，Q 和 W 通常不相等。使用 Unsymmetric solver。
        Model.Solver.Linear.Interface = 'MUMPS';
        Model.Solver.Linear.Symmetric = false; 
        Model.Solver.Linear.MumpsICNTL.i8 = 77;
        
        x_sol = linear_solve(K_step, R_step, Model);
        
        A_curr = x_sol(1:end-1);
        I_curr = x_sol(end);
        
        % 更新
        A_prev = A_curr;
        I_prev = I_curr;
        
        TimeHist.I(n) = I_curr;
        
        if mod(n, 10) == 0 || n==1
            fprintf('    Step %d: V=%.1f V, I=%.4f A\n', n, V_now, I_curr);
        end
        
        Solution.A = A_curr;
    end
    
    Solution.I_history = TimeHist.I;
    Info.FinalTime = t_now;
end