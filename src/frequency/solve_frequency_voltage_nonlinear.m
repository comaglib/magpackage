function [Solution, Info] = solve_frequency_voltage_nonlinear(Model, Coil, CircuitProps, FreqParams)
% SOLVE_FREQUENCY_VOLTAGE_NONLINEAR 非线性时谐场路耦合求解器
%
% 方法: Effective Reluctivity + Voltage Driven Coupling
% 求解: 耦合系统 [K(nu) -F(nu); W' Z] * [A; I] = [0; U]
%
% 输入:
%   CircuitProps: .R, .L_leak, .U (电压相量)
%   FreqParams:   .Frequency, .MaxIter, .Tol

    fprintf('==============================================\n');
    fprintf('   Nonlinear Freq. Solver (Voltage Driven)    \n');
    fprintf('==============================================\n');

    freq = FreqParams.Frequency;
    omega = 2 * pi * freq;
    U_phasor = CircuitProps.U;
    
    max_iter = 20; if isfield(FreqParams, 'MaxIter'), max_iter = FreqParams.MaxIter; end
    tol = 1e-4; if isfield(FreqParams, 'Tol'), tol = FreqParams.Tol; end
    
    fprintf('  - Freq: %.1f Hz, Voltage: %.2e V\n', freq, abs(U_phasor));

    % 1. 准备不变量
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    % 恒定矩阵
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    M_reg = assemble_mass_matrix(Model, 'All');
    
    % 绕组向量 W (几何相关，恒定)
    fprintf('  [Init] Assembling Winding Vector...\n');
    W_vec_A = assemble_winding_vector(Model, Coil);
    W_vec = [W_vec_A; sparse(numActiveV, 1)];
    
    % 感应源单位向量 b_ind_unit (恒定)
    CoilUnit = Coil; CoilUnit.I(:) = 1.0;
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    b_ind_unit = -1i * omega * (M_sigma * As_unit_vec);
    
    % 初始磁阻率 (Linear)
    Nu_vec = []; 
    
    % 迭代变量
    % 系统解向量 X = [A; V; I]
    
    % 自适应阻尼参数
    relax = 0.5;
    prev_diff = inf;
    
    % 2. 迭代主循环
    for iter = 1:max_iter
        
        % --- A. 组装依赖 Nu 的矩阵 ---
        % K(nu)
        K = assemble_magnetic_stiffness(Model, Nu_vec);
        
        % F(nu) = b_mag_unit(nu) + b_ind_unit
        % 注意: assemble_rhs_reduced 计算的是 b_mag
        b_mag_unit = assemble_rhs_reduced(Model, CoilUnit, Nu_vec);
        
        F_vec_A = b_mag_unit + b_ind_unit;
        F_vec = [F_vec_A; sparse(numActiveV, 1)];
        
        % --- B. 构建耦合系统 ---
        scale_K = mean(abs(nonzeros(K)));
        eps_A = 1e-6 * scale_K;
        
        K_eff = K + 1i * omega * M_sigma + eps_A * M_reg;
        C_eff = C_sigma; 
        Ct_eff = 1i * omega * C_sigma'; 
        G_eff = 1i * omega * G_sigma;
        
        % 场块: [K C; C' G]
        Block_Field = [K_eff, C_eff; Ct_eff, G_eff];
        numFieldDoFs = size(Block_Field, 1);
        
        % 电路块
        Z_circ = CircuitProps.R + 1i * omega * CircuitProps.L_leak;
        Coupling_Row = 1i * omega * W_vec';
        
        % 全局矩阵
        % [ Field   -F ]
        % [ W'*jw    Z ]
        SysK = [Block_Field, -F_vec; 
                Coupling_Row, sparse(Z_circ)];
            
        SysRHS = [sparse(numFieldDoFs, 1); U_phasor];
        
        % --- C. 边界条件与求解 ---
        if isfield(Model.Runtime,'FixedEdges'), fixed_edges=Model.Runtime.FixedEdges; else, fixed_edges=[]; end
        fixed_vals = zeros(length(fixed_edges), 1);
        
        [SysK_bc, SysRHS_bc] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, fixed_vals);
        
        Model.Solver.Linear.Interface='MUMPS'; 
        Model.Solver.Linear.Symmetric=false; % 非对称耦合
        if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
        Model.Solver.Linear.MumpsICNTL.i8=77;
        
        x_sol = linear_solve(SysK_bc, SysRHS_bc, Model);
        
        % --- D. 提取结果 ---
        A_new = x_sol(1:numEdges);
        I_new = x_sol(end);
        
        % --- E. 更新材料属性 ---
        % 需要用新的电流 I_new 来计算实际的 Bs
        % 更新 Coil 对象的电流
        CoilCurrent = Coil;
        CoilCurrent.I(:) = I_new; % 复数电流
        
        [Nu_new, max_B] = update_reluctivity_parallel(Model, A_new, CoilCurrent);
        
        % --- F. 检查收敛 ---
        if isempty(Nu_vec)
            Nu_vec = Nu_new; 
            diff = 1.0;
            fprintf('    Step 0 (Linear): I=%.2f A, Max|B|=%.4f T\n', abs(I_new), max_B);
        else
            diff = norm(Nu_new - Nu_vec) / norm(Nu_vec);
            
            % 自适应阻尼
            if diff < prev_diff
                relax = min(relax * 1.1, 1.0);
                status = 'Acc';
            else
                relax = max(relax * 0.33, 0.05);
                status = 'Damp';
            end
            
            fprintf('    Iter %2d: |I|=%.2f A, Max|B|=%.4f T, d(nu)=%.4e [%s]\n', ...
                iter, abs(I_new), max_B, diff, status);
            
            if diff < tol
                fprintf('  [Converged]\n');
                x_curr = x_sol;
                break;
            end
            
            Nu_vec = relax * Nu_new + (1-relax) * Nu_vec;
            prev_diff = diff;
        end
        
        x_curr = x_sol;
    end
    
    % 输出结果
    Solution.A = x_curr(1:numEdges);
    Solution.V_active = x_curr(numEdges+1 : numEdges+numActiveV);
    Solution.I = x_curr(end);
    
    Info.Impedance = U_phasor / Solution.I;
    Info.FinalNu = Nu_vec;
    Info.Iterations = iter;
    
    fprintf('  [Result] I = %.4f + j%.4f A\n', real(Solution.I), imag(Solution.I));
    fprintf('  [Result] Z = %.4e + j%.4e Ohm\n', real(Info.Impedance), imag(Info.Impedance));
end