function [Solution, Info] = solve_frequency_voltage(Model, Coil, CircuitProps, FreqParams)
% SOLVE_FREQUENCY_VOLTAGE 频域电压驱动求解器 (场路耦合) (日志优化版)

    fprintf('==============================================\n');
    fprintf('   Freq. Domain Solver (Voltage Driven)       \n');
    fprintf('==============================================\n');

    freq = FreqParams.Frequency;
    omega = 2 * pi * freq;
    U_phasor = CircuitProps.U;
    
    % [修改] 使用科学计数法显示电压
    fprintf('  - Freq: %.1f Hz, Voltage: %.2e V\n', freq, abs(U_phasor));

    % 1. 基础矩阵组装
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    K = assemble_magnetic_stiffness(Model);
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    
    M_reg = assemble_mass_matrix(Model, 'All');
    eps_A = 1e-6 * mean(abs(nonzeros(K))); 
    
    % 2. 构建 A-V 系统块 (K_sys)
    K_eff = K + 1i * omega * M_sigma + eps_A * M_reg;
    C_eff = C_sigma; 
    Ct_eff = 1i * omega * C_sigma'; 
    G_eff = 1i * omega * G_sigma;
    
    Block_Field = [K_eff, C_eff; Ct_eff, G_eff];
    numFieldDoFs = size(Block_Field, 1);
    
    % 3. 组装耦合向量
    fprintf('  [Freq] Assembling Coupling Vectors...\n');
    
    % 3.1 绕组向量 W
    W_vec_A = assemble_winding_vector(Model, Coil);
    W_vec = [W_vec_A; sparse(numActiveV, 1)];
    
    % 3.2 源耦合向量 F
    CoilUnit = Coil; CoilUnit.I(:) = 1.0;
    
    b_mag_unit = assemble_rhs_reduced(Model, CoilUnit);
    
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    b_ind_unit = -1i * omega * (M_sigma * As_unit_vec);
    
    F_vec_A = b_mag_unit + b_ind_unit;
    F_vec = [F_vec_A; sparse(numActiveV, 1)];
    
    % 4. 构建全局场路耦合系统
    Z_circ = CircuitProps.R + 1i * omega * CircuitProps.L_leak;
    Coupling_Row = 1i * omega * W_vec';
    
    SysK = [Block_Field, -F_vec; 
            Coupling_Row, sparse(Z_circ)];
        
    SysRHS = [sparse(numFieldDoFs, 1); U_phasor];
    
    % 5. 边界条件 (仅对 A)
    if isfield(Model.Runtime, 'FixedEdges')
        fixed_edges = Model.Runtime.FixedEdges;
    else
        fixed_edges = [];
    end
    fixed_vals = zeros(length(fixed_edges), 1);
    
    [SysK_bc, SysRHS_bc] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, fixed_vals);
    
    % 6. 求解
    fprintf('  [Freq] Solving coupled system...\n');
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = false; 
    if ~isfield(Model.Solver.Linear, 'MumpsICNTL'), Model.Solver.Linear.MumpsICNTL=struct(); end
    Model.Solver.Linear.MumpsICNTL.i8 = 77;
    
    x_sol = linear_solve(SysK_bc, SysRHS_bc, Model);
    
    % 7. 结果拆分
    Solution.A = x_sol(1:numEdges);
    Solution.V_active = x_sol(numEdges+1 : numEdges+numActiveV);
    Solution.I = x_sol(end);
    
    Info.Impedance = U_phasor / Solution.I;
    
    fprintf('  [Result] Current I = %.4f + j%.4f A (|I|=%.4f)\n', ...
        real(Solution.I), imag(Solution.I), abs(Solution.I));
    fprintf('  [Result] Equiv Z   = %.4e + j%.4e Ohm\n', ...
        real(Info.Impedance), imag(Info.Impedance));
end