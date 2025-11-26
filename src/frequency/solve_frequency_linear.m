function [Solution, Info] = solve_frequency_linear(Model, Coil, FreqParams)
% SOLVE_FREQUENCY_LINEAR 线性时谐磁场求解器 (AC Magnetic, 修正版)
% 
% 求解方程:
% [ K + jw*M_sig + eps*M_reg    C_sig ] [ A ] = [ b_mag ]
% [ jw * C_sig'                 jw*G  ] [ V ]   [ 0     ]
%
% 考虑了场路耦合的基础:
% 1. 支持复数电流输入 (幅值+相位)
% 2. 输出包含总复功率，便于后续计算阻抗 Z
% 
% 修复: 增加了感应源项 (-jw * M * As)，解决非磁性导体 RHS 为 0 的问题。
% 修复: 增大了空气域正则化系数，防止低频下 A 场发散。


    fprintf('==============================================\n');
    fprintf('   Linear Frequency Domain Solver (AC Magnetic)\n');
    fprintf('==============================================\n');

    freq = FreqParams.Frequency;
    omega = 2 * pi * freq;
    fprintf('  - Frequency: %.2f Hz\n', freq);

    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    % 1. 矩阵组装
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    K = assemble_magnetic_stiffness(Model);
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    
    % [修正] 增强正则化: 1e-6
    M_reg = assemble_mass_matrix(Model, 'All');
    scale_K = mean(abs(nonzeros(K)));
    eps_A = 1e-6 * scale_K; 
    
    % 2. 系统矩阵
    fprintf('  [Freq] Building System Matrix...\n');
    
    K_eff = K + 1i * omega * M_sigma + eps_A * M_reg;
    C_eff = C_sigma;
    Ct_eff = 1i * omega * C_sigma'; 
    G_eff  = 1i * omega * G_sigma;
    
    SysK = [K_eff, C_eff; Ct_eff, G_eff];
    
    % 3. 组装右端项
    fprintf('  [Freq] Assembling Source...\n');
    
    b_mag = assemble_rhs_reduced(Model, Coil);
    
    % 感应源
    CoilUnit = Coil; CoilUnit.I(:) = 1.0;
    As_unit = project_source_A_on_edges(Model, CoilUnit);
    I_phasor = Coil.I(1); 
    As_vec = As_unit * I_phasor;
    
    b_ind = -1i * omega * (M_sigma * As_vec);
    b_total = b_mag + b_ind;
    
    SysRHS = [b_total; sparse(numActiveV, 1)];
    
    % 4. 边界条件
    if isfield(Model, 'Runtime') && isfield(Model.Runtime, 'FixedEdges')
        fixed_edges = Model.Runtime.FixedEdges;
    else
        fixed_edges = [];
    end
    fixed_vals = zeros(length(fixed_edges), 1);
    
    [SysK_bc, SysRHS_bc] = apply_dirichlet_bc(SysK, SysRHS, fixed_edges, fixed_vals);
    
    % 5. 求解
    fprintf('  [Freq] Solving...\n');
    Model.Solver.Linear.Interface = 'MUMPS'; 
    Model.Solver.Linear.Symmetric = true; 
    if ~isfield(Model.Solver.Linear, 'MumpsICNTL')
        Model.Solver.Linear.MumpsICNTL = struct();
    end
    Model.Solver.Linear.MumpsICNTL.i8 = 77;
    
    x_sol = linear_solve(SysK_bc, SysRHS_bc, Model);
    
    Solution.A = x_sol(1:numEdges);
    Solution.V_active = x_sol(numEdges+1:end); 
    
    Info.SystemSize = length(x_sol);
    fprintf('  [Freq] Done. |A|_max = %.2e\n', max(abs(Solution.A)));
end