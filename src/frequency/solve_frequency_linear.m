function [Solution, Info] = solve_frequency_linear(Model, Coil, FreqParams)
% SOLVE_FREQUENCY_LINEAR [Final Fix: Direct Source]
% 
% 修正:
% 1. 移除电流缩放逻辑，直接使用 Coil 中的真实电流计算 RHS。
% 2. 保持空气电导率正则化逻辑。

    fprintf('==============================================\n');
    fprintf('   Linear Frequency Solver (Direct Source)    \n');
    fprintf('==============================================\n');

    freq = FreqParams.Frequency;
    omega = 2 * pi * freq;

    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    numActiveV = DofData.NumActiveNodes;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 1. 组装矩阵
    fprintf('  [Matrix] Assembling Stiffness (K)...\n');
    K = assemble_magnetic_stiffness(Model);
    
    % 使用 'Physical' 模式，包含空气微小电导率以防奇异
    fprintf('  [Matrix] Assembling Mass (M_sigma)...\n');
    M_sigma = assemble_mass_matrix(Model, 'Physical'); 
    
    fprintf('  [Matrix] Assembling Laplacian (G) & Coupling (C)...\n');
    G_sigma = assemble_scalar_laplacian(Model, DofData);
    C_sigma = assemble_coupling_matrix(Model, 'Physical', DofData);
    
    % A 场正则化
    M_reg = assemble_mass_matrix(Model, 'All');
    scale_K = mean(abs(nonzeros(K)));
    eps_A = 1e-8 * scale_K; 
    
    % 2. 系统构建
    K_eff = K + 1i * omega * M_sigma + eps_A * M_reg;
    SysK = [K_eff, C_sigma; 1i * omega * C_sigma', 1i * omega * G_sigma];
    
    % 3. 右端项
    fprintf('  [RHS] Assembling Source Terms...\n');
    
    % [CRITICAL FIX] 直接投影带真实电流的 Coil
    % 不再进行 CoilUnit * I_phasor 的操作，防止电流丢失
    As_edges = project_source_A_on_edges(Model, Coil);
    
    % 检查源项强度
    max_As = max(abs(As_edges));
    fprintf('         源场投影强度检查: |As|_max = %.4e Weber\n', max_As);
    if max_As < 1e-9
        warning('检测到源场投影几乎为零！请检查线圈电流或几何位置。');
    end
    
    % 感应源项: -jw * sigma * A_source
    b_ind = -1i * omega * (M_sigma * As_edges);
    
    % 磁化源项: TEAM 7 为非磁性，设为 0
    b_mag = sparse(numEdges, 1);
    
    SysRHS = [b_mag + b_ind; sparse(numActiveV, 1)];
    
    % 4. 边界条件
    if isfield(Model.Runtime, 'FixedEdges'), fe = Model.Runtime.FixedEdges; else, fe = []; end
    [SysK, SysRHS] = apply_dirichlet_bc(SysK, SysRHS, fe, zeros(size(fe)));
    
    % 5. 求解
    fprintf('  [Solver] Solving Linear System (Size: %d)...\n', size(SysK, 1));
    Model.Solver.Linear.Interface = 'MUMPS'; 
    Model.Solver.Linear.Symmetric = true; 
    
    x = linear_solve(SysK, SysRHS, Model);
    
    Solution.A = x(1:numEdges);
    Solution.V_active = x(numEdges+1:end); 
    
    Info.SystemSize = size(SysK, 1);
    Info.DofData = DofData;
    
    fprintf('  [Freq] Done. |A|_max = %.2e\n', max(abs(Solution.A)));
end