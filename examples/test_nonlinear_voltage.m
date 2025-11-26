function test_nonlinear_voltage()
% TEST_NONLINEAR_VOLTAGE 测试非线性场路耦合 (修正版)
% 
% 场景: 给定不同电压幅值，观察电流和等效电感的变化。
% 预期: 电压增大 -> 铁心饱和 -> 电感 L 下降 -> 电流 I 增加幅度超过线性比例。
%
% 修复: 
% 1. 内置 Robust Delaunay Mesh Generator，解决 Tet=0 问题。
% 2. 确保网格包含空气和铁心，且边界 Tag 正确。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Nonlinear Field-Circuit Coupling Test      \n');
    fprintf('==============================================\n');
    
    % 1. 生成高质量测试网格
    msh_file = 'test_nl_volt_robust.msh';
    create_delaunay_mesh_robust(msh_file); 
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    
    % 完整性检查
    if isempty(Raw.T)
        error('Fatal: Mesh file contains no tetrahedra. Mesh generation failed.');
    end
    
    Model.Mesh = build_topology(Raw);
    
    % 2. 材料
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0; 
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Steel';
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Sigma = 0; % 忽略涡流，专注磁饱和
    Model.Materials.Lib(2).Mu_r = 6000; % 初始线性值
    
    % B-H 曲线
    H_data = [0, 100, 200, 500, 1000, 5000, 10000, 50000];
    B_data = [0, 0.8, 1.2, 1.5, 1.65, 1.85, 1.95, 2.2];
    Model.Materials.Lib(2).BH_Curve.H = H_data;
    Model.Materials.Lib(2).BH_Curve.B = B_data;
    
    Model = preprocess_materials(Model);
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 线圈
    CoilParams.Center = [0, 0, 0]; 
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.4; 
    CoilParams.Radius = 0.3; 
    CoilParams.N_seg = 72;
    CoilParams.Current = 0; 
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 边界
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    if isempty(fixed_edges)
        warning('Boundary detection failed. Using natural BC.');
    else
        fprintf('  - Boundary edges identified: %d\n', length(fixed_edges));
    end
    Model.Runtime.FixedEdges = fixed_edges;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 5. 电压扫描
    voltages = [100, 10000]; % 100V (Linear), 10000V (Saturated)
    freq = 50;
    
    FreqParams.Frequency = freq;
    FreqParams.MaxIter = 25;
    FreqParams.Tol = 1e-3;
    
    CircuitProps.R = 1e-3;
    CircuitProps.L_leak = 0;
    
    Results = [];
    
    fprintf('\n--- Voltage Sweep (50 Hz) ---\n');
    
    for U_mag = voltages
        fprintf('\n>> Testing Voltage U = %.0f V\n', U_mag);
        CircuitProps.U = U_mag;
        
        % 调用非线性电压驱动求解器
        [Sol, Info] = solve_frequency_voltage_nonlinear(Model, Coil, CircuitProps, FreqParams);
        
        % 提取结果
        I_mag = abs(Sol.I);
        Z = Info.Impedance;
        L_eff = imag(Z) / (2*pi*freq);
        
        Res.U = U_mag;
        Res.I = I_mag;
        Res.L = L_eff;
        Results = [Results, Res];
    end
    
    % 6. 验证
    fprintf('\n--- Comparison ---\n');
    fprintf('Low Voltage (%.0fV)  : I = %.2f A, L = %.4e H\n', Results(1).U, Results(1).I, Results(1).L);
    fprintf('High Voltage (%.0fV) : I = %.2f A, L = %.4e H\n', Results(2).U, Results(2).I, Results(2).L);
    
    % 检查 1: 非线性增加
    I_ratio = Results(2).I / Results(1).I;
    U_ratio = Results(2).U / Results(1).U; % 100
    fprintf('Current Ratio: %.2f (Voltage Ratio: %.2f)\n', I_ratio, U_ratio);
    
    % 检查 2: 电感下降
    L_ratio = Results(2).L / Results(1).L;
    fprintf('Inductance Ratio: %.2f\n', L_ratio);
    
    if I_ratio > U_ratio * 1.2 || L_ratio < 0.8
        fprintf('[PASS] Saturation detected (Current spiked, Inductance dropped).\n');
    else
        fprintf('[WARN] Saturation not significant. Check BH curve or Voltage levels.\n');
    end
end
