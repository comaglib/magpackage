function solve_nonlinear_ac()
% SOLVE_NONLINEAR_AC 测试非线性时谐求解器
%
% 案例: 铁心线圈在 AC 激励下的饱和特性
% 验证: 随着电流增大，铁心进入饱和，等效电感 L 应下降。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Nonlinear AC Magnetic Test (Saturation)    \n');
    fprintf('==============================================\n');
    
    % 1. 网格 (铁心+空气)
    msh_file = 'test_nonlinear_ac.msh';
    create_iron_air_mesh(msh_file); % 复用之前的铁心网格
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 2. 材料 (定义 BH 曲线)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Steel';
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Sigma = 0; % 忽略涡流，专注磁饱和
    % 定义 BH 曲线
    H_data = [0, 100, 200, 500, 1000, 5000, 10000, 50000];
    B_data = [0, 0.8, 1.2, 1.5, 1.65, 1.85, 1.95, 2.2];
    Model.Materials.Lib(2).BH_Curve.H = H_data;
    Model.Materials.Lib(2).BH_Curve.B = B_data;
    
    Model = preprocess_materials(Model); % 预处理
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.5; CoilParams.Radius = 0.2; 
    CoilParams.N_seg = 36;
    % Current 将在循环中设定
    
    % 4. 边界
    Model.Runtime.FixedEdges = identify_boundary_edges(Model.Mesh, Raw, 100);
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 5. 电流扫描 (Linear vs Saturated)
    currents = [100, 50000]; % 小电流(线性区), 大电流(饱和区)
    
    fprintf('\n--- Saturation Test (Frequency = 50 Hz) ---\n');
    
    L_results = [];
    
    for I_mag = currents
        fprintf('\n>> Testing Current I = %.0f A\n', I_mag);
        
        CoilParams.Current = I_mag;
        Coil = create_racetrack_coil(CoilParams);
        
        FreqParams.Frequency = 50;
        FreqParams.MaxIter = 25;
        FreqParams.Tol = 1e-3;
        
        % 调用非线性求解器
        [Sol, Info] = solve_frequency_nonlinear(Model, Coil, FreqParams);
        
        % 提取电感
        % 注意: extract_circuit_parameters 是基于能量法的
        CircuitParams = extract_circuit_parameters(Model, Sol, Coil, 50);
        
        L_results(end+1) = CircuitParams.L;
    end
    
    % 6. 验证
    fprintf('\n--- Result Comparison ---\n');
    fprintf('L_linear (100A)   = %.4e H\n', L_results(1));
    fprintf('L_sat    (50000A) = %.4e H\n', L_results(2));
    
    ratio = L_results(2) / L_results(1);
    fprintf('Inductance Ratio (Sat/Lin) = %.2f\n', ratio);
    
    if ratio < 0.5
        fprintf('[PASS] Inductance dropped significantly due to saturation.\n');
    else
        fprintf('[WARN] Saturation effect weak. Check B-H curve or current level.\n');
    end
end

function create_iron_air_mesh(filename)
    % 简单的铁心网格 (ID X Y Z 格式)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n');
    fprintf(fid, '1 0 0 0\n2 0.1 0 0\n3 0 0.1 0\n4 0 0 0.1\n5 0.2 0.2 0.2\n');
    fprintf(fid, '$EndNodes\n');
    fprintf(fid, '$Elements\n3 3 1 3\n');
    fprintf(fid, '2 100 2 1\n1 2 3 5\n');
    fprintf(fid, '3 2 4 1\n2 1 2 3 4\n'); % Iron
    fprintf(fid, '3 1 4 1\n3 2 3 4 5\n'); % Air
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end