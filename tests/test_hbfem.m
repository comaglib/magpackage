function test_hbfem()
% TEST_HBFEM 测试谐波平衡法求解器 (验证版)
%
% 功能:
% 1. 验证高次谐波产生 (非线性材料)。
% 2. 验证 TimeSteps 自适应逻辑。
% 3. 验证能量守恒 (Parseval's Theorem)。
%
% 修复: 解决了材料库结构体字段不一致导致的赋值报错。

    % 确保路径
    if ~exist('solve_hbfem', 'file')
        addpath(genpath(pwd)); 
    end

    fprintf('==============================================\n');
    fprintf('   HBFEM Verification: Nonlinear + Parseval   \n');
    fprintf('==============================================\n');
    
    % 1. 网格与材料
    msh_file = 'test_hbfem.msh';
    if exist('create_test_mesh_direct', 'file')
        create_test_mesh_direct(msh_file); 
    end
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % --- [修复] 统一初始化材料库结构体 ---
    % 必须确保 Lib(1) 和 Lib(2) 拥有完全相同的字段
    
    % Material 1: Air (Linear)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1;
    Model.Materials.Lib(1).Sigma = 0;
    Model.Materials.Lib(1).BH_Curve = []; % 占位，保持结构一致
    
    % Material 2: Steel (Nonlinear)
    Model.Materials.Lib(2).Name = 'Steel';
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Mu_r = 1000;   % 初始线性磁导率估计
    Model.Materials.Lib(2).Sigma = 0;
    
    % 定义 B-H 曲线
    Model.Materials.Lib(2).BH_Curve.H = [0, 100, 200, 500, 1000, 5000, 10000]; 
    Model.Materials.Lib(2).BH_Curve.B = [0, 0.8, 1.2, 1.5, 1.65, 1.85, 1.95];
    
    Model = preprocess_materials(Model);
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 2. 线圈 (高饱和电流)
    CoilParams.Center = [0,0,0]; 
    CoilParams.Normal = [0,0,1];
    CoilParams.Length = 0.4; 
    CoilParams.Radius = 0.3; 
    CoilParams.N_seg = 72; 
    CoilParams.Turns = 50;
    
    Current_Amp = 4000; % 4000A 强制饱和
    CoilParams.Current = Current_Amp; 
    Coil = create_racetrack_coil(CoilParams);
    Coil.HarmonicCurrents = [Current_Amp; 0; 0]; 
    
    % 3. 参数设置 (测试自适应)
    HBFEMParams.Frequency = 50;
    HBFEMParams.Harmonics = [1, 3, 5]; 
    HBFEMParams.TimeSteps = []; % 留空以触发自适应调整 (Auto-Setup)
    
    if isfield(Raw, 'FaceTags')
        fe = identify_boundary_edges(Model.Mesh, Raw, 100); 
    else
        fe = []; 
    end
    Model.Runtime.FixedEdges = fe;
    
    % 4. 求解
    [Sol, Info] = solve_hbfem(Model, Coil, HBFEMParams);
    
    % 5. 结果分析
    mag_1 = norm(Sol.Harmonics(:, 1));
    mag_3 = norm(Sol.Harmonics(:, 2));
    fprintf('\nHarmonics (Norm): 1st=%.2e, 3rd=%.2e (Distortion=%.2f%%)\n', ...
        mag_1, mag_3, (mag_3/mag_1)*100);
    
    % --- 6. 能量守恒验证 (Parseval) ---
    fprintf('\n----------------------------------------------\n');
    fprintf('   Energy Conservation Check (Parseval)       \n');
    fprintf('----------------------------------------------\n');
    
    % 随机选取一个非零解的 Edge 进行验证
    vals = vecnorm(Sol.Harmonics, 2, 2);
    [~, edge_idx] = max(vals); % 选能量最大的边
    
    if isempty(edge_idx) || vals(edge_idx) == 0
        fprintf('[WARN] Solution is zero. Cannot check Parseval.\n');
    else
        coeffs_A = Sol.Harmonics(edge_idx, :); % [1, 3, 5]
        
        % 重建时域
        AFT = aft_module();
        % 注意: 需要使用求解器实际计算出的 TimeSteps
        Info_Check = AFT.prepare_indices(HBFEMParams.Harmonics, Info.TimeSteps);
        val_t = AFT.freq2time(coeffs_A, Info_Check);
        
        % 计算能量
        % Time Domain Energy: sum(x^2)/N (均值)
        E_time = mean(val_t.^2);
        
        % Freq Domain Energy: sum(|c_k|^2) (Parseval)
        % 对于实信号单边谱(不含DC): Energy = 2 * sum(|c_k|^2)
        % 如果有 DC，则为 |c_0|^2 + 2*sum(|c_k|^2)
        E_freq = 2 * sum(abs(coeffs_A).^2);
        
        err_parseval = abs(E_time - E_freq) / (E_freq + 1e-20);
        
        fprintf('  Edge ID      : %d\n', edge_idx);
        fprintf('  Time Energy  : %.6e (Mean Square)\n', E_time);
        fprintf('  Freq Energy  : %.6e (Sum Squares * 2)\n', E_freq);
        fprintf('  Error        : %.2e %%\n', err_parseval * 100);
        
        if err_parseval < 1e-10
            fprintf('[PASS] Energy is conserved between domains.\n');
        else
            fprintf('[FAIL] Parseval check failed.\n');
        end
    end
    
    if exist(msh_file,'file'), try delete(msh_file); catch; end; end
end