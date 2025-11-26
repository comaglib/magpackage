function test_biot_savart()
% TEST_BIOT_SAVART 验证解析线圈源计算模块
% 1. 生成跑道型线圈
% 2. 计算轴线上的 A 场分布 (应为 0, 因为 A 沿切向)
% 3. 计算偏轴点的 A 场，验证对称性

    addpath(genpath('src'));
    fprintf('Starting Biot-Savart Source Test...\n');
    
    % 1. 定义线圈 (圆形: Length=0)
    Params.Center = [0, 0, 0];
    Params.Normal = [0, 0, 1];
    Params.Length = 0;      % 纯圆
    Params.Radius = 1.0;    % R = 1m
    Params.Current = 1000;  % I = 1000A
    Params.N_seg = 72;      % 离散精度
    
    Coil = create_racetrack_coil(Params);
    
    fprintf('  - Generated Coil: %d segments\n', size(Coil.P1, 2));
    
    % 2. 定义测试点
    % Case A: Z轴上的点 (0,0,z)
    % 对于圆环，A 在轴线上应为 0 (因为环形电流产生的 A 只有 phi 分量，轴线上 r=0)
    z = linspace(-1, 1, 5);
    Pts_Axis = [zeros(1,5); zeros(1,5); z];
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    A_axis = compute_biot_savart_A(Coil, Pts_Axis);
    
    norm_A_axis = sqrt(sum(A_axis.^2, 1));
    fprintf('  - On-Axis Field Norm (Expect ~0): %e\n', max(norm_A_axis));
    
    if max(norm_A_axis) < 1e-6
        fprintf('[PASS] Vector potential is zero on the symmetry axis.\n');
    else
        warning('[WARN] Non-zero potential on axis (Discretization error?).');
    end
    
    % Case B: XY平面上的点 (R/2, 0, 0) 和 (-R/2, 0, 0)
    % 对称性：Ay 应相反，Ax 应为 0 (电流主要沿 Phi)
    Pts_Sym = [0.5, -0.5; 0, 0; 0, 0]; 
    A_sym = compute_biot_savart_A(Coil, Pts_Sym);
    
    disp('  - Symmetry Check Points:');
    disp(A_sym);
    
    % 检查 Ay 的符号反转
    if abs(A_sym(2,1) + A_sym(2,2)) < 1e-6 && abs(A_sym(2,1)) > 1e-6
        fprintf('[PASS] Antisymmetry verified (Ay(x) = -Ay(-x)).\n');
    else
        warning('[FAIL] Symmetry check failed.');
    end
    
    % 3. 性能测试 (模拟大网格)
    fprintf('  - Performance Benchmark (100k points)...\n');
    Pts_Large = rand(3, 100000);
    tic;
    compute_biot_savart_A(Coil, Pts_Large);
    t_bench = toc;
    fprintf('    100k points calculated in %.4f s\n', t_bench);
    
    fprintf('Biot-Savart Test Finished.\n');
end