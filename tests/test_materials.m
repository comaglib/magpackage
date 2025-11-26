function test_materials()
% TEST_MATERIALS 验证材料非线性处理模块

    addpath(genpath('src'));
    fprintf('Starting Material Module Test...\n');
    
    % 1. 定义一个典型的软磁 B-H 曲线 (Froehlich模型近似)
    % B = H / (a + b H)  => H = a B / (1 - b B)
    % 设 Saturation Bs = 2.0 T, Initial Mu_r = 1000
    mu0 = 4*pi*1e-7;
    mu_i = 1000 * mu0;
    Bs = 2.0;
    
    % 生成数据点
    B_pts = linspace(0, 1.9, 50); % 不要碰到极点
    H_pts = B_pts ./ (mu_i * (1 - B_pts/Bs));
    
    Model.Materials.Lib(1).Name = 'TestIron';
    Model.Materials.Lib(1).BH_Curve.B = B_pts;
    Model.Materials.Lib(1).BH_Curve.H = H_pts;
    
    % 2. 执行预处理
    Model = preprocess_materials(Model);
    MatInfo = Model.Materials.Lib(1);
    
    % 3. 测试查询
    B_test = [0.1, 1.0, 1.5, 2.5 1000]'; % 2.5 1000 是外推区域
    B_sq = B_test.^2;
    
    [nu, dnu] = eval_material_nu(B_sq, MatInfo);
    
    fprintf('\nTesting BH Evaluation:\n');
    fprintf('%10s | %10s | %10s | %10s\n', 'B (T)', 'nu (m/H)', 'dnu/dB2', 'mu_r_eff');
    for i = 1:length(B_test)
        mu_eff = 1 / (nu(i) * mu0);
        fprintf('%10.2f | %10.1f | %10.1f | %10.1f\n', ...
            B_test(i), nu(i), dnu(i), mu_eff);
    end
    
    % 4. 验证导数一致性 (数值差分)
    delta = 1e-5;
    B2_center = 1.0^2;
    [nu_c, dnu_c] = eval_material_nu(B2_center, MatInfo);
    [nu_p, ~]     = eval_material_nu(B2_center + delta, MatInfo);
    [nu_m, ~]     = eval_material_nu(B2_center - delta, MatInfo);
    
    dnu_num = (nu_p - nu_m) / (2 * delta);
    
    err = abs(dnu_c - dnu_num) / abs(dnu_c);
    fprintf('\nDerivative Check (at B=1T):\n');
    fprintf('  Analytical dnu: %e\n', dnu_c);
    fprintf('  Numerical  dnu: %e\n', dnu_num);
    fprintf('  Rel Error:      %e\n', err);
    
    if err < 1e-3
        fprintf('[PASS] Material derivative is consistent.\n');
    else
        error('[FAIL] Derivative computation inaccurate.');
    end
end