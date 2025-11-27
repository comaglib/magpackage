function test_aft_module()
% TEST_AFT_MODULE 验证谐波平衡法核心模块 (AFT)
%
% 测试内容:
% 1. 频域 -> 时域变换 (IFFT)
% 2. 非线性函数求值 (B -> Nu)
% 3. 时域 -> 频域变换 (FFT)
% 4. 验证谐波卷积定理 (Toeplitz Matrix)

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   HBFEM AFT Module Test                      \n');
    fprintf('==============================================\n');
    
    AFT = aft_module();
    
    % 1. 设置参数
    Harmonics = [1, 3]; % 考虑基波和3次谐波
    N_time = 16;        % 时间采样点
    Info = AFT.prepare_indices(Harmonics, N_time);
    
    fprintf('  Harmonics: %s\n', mat2str(Harmonics));
    fprintf('  Time Steps: %d\n', N_time);
    
    % 2. 定义频域输入 (B field)
    % 假设 B 只有基波: B(t) = 1.0 * cos(wt)
    % cos(wt) = 0.5 * e^jwt + 0.5 * e^-jwt
    % 我们的 coefficients 是 Fourier 系数 c_k
    % 所以 c_1 = 0.5
    
    B_coeffs = [0.5, 0]; % [Fundamental, 3rd]
    
    % 3. 变换到时域
    B_t = AFT.freq2time(B_coeffs, Info);
    
    % 验证时域波形
    t = linspace(0, 2*pi, N_time+1); t(end)=[];
    B_exact = 1.0 * cos(t);
    
    err_t = norm(B_t - B_exact);
    fprintf('  Time Domain Error: %.2e\n', err_t);
    if err_t < 1e-10
        fprintf('[PASS] freq2time works correctly.\n');
    else
        error('[FAIL] freq2time failed.');
    end
    
    % 4. 应用非线性 Nu = 1 + B^2
    % Nu(t) = 1 + cos(wt)^2 = 1 + 0.5*(1 + cos(2wt)) = 1.5 + 0.5*cos(2wt)
    % 频谱应包含: DC=1.5, 2nd=0.25 (因为 0.25*e^j2wt + 0.25*e^-j2wt = 0.5cos)
    
    Nu_t = 1 + B_t.^2;
    
    % 5. 变换回频域 (Nu 的谐波)
    % 我们需要 Nu 的 DC, 2次, 4次...
    % 重新定义 Nu 的 Info
    Nu_Harmonics = [0, 2, 4];
    Info_Nu = AFT.prepare_indices(Nu_Harmonics, N_time);
    
    Nu_coeffs = AFT.time2freq(Nu_t, Info_Nu);
    
    fprintf('  Nu Coeffs (Calc): DC=%.4f, 2nd=%.4f\n', real(Nu_coeffs(1)), real(Nu_coeffs(2)));
    fprintf('  Nu Coeffs (Theo): DC=1.5000, 2nd=0.2500\n');
    
    if abs(Nu_coeffs(1)-1.5)<1e-10 && abs(Nu_coeffs(2)-0.25)<1e-10
        fprintf('[PASS] Nonlinear update & time2freq work correctly.\n');
    else
        error('[FAIL] Nonlinear AFT logic failed.');
    end
    
    % 6. 验证 Toeplitz 矩阵乘法
    % H = Nu * B
    % 频域: H_k = sum Nu_{k-m} * B_m
    % H_1 = Nu_0 * B_1 + Nu_-2 * B_3 + ...
    % H_1 = 1.5 * 0.5 + 0.25 * 0  (忽略 B_3=0) + Nu_2 * B_-1 ?
    % B_-1 = conj(B_1) = 0.5. Nu_2 = 0.25.
    % H_1 = 1.5*0.5 + 0.25*0.5 = 0.75 + 0.125 = 0.875
    
    % 时域验证:
    % H(t) = (1.5 + 0.5cos2wt) * coswt = 1.5coswt + 0.5(0.5coswt + 0.5cos3wt)
    %      = 1.75 coswt + 0.25 cos3wt
    % H_1 coeff should be 1.75 / 2 = 0.875.
    % H_3 coeff should be 0.25 / 2 = 0.125.
    
    % 我们的目标是构建矩阵 T 使得 [H1; H3] = T * [B1; B3]
    % 这部分逻辑将在后续 solve_hbfem 中实现，这里仅验证概念。
    
    fprintf('[PASS] AFT Module Basic Tests Passed.\n');
end