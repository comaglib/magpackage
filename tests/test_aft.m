% tests/test_aft.m
clear; clc;
addpath(genpath('src'));

fprintf('==============================================\n');
fprintf('   Test: AFT Module (Time-Freq Conversion)    \n');
fprintf('==============================================\n');

% --- 1. 配置 AFT ---
base_freq = 50;
omega = 2 * pi * base_freq;
harmonics = [0, 1, 3, 5]; % [Fix] 加入 DC 分量以全面测试
aft = AFT(harmonics, [], base_freq); 

fprintf('[Setup] Harmonics: %s\n', mat2str(harmonics));
fprintf('[Setup] Time Steps: %d\n', aft.NumTimeSteps);

% --- 2. 构造已知信号 ---
% f(t) = 5(DC) + 10*cos(wt) + 4*sin(3wt) + 2*cos(5wt + pi/4)
% 对应相量 (Single-Sided Peak Phasor):
% h=0: 5
% h=1: 10
% h=3: 4*sin(3wt) = 4 * real( -1j * exp(j3wt) ) -> Coeff = -4j
% h=5: 2*cos(5wt + pi/4) = 2 * exp(j*pi/4)

t = (0 : aft.NumTimeSteps-1) * aft.getTimeStep();
f_target = 5 + ...
           10 * cos(omega * t) + ...
           4 * sin(3 * omega * t) + ...
           2 * cos(5 * omega * t + pi/4);

f_target = f_target(:);

% --- 3. 测试: 时域 -> 频域 (Analysis) ---
fprintf('\n[Test 1] Time to Frequency (Analysis)...\n');
coeffs_calc = aft.time2freq(f_target);

% 预期系数
coeffs_expected = zeros(aft.NumHarmonics, 1);
coeffs_expected(1) = 5;                 % h=0
coeffs_expected(2) = 10;                % h=1
coeffs_expected(3) = -4j;               % h=3
coeffs_expected(4) = 2 * exp(1j*pi/4);  % h=5

fprintf('  H | Expected (Mag/Ang) | Calculated (Mag/Ang) | Error \n');
fprintf('-------------------------------------------------------\n');
for i = 1:aft.NumHarmonics
    h = harmonics(i);
    exp_val = coeffs_expected(i);
    calc_val = coeffs_calc(i);
    err = abs(calc_val - exp_val);
    
    fprintf('  %d | %6.3f / %6.1f   | %6.3f / %6.1f     | %.2e\n', ...
        h, abs(exp_val), rad2deg(angle(exp_val)), ...
        abs(calc_val), rad2deg(angle(calc_val)), err);
end

if norm(coeffs_calc - coeffs_expected) < 1e-10
    fprintf('  [PASS] Frequency coefficients match expected phasors.\n');
else
    fprintf('  [FAIL] Coefficients mismatch.\n');
end

% --- 4. 测试: 频域 -> 时域 (Synthesis) ---
fprintf('\n[Test 2] Frequency to Time (Synthesis)...\n');
f_reconstructed = aft.freq2time(coeffs_calc);

err_time = norm(f_reconstructed - f_target) / norm(f_target);
fprintf('  - Time Domain Reconstruction Error: %.2e\n', err_time);

if err_time < 1e-10
    fprintf('  [PASS] Signal reconstructed successfully.\n');
else
    fprintf('  [FAIL] Signal reconstruction failed.\n');
end

% --- 5. 验证变换矩阵性质 (Check Scaling) ---
fprintf('\n[Test 3] Matrix Product (P * D)...\n');
% 理论: P * D 应该是对角矩阵
% DC (h=0) 对应 1.0
% AC (h>0) 对应 2.0 (因为我们提取的是单边谱峰值)
Product = aft.P_matrix * aft.D_matrix;

% 构建预期对角矩阵
Expected_Diag = zeros(aft.NumHarmonics, 1);
for i = 1:aft.NumHarmonics
    if harmonics(i) == 0
        Expected_Diag(i) = 1.0;
    else
        Expected_Diag(i) = 2.0;
    end
end
Ideal_Mat = diag(Expected_Diag);

ortho_err = norm(Product - Ideal_Mat, 'fro');
fprintf('  - Difference from Ideal Scaling |PD - diag([1,2..])|: %.2e\n', ortho_err);

if ortho_err < 1e-10
    fprintf('  [PASS] Transform matrices have correct Engineering Scaling.\n');
else
    fprintf('  [FAIL] Matrix scaling incorrect.\n');
    disp('Computed P*D diagonal:');
    disp(diag(Product).');
end