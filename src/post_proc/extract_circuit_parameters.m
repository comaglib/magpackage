function CircuitParams = extract_circuit_parameters(Model, Solution, Coil, freq)
% EXTRACT_CIRCUIT_PARAMETERS 提取线圈等效电路参数 (R, L, Z)
% 基于能量法: 
%   P_loss = 1/2 * R * |I|^2  => R
%   W_mag  = 1/4 * L * |I|^2  => L (线性近似)
%
% 输入:
%   Model, Solution: 求解结果
%   Coil: 线圈定义 (用于获取电流幅值 I)
%   freq: 工作频率 (Hz)
%
% 输出:
%   CircuitParams: 结构体 (.R, .L, .Z, .X)

    fprintf('  [Post] Extracting Circuit Parameters (Energy Method)...\n');
    
    % 1. 获取电流幅值
    % 假设 Coil.I 是相量数组，取第一个元素的模
    I_peak = abs(Coil.I(1));
    if I_peak < 1e-12
        warning('Coil current is zero. Cannot extract parameters.');
        CircuitParams = struct('R',0,'L',0,'Z',0);
        return;
    end
    
    % 2. 计算全域能量与损耗
    % 欧姆损耗 (平均值)
    % 注意: compute_ohmic_loss 内部已经使用了 1/2 系数 (针对峰值输入)
    A = Solution.A;
    V = Solution.V_active; % 这里假设 solve_frequency_linear 返回的是 V_active
    % 如果 Solution.V 是全节点，则需要截取? 
    % 不，solve_frequency_linear 中 Solution.V = x(numEdges+1:end)，这正是 V_active
    % 所以直接传入即可。
    
    [P_loss_avg, ~] = compute_ohmic_loss(Model, A, V, [], freq);
    
    % 磁场储能 (平均值)
    W_mag_avg = compute_magnetic_energy(Model, A, Coil);
    
    % 3. 导出参数
    % P_avg = 0.5 * R_ac * I_peak^2
    R_ac = (2 * P_loss_avg) / (I_peak^2);
    
    % W_avg = 0.25 * L_ac * I_peak^2  (注意系数 1/4)
    % 或者用无功功率 Q = 2 * omega * (Wm - We) = 0.5 * X * I^2
    % X = 4 * omega * Wm / I^2
    % L = X / omega = 4 * Wm / I^2
    
    L_ac = (4 * W_mag_avg) / (I_peak^2);
    
    omega = 2 * pi * freq;
    X_ac = omega * L_ac;
    Z_ac = R_ac + 1i * X_ac;
    
    CircuitParams.R = R_ac;
    CircuitParams.L = L_ac;
    CircuitParams.X = X_ac;
    CircuitParams.Z = Z_ac;
    
    fprintf('    Freq = %.1f Hz\n', freq);
    fprintf('    Current (Peak) = %.2f A\n', I_peak);
    fprintf('    Total Loss     = %.4e W\n', P_loss_avg);
    fprintf('    Stored Energy  = %.4e J\n', W_mag_avg);
    fprintf('    ------------------------\n');
    fprintf('    R_ac = %.4e Ohm\n', R_ac);
    fprintf('    L_ac = %.4e H\n', L_ac);
    fprintf('    Z_ac = %.4e + j%.4e Ohm\n', real(Z_ac), imag(Z_ac));
end