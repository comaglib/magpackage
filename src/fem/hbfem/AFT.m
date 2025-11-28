classdef AFT < handle
    % AFT 谐波-时域转换模块 (Alternating Frequency Time)
    % 负责 FFT/IFFT 变换以及谐波索引管理
    
    properties
        Harmonics   % 谐波次数数组 [0, 1, 3, 5...]
        NumHarmonics % 谐波数量 K
        NumTimeSteps % 时域采样点数 N_t (通常 > 2*max(H)+1)
        
        % 变换矩阵 (预计算以加速)
        D_matrix    % [N_t x K] 频域到时域变换矩阵 (DFT逆变换)
        P_matrix    % [K x N_t] 时域到频域变换矩阵 (DFT正变换)
        
        % 基础频率
        BaseFreq
    end
    
    methods
        function obj = AFT(harmonics, n_time_steps, base_freq)
            obj.Harmonics = harmonics(:)'; % 行向量
            obj.NumHarmonics = length(harmonics);
            
            % 自动确定最少时间步数 (Nyquist)
            max_h = max(harmonics);
            min_steps = 2 * max_h + 1;
            
            if nargin < 2 || isempty(n_time_steps)
                % 默认取 2 的幂次以优化 FFT，或者至少满足 Nyquist
                obj.NumTimeSteps = max(32, 2^nextpow2(min_steps));
            else
                if n_time_steps < min_steps
                    warning('TimeSteps (%d) is insufficient for Harmonics max(%d). Adjusted to %d.', ...
                        n_time_steps, max_h, min_steps);
                    obj.NumTimeSteps = min_steps;
                else
                    obj.NumTimeSteps = n_time_steps;
                end
            end
            
            if nargin < 3, base_freq = 50; end
            obj.BaseFreq = base_freq;
            
            obj.precomputeMatrices();
        end
        
        function precomputeMatrices(obj)
            % 生成 DFT 和 IDFT 矩阵
            % 谐波表示: A(t) = sum( A_k_cos * cos(wkt) + A_k_sin * sin(wkt) )
            % 为了通用性，我们使用复数形式或者三角形式。
            % 这里的实现采用实数三角形式 (MagPackage v1.0 风格)，
            % 系数排列通常为: [DC, Cos1, Sin1, Cos3, Sin3 ...]
            
            % 但为了与标准 FEM 兼容，我们这里采用 **复数形式 (Complex Phasor)**
            % A_k 是复数，A(t) = Re( sum( A_k * exp(j*w*k*t) ) )
            % 注意：DC 分量通常是实数，但为了统一，我们统一处理为复数向量
            
            N = obj.NumTimeSteps;
            K = obj.NumHarmonics;
            H = obj.Harmonics;
            
            % 时间点 t_n = (0 : N-1) * T / N
            % theta_n = 2*pi * n / N
            theta = (0:N-1)' * (2*pi / N);
            
            % D Matrix (Freq -> Time)
            % f(t) = Re [ sum_{k} X_k * exp(j * h_k * theta) ]
            % 这里有一个系数问题。通常 FEM 中的 Phasor 是有效值还是峰值？
            % 假设 X_k 是峰值相量。
            % 特别注意 DC 分量 (h=0) 不需要 * 2 (如果展开的话)，但在复数形式下：
            % Real Signal f(t) <-> Two-sided spectrum.
            % 我们这里采用 **单边谱 (Single-Sided) + 复数** 的工程表示法。
            % f(t) = X_0 + sum_{k>0} Re( X_k * exp(j*w*t) ) * 2 ??? 
            % 为了简化，我们定义 X_k 为 "Complex Peak Value"。
            % A(t) = real( X_0 * exp(0) + X_1 * exp(j*w*t) + ... )
            % 这种定义下，X_k 的模就是波形的峰值。
            
            obj.D_matrix = zeros(N, K); % Complex matrix
            for i = 1:K
                h = H(i);
                % exp(j * h * theta)
                obj.D_matrix(:, i) = exp(1j * h * theta);
            end
            
            % P Matrix (Time -> Freq)
            % 利用正交性: (1/N) * sum( f(t) * exp(-j*h*theta) )
            % 对于 DC: factor = 1/N
            % 对于 AC: factor = 2/N (提取单边谱峰值)
            
            obj.P_matrix = zeros(K, N);
            for i = 1:K
                h = H(i);
                factor = 2.0 / N;
                if h == 0, factor = 1.0 / N; end
                
                % exp(-j * h * theta)
                obj.P_matrix(i, :) = factor * exp(-1j * h * theta').'; 
            end
        end
        
        function val_t = freq2time(obj, coeffs_freq)
            % 输入: coeffs_freq [K x M] (每列是一个变量的谐波系数)
            % 输出: val_t [N x M] (每列是时域波形)
            % 公式: Val = Real( D * Coeffs )
            
            % 矩阵乘法直接完成求和
            val_complex = obj.D_matrix * coeffs_freq;
            val_t = real(val_complex);
        end
        
        function coeffs_freq = time2freq(obj, val_t)
            % 输入: val_t [N x M]
            % 输出: coeffs_freq [K x M]
            % 公式: Coeffs = P * Val
            
            % 矩阵乘法完成积分/投影
            coeffs_freq = obj.P_matrix * val_t;
        end
        
        function dt = getTimeStep(obj)
            dt = (1.0 / obj.BaseFreq) / obj.NumTimeSteps;
        end
    end
end