classdef AFT < handle
    % AFT 谐波-时域转换模块 (Alternating Frequency Time)
    % 负责 FFT/IFFT 变换以及谐波索引管理
    
    properties
        Harmonics   % 谐波次数数组 [0, 1, 3, 5...]
        NumHarmonics % 谐波数量 K
        NumTimeSteps % 时域采样点数 N_t (通常 > 2*max(H)+1)
        
        % 变换矩阵 (预计算以加速)
        D_matrix    % [N_t x K] 频域到时域变换矩阵 (IDFT)
        P_matrix    % [K x N_t] 时域到频域变换矩阵 (DFT)
        
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
                if ~isscalar(n_time_steps)
                    % [修复] 如果 n_time_steps 是向量 (参数错误)，仅使用第一个元素进行警告，并使用默认值
                    ts_scalar = n_time_steps(1);
                    warning('AFT:InputError', 'TimeSteps argument is a vector (%.0f,...). Check argument order. Using default steps.', ts_scalar);
                    obj.NumTimeSteps = max(32, 2^nextpow2(min_steps));
                elseif n_time_steps < min_steps
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
            N = obj.NumTimeSteps;
            K = obj.NumHarmonics;
            H = obj.Harmonics;
            
            theta = (0:N-1)' * (2*pi / N);
            
            % D Matrix (Freq -> Time)
            obj.D_matrix = zeros(N, K); % Complex matrix
            for i = 1:K
                h = H(i);
                obj.D_matrix(:, i) = exp(1j * h * theta);
            end
            
            % P Matrix (Time -> Freq)
            obj.P_matrix = zeros(K, N);
            for i = 1:K
                h = H(i);
                factor = 2.0 / N;
                if h == 0, factor = 1.0 / N; end
                
                obj.P_matrix(i, :) = factor * exp(-1j * h * theta').'; 
            end
        end
        
        function val_t = freq2time(obj, coeffs_freq)
            val_complex = obj.D_matrix * coeffs_freq;
            val_t = real(val_complex);
        end
        
        function coeffs_freq = time2freq(obj, val_t)
            coeffs_freq = obj.P_matrix * val_t;
        end
        
        function dt = getTimeStep(obj)
            dt = (1.0 / obj.BaseFreq) / obj.NumTimeSteps;
        end
    end
end