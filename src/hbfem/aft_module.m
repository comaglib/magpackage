function AFT = aft_module()
% AFT_MODULE 谐波平衡法核心模块 (修正版)
    AFT.prepare_indices = @prepare_indices;
    AFT.freq2time = @freq2time;
    AFT.time2freq = @time2freq;
end

function Info = prepare_indices(Harmonics, N_time_steps)
    if N_time_steps < 2 * max(Harmonics) + 1
        warning('AFT: N_time_steps (%d) is too small.', N_time_steps);
    end
    
    Info.Harmonics = Harmonics;
    Info.NumHarmonics = length(Harmonics);
    Info.Nt = N_time_steps;
    
    Info.PosIndices = zeros(1, Info.NumHarmonics);
    Info.NegIndices = zeros(1, Info.NumHarmonics);
    
    for k = 1:Info.NumHarmonics
        h = Harmonics(k);
        % MATLAB FFT: Index 1 is DC, Index 2 is 1st harmonic...
        Info.PosIndices(k) = h + 1;
        
        if h == 0
            % DC (0Hz) 特殊处理：映射到自身 (Index 1)
            % 在 freq2time 中，DC 不需要共轭镜像，或者镜像位置仍指向 1
            Info.NegIndices(k) = 1; 
        else
            Info.NegIndices(k) = N_time_steps - h + 1;
        end
    end
end

function val_t = freq2time(coeffs, Info)
    [M, ~] = size(coeffs);
    Nt = Info.Nt;
    fft_vec = zeros(M, Nt);
    
    idx_pos = Info.PosIndices;
    idx_neg = Info.NegIndices;
    
    % 填充频谱
    fft_vec(:, idx_pos) = coeffs;
    
    % 填充负频率 (共轭)，保证时域为实数
    % 对于 DC (h=0)，idx_pos=1, idx_neg=1。
    % 重复赋值无害，conj(Real_DC) = Real_DC。
    % 注意：通常 DC 分量在 coeffs 里存的是全幅值，不需要除以2再镜像。
    % 这里假设 coeffs 存储的是单边谱 (Single-Sided) 的复振幅 X_k (k>0)
    % 对于 k=0 (DC)，直接放入 Index 1。
    
    % 区分处理：
    for k = 1:Info.NumHarmonics
        h = Info.Harmonics(k);
        idx_p = Info.PosIndices(k);
        idx_n = Info.NegIndices(k);
        
        c = coeffs(:, k);
        
        if h == 0
            % DC 分量直接赋值
            fft_vec(:, idx_p) = c;
        else
            % AC 分量：假设输入是双边谱系数的一半 (或单边谱系数)
            % 这里的定义要与 time2freq 匹配。
            % time2freq 输出的是 fft()/N。
            % 那么正频率处是 c_k, 负频率处是 conj(c_k).
            % 只有这样 ifft(sum c_k e^...) * N 才能还原 x(t).
            fft_vec(:, idx_p) = c;
            fft_vec(:, idx_n) = conj(c);
        end
    end
    
    val_t = ifft(fft_vec, Nt, 2, 'symmetric') * Nt; 
end

function coeffs = time2freq(val_t, Info)
    [~, Nt] = size(val_t);
    fft_vec = fft(val_t, [], 2) / Nt; % Normalize by N
    idx_pos = Info.PosIndices;
    coeffs = fft_vec(:, idx_pos);
end