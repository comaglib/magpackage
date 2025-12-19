function y = smooth_square_wave(t, freq, range, duty, tr, tf)
% SMOOTH_SQUARE_WAVE 生成带平滑过渡的矩形波 (C1 连续)
% 
% 输入参数:
%   t     - 时间向量
%   freq  - 频率 (Hz)
%   range - [V_low, V_high] 低电平和高电平的数值，例如 [0, 5] 或 [-1, 1]
%   duty  - 占空比 (0 到 100)，定义为上升沿开始到下降沿开始的时间比例
%   tr    - 上升时间 (s)
%   tf    - 下降时间 (s)
% 
% 输出:
%   y     - 波形数值
% 
% 示例:
%   t = linspace(0, 0.01, 1000);
%   y = smooth_square_wave(t, 100, [0, 5], 50, 1e-3, 1e-3);
%   plot(t, y);

    % 1. 基础参数计算
    T = 1 / freq;            % 周期
    w = (duty / 100) * T;    % 脉宽 (从上升沿开始到下降沿开始)
    V_low = range(1);
    V_high = range(2);
    H = V_high - V_low;      % 波幅高度

    % 2. 参数合法性检查
    if tr + tf > T
        error('错误: 上升时间 + 下降时间 超过了周期长度！');
    end
    if w < tr
        warning('警告: 占空比时间小于上升时间，波形可能无法达到高电平。');
    end
    if (T - w) < tf
        warning('警告: 关断时间小于下降时间，波形可能无法回到低电平。');
    end

    % 3. 将时间映射到单个周期内 [0, T)
    t_cycle = mod(t, T);

    % 4. 初始化输出
    y = zeros(size(t));

    % 5. 分段计算 (使用逻辑索引向量化操作)
    
    % --- A. 上升区 (Rising) ---
    % 范围: 0 <= t < tr
    % 归一化时间 x 从 0 -> 1
    idx_rise = (t_cycle >= 0) & (t_cycle < tr);
    if any(idx_rise)
        x = t_cycle(idx_rise) / tr;
        % Smoothstep 函数: 3x^2 - 2x^3
        scale = 3*x.^2 - 2*x.^3; 
        y(idx_rise) = V_low + H * scale;
    end

    % --- B. 高电平区 (High Stable) ---
    % 范围: tr <= t < w
    idx_high = (t_cycle >= tr) & (t_cycle < w);
    if any(idx_high)
        y(idx_high) = V_high;
    end

    % --- C. 下降区 (Falling) ---
    % 范围: w <= t < w + tf
    % 归一化时间 x 从 0 -> 1 (对应从高到低)
    idx_fall = (t_cycle >= w) & (t_cycle < w + tf);
    if any(idx_fall)
        x = (t_cycle(idx_fall) - w) / tf;
        % Smoothstep 函数
        scale = 3*x.^2 - 2*x.^3;
        % 从高向低过渡: V_high - H * scale
        y(idx_fall) = V_high - H * scale;
    end

    % --- D. 低电平区 (Low Stable) ---
    % 范围: w + tf <= t < T
    idx_low = (t_cycle >= w + tf);
    if any(idx_low)
        y(idx_low) = V_low;
    end
    
end