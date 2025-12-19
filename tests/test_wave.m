% --- 测试脚本 ---
clear; clc; close all;

% 1. 设置参数
fs = 1e2;               % 采样率 1MHz
T_sim = 1.4;            % 仿真时长 2ms
t = 0 : 1/fs : T_sim;   % 时间向量

freq = 0.5;             % 信号频率 1kHz
levels = [0, 100];      % 电压范围 -5V 到 5V
duty = 40;              % 占空比 50%
tr = 0.2;               % 上升时间 100us
tf = 0.2;               % 下降时间 100us

% 2. 生成波形
y = smooth_square_wave(t, freq, levels, duty, tr, tf);

% 3. 计算导数 (验证连续性)
dy = gradient(y) ./ gradient(t);

% 4. 绘图
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

subplot(2,1,1);
plot(t*1e3, y, 'LineWidth', 2);
grid on;
title('生成的 C1 连续矩形波');
ylabel('幅值 (V)');
ylim([levels(1)-1, levels(2)+1]);

subplot(2,1,2);
plot(t*1e3, dy, 'r', 'LineWidth', 1.5);
grid on;
title('一阶导数 (dy/dt)');
xlabel('时间 (ms)');
ylabel('斜率 (V/s)');

% 局部放大查看过渡区
axes('Position',[.6 .7 .25 .15])
box on
idx_zoom = (t < tr*1.5);
plot(t(idx_zoom)*1e3, y(idx_zoom), 'LineWidth', 2);
title('上升沿局部放大 (Smoothstep)');
grid on;