% --- 验证脚本 ---
clc; clear;

% 定义函数
f = @(t) 1500 * sin(2 * pi * 50 * t) .* (t < 0.01);

% 生成测试时间序列 (0 到 0.02s)
t_vec = linspace(0, 0.02, 1000);
y_vec = f(t_vec);

% 绘图
figure('Color', 'w');
plot(t_vec, y_vec, 'LineWidth', 2, 'Color', 'b');
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Excitation Function Waveform');
xline(0.01, '--r', 'Cut-off at 0.01s'); % 标记切断点
axis tight;