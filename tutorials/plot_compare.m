% load("C:\Users\haoha\code\magpackage\data\R-10-step-1e4-tol-1e3-2.mat")
% load("C:\Users\haoha\code\magpackage\data\R-10-step-1e3-tol-1e3-2.mat")
% 改名字为 'info_bdf_1e4' 和 'info_bdf_1e3'

figure;
subplot(2,1,1);
plot(info_bdf_1e4.times,info_bdf_1e4.CurrentHistory);
hold on;
plot(info_bdf_1e3.times,info_bdf_1e3.CurrentHistory);
legend('step-1e4-tol-1e3','step-1e3-tol-1e2');
grid on;
subplot(2,1,2);
plot(info_bdf_1e4.times,info_bdf_1e4.ProbeB_History);
hold on;
plot(info_bdf_1e3.times,info_bdf_1e3.ProbeB_History);
legend('step-1e4-tol-1e3','step-1e3-tol-1e2');
grid on;
subplot(2,1,1);
xlabel('Time (s)'); ylabel('Current (A)');
title(sprintf('Inrush Current'));
subplot(2,1,2);
xlabel('Time (s)'); ylabel('|B| at Origin (T)');
title(sprintf('B-Field at Core Center'));