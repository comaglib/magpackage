% run_ecore_dftfem.m - E-Core Transformer DFT-FEM Benchmark
%
% 描述:
%   该脚本使用离散傅里叶变换有限元法 (DFTFEM) 求解 E 型变压器的稳态响应。
%   相比 HBFEM，DFTFEM 求解全频谱（连续谐波），矩阵更稠密，适合强非线性
%   或波形畸变严重的情况。
%
% 配置变更:
%   1. 求解器: 使用 DFTFEMCoupledSolver。
%   2. 谐波: 使用 0~7 次全连续谐波，以展示全频谱求解能力。
%   3. 内存: 增加了 MUMPS 内存分配以应对稠密耦合矩阵。

clear; clc; tic;
fprintf('=========================================================\n');
fprintf('   E-Core DFT-FEM Analysis (Full Spectrum 0-7)           \n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格加载 ---
% (与 run_ecore_hbfem 保持一致)
meshFile = 'data/meshData/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

tags_primary   = [5, 6, 8, 9];
tags_secondary = [3, 4, 7, 10];
tag_core       = 2;
tag_air        = 1;
NEW_TAG_PRIM = 50;
NEW_TAG_SEC  = 60;

mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;
mesh.RegionTags(ismember(mesh.RegionTags, tags_secondary)) = NEW_TAG_SEC;

%% --- 2. 材料与物理属性配置 ---
% (与 run_ecore_hbfem 保持一致)
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 2.1 铁芯材料 (非线性)
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

% 2.2 线性材料
matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;
matLib(NEW_TAG_SEC)  = matLinear;

% 2.3 有限元空间
space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

%% --- 3. 线圈设置 ---
% (与 run_ecore_hbfem 保持一致)
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);
winding = Winding('Primary', NEW_TAG_PRIM, 300, 100, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

%% --- 4. DFT-FEM 求解器配置 ---
assembler = Assembler(mesh, dofHandler);

% 4.1 定义全频谱谐波
% DFTFEM 的特点是计算连续频谱。这里取 0 到 7 次谐波。
% 0=DC, 1=Fund, 2,3...
harmonics = 0:7;
baseFreq = 50;
n_time_steps = 64;
aft = AFT(harmonics, n_time_steps, baseFreq);

fprintf('   -> Harmonics: Full Spectrum [0 .. 7] (Total %d)\n', length(harmonics));

R_circuit = 100;

% 4.2 创建 DFTFEM 求解器
solver = DFTFEMCoupledSolver(assembler, aft, winding, R_circuit);
solver.Tolerance = 1e-3;
solver.MaxIter = 50;

% [关键配置]
% 1. 非对称模式: DFTFEM 耦合矩阵高度非对称且复数。
% 2. 增加内存: 全频谱耦合导致矩阵非零元大幅增加。
solver.LinearSolver.MumpsSymmetry = 0; 
solver.LinearSolver.MumpsICNTL.i14 = 400; % 400% 内存预分配

% 4.3 定义电压激励 (全频谱向量)
% V(t) = 76 * sin(wt) -> 相量 V = -76j (在基频 h=1 处)
% 注意：harmonics = [0, 1, 2, ...]，因此 h=1 对应索引 2
idx_fund = find(harmonics == 1);
V_spectrum = zeros(length(harmonics), 1);
V_spectrum(idx_fund) = -76j; 

fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

% 4.4 动态正则化 (与 HBFEM 保持一致以确保 Sigma=0 稳定)
fprintf('[Auto-Scale] Estimating Stiffness Norm...\n');
nu_vec = ones(mesh.NumElements, 1) * (1/(4*pi*1e-7)); 
K_sample = assembler.assembleStiffness(space_A, nu_vec);
K_norm = norm(K_sample, 1);
M_reg = assembler.assembleMass(space_A);
M_norm = norm(M_reg, 1);
target_ratio = 1e-7; 
epsilon = (K_norm * target_ratio) / M_norm;
fprintf('   -> Calculated Robust Epsilon = %.4e\n', epsilon);
solver.MatrixK_Add = epsilon * M_reg;

%% --- 5. 执行求解 ---
fprintf('[Step 5] Solving DFT-FEM...\n');
try
    [X_sol, I_sol, info] = solver.solve(space_A, matLib, V_spectrum, fixedDofs_A);
    
    % 结果检查
    I_fund_mag = abs(I_sol(idx_fund));
    if I_fund_mag > 1e10 || isnan(I_fund_mag)
        error('Solution Exploded. Matrix is singular.');
    elseif I_fund_mag < 1e-6
        warning('Solution is Zero.');
    end
    
catch ME
    fprintf('Error during solve: %s\n', ME.message);
    return;
end

%% --- 6. 后处理与可视化 ---
fprintf('[Step 6] Post-Processing...\n');

% 6.1 重构时域电流波形
num_time_pts = 100;
t_cycle = linspace(0, 1/baseFreq, num_time_pts);
I_time = zeros(size(t_cycle));

for k = 1:length(harmonics)
    h = harmonics(k);
    omega = 2 * pi * baseFreq * h;
    phasor = I_sol(k);
    I_time = I_time + real(phasor * exp(1j * omega * t_cycle));
end

% 6.2 打印前几个主要谐波
fprintf('\n--- Current Harmonics (Top Components) ---\n');
for k = 1:length(harmonics)
    h = harmonics(k);
    % 仅打印幅值显著的分量
    if abs(I_sol(k)) > 1e-3
        fprintf('   Order %d:  %.4f A  (Phase: %.1f deg)\n', ...
            h, abs(I_sol(k)), angle(I_sol(k)) * 180/pi);
    end
end

% 6.3 绘图
figure('Name', 'E-Core DFT-FEM Analysis', 'Position', [150, 150, 1000, 500]);

% 子图1: 电流波形
subplot(1, 2, 1);
plot(t_cycle*1000, I_time, 'r-', 'LineWidth', 2); grid on;
xlabel('Time (ms)'); ylabel('Current (A)');
title('Steady-State Current (Full Spectrum)');
legend('DFT-FEM Current');

% 子图2: 峰值时刻 B 场
[~, max_idx] = max(abs(I_time));
t_peak = t_cycle(max_idx);

A_peak = zeros(dofHandler.NumGlobalDofs, 1);
for k = 1:length(harmonics)
    h = harmonics(k);
    phasor = exp(1j * 2*pi*baseFreq*h * t_peak);
    A_peak = A_peak + real(X_sol(:, k) * phasor);
end

post = PostProcessor(assembler);
B_elems = post.computeElementB(A_peak, 'Nedelec_P1');
B_mag = sqrt(sum(abs(B_elems).^2, 1));
B_nodes = post.mapElementsToNodes(B_mag);
viz = Visualizer(post);
B_core = post.combineRegionData(B_nodes, tag_core);

subplot(1, 2, 2);
viz.plotFieldOnSurface(tag_core, B_core(:), 'FaceAlpha', 1.0, 'EdgeColor', 'none');
title(sprintf('|B| at Peak (t=%.2f ms)', t_peak*1000));
colorbar; view(3); axis equal;

toc;