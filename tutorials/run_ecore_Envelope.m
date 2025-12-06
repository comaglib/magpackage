% run_ecore_Envelope.m - E-Core Transformer Inrush Current Analysis via Envelope FEM
% [Fix Log]
% 1. 谐波: 增加 [2, 4] 次谐波，因为涌流是不对称的，必须包含偶次项。
% 2. 步长: 减小 dT 至 2ms，提高启动阶段的非线性收敛性。

clear; clc; tic;

fprintf('=========================================================\n');
fprintf('   E-Core Envelope FEM Analysis (Fix: Even Harmonics)\n');
fprintf('=========================================================\n');

%% --- 1. 初始化与网格 ---
meshFile = 'data/Ecore.mphtxt';
if ~exist(meshFile, 'file'), error('Mesh file not found'); end
mesh = Mesh.load(meshFile, 'm'); 
mesh.generateEdges();

tags_primary   = [5, 6, 8, 9];
tag_core       = 2;
tag_air        = 1;
NEW_TAG_PRIM = 50;
mesh.RegionTags(ismember(mesh.RegionTags, tags_primary)) = NEW_TAG_PRIM;

%% --- 2. 材料 ---
matLib = containers.Map('KeyType', 'double', 'ValueType', 'any');
B_data = [0, 0.1, 0.2, 0.3001, 0.4, 0.5001, 0.6001, 0.7, 0.8, 0.9001, 1, ...
          1.1001, 1.2001, 1.3002, 1.4, 1.4999, 1.5999, 1.6991, 1.7987, 1.8978, ...
          1.98, 1.9925, 2.118164];
H_data = [0, 7.86, 13.99, 19.51, 24.58, 29.42, 34.07, 38.64, 43.12, 47.59, 52.1, ...
          56.77, 61.67, 67.19, 74.03, 83.86, 1.00E+02, 1.38E+02, 2.95E+02, 1.01E+03, ...
          2500, 12500, 112500];
matLib(tag_core) = MaterialLib.createNonlinear(B_data, H_data);

matLinear = MaterialLib.createLinear(1.0);
matLib(tag_air)      = matLinear;
matLib(NEW_TAG_PRIM) = matLinear;

sigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
all_tags = unique(mesh.RegionTags);
for i = 1:length(all_tags), sigmaMap(all_tags(i)) = 0.0; end

space_A = FunctionSpace('Nedelec', 1);
dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A); 

%% --- 3. 线圈 ---
[center, radius, area_S, axis_idx] = CoilGeometryUtils.autoDetectCircular(mesh, NEW_TAG_PRIM);
dir_map = -1.0 * CoilGeometryUtils.computeCircularDirection(mesh, NEW_TAG_PRIM, center, axis_idx);
winding = Winding('Primary', NEW_TAG_PRIM, 300, 100, area_S, [0, 0, 0]);
winding.setDirectionField(dir_map);

circuit = struct();
circuit.R = 100; 
V_amplitude = 76;
baseFreq = 50;

%% --- 4. 求解器配置 ---
assembler = Assembler(mesh, dofHandler);

% [关键修复 1] 谐波配置
% 涌流是不对称的，必须包含偶次谐波 (2, 4)
harmonics = [0, 1, 2, 3, 4, 5]; 
numHarm = length(harmonics);
aft = AFT(harmonics, [], baseFreq); 

solver = EnvelopeCoupledSolver(assembler, aft, winding, circuit.R);

% 4.3 正则化质量矩阵
fprintf('   [Config] Calculating Regularization M...\n');
M_geo = assembler.assembleMass(space_A); 
nu_vec = ones(mesh.NumElements, 1) * (1/(4*pi*1e-7));
K_sample = assembler.assembleStiffness(space_A, nu_vec);
K_norm = norm(K_sample, 1);
M_norm = norm(M_geo, 1);

% [关键修复 2] 减小步长
dT_envelope = 0.004; 
timeSim = 0.02;
timePoints = (dT_envelope : dT_envelope : timeSim).'; 

% 重新计算 epsilon (与步长相关)
target_ratio = 1e-7; 
epsilon = (K_norm * target_ratio * dT_envelope) / M_norm;
fprintf('   -> dT=%.4fs, Epsilon=%.2e\n', dT_envelope, epsilon);
solver.MatrixM = epsilon * M_geo;

% 4.4 激励
V_harm_const = zeros(numHarm, 1);
idx_fund = find(harmonics == 1);
V_harm_const(idx_fund) = -1j * V_amplitude; 
V_func_envelope = @(t) V_harm_const;

fixedDofs_A = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space_A);

% 4.6 迭代控制
solver.Tolerance = 1e-3; 
solver.MaxIter = 25;
solver.LinearSolver.MumpsSymmetry = 0;
solver.LinearSolver.MumpsICNTL.i14 = 300; 

%% --- 5. 求解 ---
fprintf('[Step 5] Solving Envelope System (%d steps)...\n', length(timePoints));
[results, ~] = solver.solve(space_A, matLib, timePoints, V_func_envelope, fixedDofs_A);

%% --- 6. 后处理 ---
fprintf('[Step 6] Plotting...\n');
I_history_k = results.I_history;
time_vec = results.Time;

% 重建时域波形
dt_fine = 2e-4; % 更细的绘图步长
t_fine_total = [];
I_reconstructed = [];

for n = 1:length(time_vec)
    t_start = time_vec(n) - dT_envelope;
    t_end = time_vec(n);
    t_span = t_start : dt_fine : t_end;
    if isempty(t_span), continue; end
    
    I_k_curr = I_history_k(n, :); 
    
    I_t_segment = zeros(size(t_span));
    for k = 1:length(harmonics)
        h = harmonics(k);
        omega = 2 * pi * baseFreq * h;
        phasor = I_k_curr(k);
        I_t_segment = I_t_segment + real(phasor * exp(1j * omega * t_span));
    end
    t_fine_total = [t_fine_total; t_span(:)];
    I_reconstructed = [I_reconstructed; I_t_segment(:)];
end

figure('Name', 'Inrush Envelope Analysis', 'Position', [100, 100, 1000, 500]);
subplot(1, 2, 1);
plot(t_fine_total, I_reconstructed, 'r-', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Current (A)');
title('Inrush Current (Reconstructed)');

subplot(1, 2, 2);
% 绘制直流偏磁和二次谐波
idx_dc = find(harmonics == 0);
idx_2nd = find(harmonics == 2);
plot(time_vec, abs(I_history_k(:, idx_dc)), 'b-', 'LineWidth', 1.5); hold on;
plot(time_vec, abs(I_history_k(:, idx_2nd)), 'g--', 'LineWidth', 1.5);
legend('DC Component', '2nd Harmonic');
grid on; xlabel('Envelope Time (s)'); ylabel('Magnitude (A)');
title('Harmonic Envelope Decay');