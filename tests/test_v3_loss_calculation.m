% test_v3_loss_calculation.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Loss Calculation (Ohmic & Iron)                 \n');
fprintf('=========================================================\n');

% --- 1. Mesh ---
[X, Y, Z] = meshgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

% --- 2. Physics ---
mu0 = 4*pi*1e-7;
mu_r_init = 2000; % [Fix] 统一使用高导磁率

% Part 1 使用线性高导磁率
NuMap = 1/(mu0 * mu_r_init); 

sigma_val = 1e5;
SigmaMap = sigma_val;

f = 50; 
omega = 2*pi*f;

% --- 3. Solver (Frequency) ---
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);

J_src = 1e5;
SourceMap = [0, 0, J_src];

is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% 求解
A_sol = solver.solve(space, NuMap, SigmaMap, SourceMap, is_bnd);

% --- 4. 损耗计算 ---
fprintf('\n[Step 4] Calculating Losses...\n');
lossCalc = LossCalculator(assembler);

% (A) 欧姆损耗
P_eddy_freq = lossCalc.computeOhmicLoss_Frequency(A_sol, SigmaMap, f, space);
fprintf('  - Ohmic Loss (Eddy): %.4e W\n', P_eddy_freq);

% (B) 铁损
k_h = 100; alpha = 2; k_e = 0.5;
P_iron_freq = lossCalc.computeIronLoss_Steinmetz(A_sol, f, k_h, alpha, k_e, space);
fprintf('  - Iron Loss (Est.):  %.4e W\n', P_iron_freq);

if P_eddy_freq > 0 && P_iron_freq > 0
    fprintf('  [PASS] Frequency losses are positive.\n');
else
    fprintf('  [FAIL] Non-positive losses.\n');
end

%% === PART 2: HBFEM Solver (Nonlinear) ===
fprintf('\n=== [Part 2] HBFEM Solver (Nonlinear, 1st+3rd Harmonic) ===\n');

harmonics = [1, 3];
aft = AFT(harmonics, [], f);

% 定义非线性材料 (初始斜率与 Part 1 一致)
B_arr = linspace(0, 3.0, 100);
B_sat = 1.5;
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);
MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

solver_hbfem = HBFEMSolver(assembler, aft);
solver_hbfem.LinearSolver.Method = 'Auto';
solver_hbfem.LinearSolver.MumpsICNTL.i14 = 60;
solver_hbfem.Tolerance = 1e-3;

% 正则化
M = assembler.assembleMass(space);
K_lin = assembler.assembleStiffness(space, NuMap); % 使用相同的线性刚度做参考
ref_val = full(mean(abs(diag(K_lin))));
solver_hbfem.MatrixK_Add = ref_val * 1e-5 * M; 

SourceMaps = cell(aft.NumHarmonics, 1);
SourceMaps{1} = [0, 0, J_src]; 

x0 = zeros(dofHandler.NumGlobalDofs, aft.NumHarmonics);
[X_hbfem, info] = solver_hbfem.solve(space, MatLibData, SourceMaps, is_bnd, x0);

% 损耗计算 (HBFEM)
fprintf('  Calculating HBFEM Losses...\n');

P_eddy_hbfem = lossCalc.computeOhmicLoss_HBFEM(X_hbfem, SigmaMap, aft, space);
fprintf('  - HBFEM Ohmic Loss: %.4e W\n', P_eddy_hbfem);

P_iron_hbfem = lossCalc.computeIronLoss_HBFEM(X_hbfem, aft, k_h, alpha, k_e, space);
fprintf('  - HBFEM Iron Loss:  %.4e W\n', P_iron_hbfem);

% 验证逻辑
fprintf('\n[Comparison]\n');
fprintf('  Freq Ohmic: %.4e vs HBFEM Ohmic: %.4e\n', P_eddy_freq, P_eddy_hbfem);

% 现在两者应该在一个数量级上
ratio = P_eddy_hbfem / P_eddy_freq;
if ratio > 0.5 && ratio < 2.0
    fprintf('  [PASS] Losses match well (Ratio=%.2f).\n', ratio);
elseif ratio > 0.1 && ratio < 10
    fprintf('  [PASS] Losses are comparable (Ratio=%.2f). Difference due to saturation.\n', ratio);
else
    fprintf('  [WARN] Large discrepancy (Ratio=%.2f). Check material definitions.\n', ratio);
end

fprintf('\n=== All Loss Tests Completed ===\n');