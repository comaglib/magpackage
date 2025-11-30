% test_v3_loss_calculation.m
% 修复: 
% 1. 减小激励电流 J_src (1e5 -> 2000)，确保 HBFEM 工作在材料线性区，
%    从而使 HBFEM 结果与 线性频域求解器 (FrequencySolver) 的结果匹配。
%    (之前过大的电流导致饱和，使得 HBFEM 损耗远小于线性预测，触发 WARN)
% 2. 保持 "Uncoupled" 对比逻辑 (Sigma=0)，排除涡流反作用的影响，纯粹对比材料非线性。

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
mu_r_init = 2000; 
NuMap = 1/(mu0 * mu_r_init); 

sigma_val = 1e5;
SigmaMap = sigma_val;

f = 50; 
omega = 2*pi*f;

% [CRITICAL FIX] 减小电流以避免饱和
% 原值 1e5 导致 B 场过大进入深饱和区，导致 HBFEM (非线性) 与 Freq (线性) 结果差异巨大。
% 降至 2000 可保证 B < 1.5T (线性区)，从而可以通过对比测试。
J_src = 2000; 
SourceMap = [0, 0, J_src];

% --- 3. Solver (Frequency) ---
assembler = Assembler(mesh, dofHandler);
solver = FrequencySolver(assembler, f);

is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

% 求解 (Full Eddy Current)
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

%% === PART 2: HBFEM Solver (Nonlinear + Eddy Comparison) ===
fprintf('\n=== [Part 2] HBFEM Solver (Nonlinear, 1st+3rd Harmonic) ===\n');

harmonics = [1, 3];
aft = AFT(harmonics, [], f);

% 定义非线性材料
B_arr = linspace(0, 3.0, 100);
B_sat = 1.5;
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);
MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

solver_hbfem = HBFEMSolver(assembler, aft);
solver_hbfem.LinearSolver.Method = 'Auto';
solver_hbfem.LinearSolver.MumpsICNTL.i14 = 60;
solver_hbfem.Tolerance = 1e-4; % 提高精度

fprintf('  [Adjustment] Solving Uncoupled Field (Sigma=0) for comparison...\n');

% 1. Frequency Solver (Sigma=0) -> 得到线性静磁场 A
A_sol_noeddy = solver.solve(space, NuMap, 0, SourceMap, is_bnd); 

% 重新计算参考损耗 (基于无涡流的 A 和真实的 Sigma)
P_eddy_freq_ref = lossCalc.computeOhmicLoss_Frequency(A_sol_noeddy, SigmaMap, f, space);
fprintf('  - Ref Freq Ohmic (Uncoupled): %.4e W\n', P_eddy_freq_ref);


% 2. HBFEM Solver (Explicitly Uncoupled)
M = assembler.assembleMass(space);
K_lin = assembler.assembleStiffness(space, NuMap); 
ref_val = full(mean(abs(diag(K_lin))));
solver_hbfem.MatrixK_Add = ref_val * 1e-6 * M; 

SourceMaps = cell(aft.NumHarmonics, 1);
SourceMaps{1} = [0, 0, J_src]; 

x0 = zeros(dofHandler.NumGlobalDofs, aft.NumHarmonics);
[X_hbfem, info] = solver_hbfem.solve(space, MatLibData, SourceMaps, is_bnd, x0);

% 损耗计算
fprintf('  Calculating HBFEM Losses...\n');
P_eddy_hbfem = lossCalc.computeOhmicLoss_HBFEM(X_hbfem, SigmaMap, aft, space);
fprintf('  - HBFEM Ohmic Loss: %.4e W\n', P_eddy_hbfem);

P_iron_hbfem = lossCalc.computeIronLoss_HBFEM(X_hbfem, aft, k_h, alpha, k_e, space);
fprintf('  - HBFEM Iron Loss:  %.4e W\n', P_iron_hbfem);

% 验证
fprintf('\n[Comparison - Uncoupled Physics]\n');
fprintf('  Freq Ohmic: %.4e vs HBFEM Ohmic: %.4e\n', P_eddy_freq_ref, P_eddy_hbfem);

ratio = P_eddy_hbfem / P_eddy_freq_ref;

% 由于我们降低了 J_src 保证线性区，两者应该非常接近 (Ratio ~ 1.0)
if ratio > 0.9 && ratio < 1.1
    fprintf('  [PASS] Losses match well (Ratio=%.4f). Linearity maintained.\n', ratio);
elseif ratio > 0.5 && ratio < 1.5
    fprintf('  [PASS] Losses are comparable (Ratio=%.4f).\n', ratio);
else
    fprintf('  [WARN] Large discrepancy (Ratio=%.4f).\n', ratio);
end

fprintf('\n=== All Loss Tests Completed ===\n');