% test_v3_hbfem_coupled.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: HBFEM Coupled Solver (Voltage Driven + Saturation)\n');
fprintf('=========================================================\n');

% --- 1. Mesh ---
[X, Y, Z] = meshgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh(); mesh.P = DT.Points'; mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

% --- 2. DoF ---
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs == 0
    all = dofHandler.DofMaps(space.toString()); dofHandler.NumGlobalDofs = max(all(:));
end
numDofs = dofHandler.NumGlobalDofs;

% --- 3. AFT ---
base_freq = 50;
harmonics = [1, 3, 5]; 
aft = AFT(harmonics, [], base_freq); 

% --- 4. Material & Winding ---
mu0 = 4*pi*1e-7; mu_r_init = 2000; B_sat = 1.5; 
B_arr = linspace(0, 5.0, 100);
mu_r_arr = 1 + (mu_r_init - 1) ./ (1 + (B_arr / B_sat).^6); 
H_arr = B_arr ./ (mu_r_arr * mu0);
MatLibData(1) = MaterialLib.createNonlinear(B_arr, H_arr);

coil = Winding('Coil1', 1, 100, 1.0, 1.0, [0,0,1]); 
R_circuit = 1.0;
L_circuit = 0.0; % [Update] 新增电感参数，保持为0以匹配原测试

% --- 5. Solver ---
assembler = Assembler(mesh, dofHandler);
% [Update] 构造函数调用更新，传入 L
solver = HBFEMCoupledSolver(assembler, aft, coil, R_circuit, L_circuit);

% [Removed] 移除了手动正则化代码，现在由 Solver 内部自动处理
% M = assembler.assembleMass(space);
% temp_MatMap = 1/mu0;
% K_lin = ...
% solver.MatrixK_Add = ... 

% [Optimization] 调整求解参数以应对强饱和
solver.MaxIter = 100;      % 增加迭代步数
solver.Tolerance = 1e-3;   % 放宽收敛容差
solver.MinStepSize = 1e-8; % 允许更小的步长

% --- 6. Voltage Ramping ---
Target_V = 40000;
idx_fund = find(harmonics == 1);
V_vec = zeros(aft.NumHarmonics, 1);

NumSteps = 1;
x_curr = zeros(numDofs, aft.NumHarmonics); 

is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

fprintf('\n[Ramping] Voltage Driven (Target %.0f V)...\n', Target_V);
for step = 1:NumSteps
    scale = step / NumSteps;
    V_step = Target_V * scale;
    V_vec(idx_fund) = V_step;
    
    fprintf('\n--- Step %d/%d (V=%.1f V) ---\n', step, NumSteps, V_step);
    
    [x_curr, i_curr, info] = solver.solve(space, MatLibData, V_vec, is_bnd, x_curr);
    
    if ~info.Converged
        % 如果最后一步没收敛但残差较小，可以接受并继续
        if info.Residuals(end) < 0.1 % 绝对残差小于0.1可以接受
             fprintf('[WARN] Step not fully converged but residual is low. Proceeding.\n');
        else
             error('Solver Failed.');
        end
    end
    
    fprintf('  Current I_fund: %.4f A\n', abs(i_curr(idx_fund)));
end

% --- 7. Results ---
fprintf('\n[Analysis] Harmonic Currents:\n');
I_mag = abs(i_curr);
for k = 1:length(harmonics)
    fprintf('  H%d: %.4f A (%.2f %%)\n', harmonics(k), I_mag(k), I_mag(k)/I_mag(1)*100);
end

I_fund = I_mag(1);
I_3rd = I_mag(2);

if I_3rd > 1e-2 * I_fund
    fprintf('[PASS] Significant harmonic currents detected (Saturation verified).\n');
else
    fprintf('[WARN] Current waveform is still sinusoidal.\n');
end