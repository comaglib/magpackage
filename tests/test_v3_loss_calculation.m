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
NuMap = [1/mu0];

sigma_val = 1e5;
SigmaMap = [sigma_val];

f = 50; 
omega = 2*pi*f;

% --- 3. Solver ---
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

% (A) 欧姆损耗 (Eddy Current)
P_eddy = lossCalc.computeOhmicLoss_Frequency(A_sol, SigmaMap, f, space);
fprintf('  - Ohmic Loss (Eddy): %.4e W\n', P_eddy);

% (B) 铁损 (Steinmetz)
k_h = 100; alpha = 2; k_e = 0.5;
P_iron = lossCalc.computeIronLoss_Steinmetz(A_sol, f, k_h, alpha, k_e, space);
fprintf('  - Iron Loss (Est.):  %.4e W\n', P_iron);

% --- 5. 验证 ---
% 检查数值是否为正值且合理
if P_eddy > 0
    fprintf('[PASS] Ohmic loss is positive.\n');
else
    fprintf('[FAIL] Ohmic loss is non-positive!\n');
end

if P_iron > 0
    fprintf('[PASS] Iron loss is positive.\n');
else
    fprintf('[FAIL] Iron loss is non-positive!\n');
end

% 能量平衡检查 (定性)
% 欧姆损耗应与 J_src 做功有关
% P_in = 0.5 * real( Integral( J_src * E* ) ) 
% E = -j*w*A
% P_in = 0.5 * real( Integral( J_src * (j*w*A') ) )
% 注意: 这里忽略了边界项
F_vec = assembler.assembleSource(space, SourceMap);
P_input = 0.5 * real( 1j * omega * (A_sol' * F_vec) ); 
% 注意：这是电源注入到磁场的功率。
% 在涡流问题中，电源功率 = 欧姆损耗 + 磁场储能变化(虚部)
% 但由于 A 和 J 的相位关系，这项计算比较复杂。
% 我们主要验证 P_eddy 和 P_iron 的计算过程本身是否不再报错。

fprintf('  - (Ref) Input Power calc check: %.4e W (Should match Ohmic if pure loss)\n', abs(P_input));