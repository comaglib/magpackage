% test_v3_energy.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Magnetic Energy Calculation (v3.0)              \n');
fprintf('=========================================================\n');

% --- 1. Mesh ---
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,1,6));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh(); mesh.P = DT.Points'; mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

% --- 2. Material ---
mu0 = 4*pi*1e-7;
MatLibLinear(1) = MaterialLib.createLinear(1000); % mu_r = 1000

% --- 3. Solver ---
assembler = Assembler(mesh, dofHandler);
% 伪造一个 A 解
A_sol = rand(dofHandler.NumGlobalDofs, 1);

% --- 4. 验证 ---
fprintf('[Test 1] Linear Material Energy Check...\n');
post = PostProcessor(assembler);

% Matrix Method (Reference)
K = assembler.assembleStiffness(space, 1/(1000*mu0));
W_matrix = 0.5 * A_sol' * K * A_sol;

% Integral Method (Target)
W_integral = post.computeMagneticEnergy(A_sol, MatLibLinear);

fprintf('  - Matrix Method:   %.6e J\n', W_matrix);
fprintf('  - Integral Method: %.6e J\n', W_integral);

err = abs(W_matrix - W_integral) / W_matrix;
fprintf('  - Relative Error:  %.2e\n', err);

if err < 1e-3
    fprintf('[PASS] Energy calculation matches for linear case.\n');
else
    fprintf('[FAIL] Energy mismatch.\n');
end

% Non-linear Check
fprintf('\n[Test 2] Nonlinear Material Energy Check...\n');
B_arr = linspace(0, 3.0, 50);
H_arr = B_arr ./ (1000 * mu0); % Linear data wrapped in nonlinear struct
MatLibNonlin(1) = MaterialLib.createNonlinear(B_arr, H_arr);

W_nonlin = post.computeMagneticEnergy(A_sol, MatLibNonlin);
fprintf('  - Nonlinear Integral: %.6e J\n', W_nonlin);

err_nl = abs(W_nonlin - W_matrix) / W_matrix;
if err_nl < 1e-2
    fprintf('[PASS] Nonlinear integration matches linear value.\n');
end