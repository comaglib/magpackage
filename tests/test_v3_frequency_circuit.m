% test_v3_frequency_circuit.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Frequency Coupled Solver (RL Impedance)         \n');
fprintf('=========================================================\n');

% --- 1. 网格与模型 ---
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,1,6), linspace(0,2,11));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); 
mesh.NumNodes = size(mesh.P, 2);
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs==0
    all = dofHandler.DofMaps(space.toString()); dofHandler.NumGlobalDofs = max(all(:));
end

% 线性材料 (空气)
mu0 = 4*pi*1e-7;
NuMap = 1/mu0;
SigmaMap = 0; % 无涡流损耗，纯电感行为

% --- 2. 预计算静态电感 L (作为基准) ---
fprintf('[Step 1] Calculating Static Inductance...\n');
assembler = Assembler(mesh, dofHandler);
coil = Winding('Coil1', 1, 100, 1.0, 1.0, [0,0,1]); % N=100, R=1
C = assembler.assembleWinding(space, coil);
K = assembler.assembleStiffness(space, NuMap);

is_bnd = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);
[K_stat, F_stat] = BoundaryCondition.applyDirichlet(K, C*1.0, is_bnd);
ls = LinearSolver('Auto');
A_stat = ls.solve(K_stat, F_stat);
L_ref = C' * A_stat;
fprintf('  - Reference Inductance L: %.6e H\n', L_ref);

% --- 3. 频域耦合求解 ---
fprintf('[Step 2] Running Frequency Coupled Solver...\n');
f = 50;
omega = 2*pi*f;
R_circuit = 1.0;
V_mag = 10.0;

solver = FrequencyCoupledSolver(assembler, coil, R_circuit, f);

% 求解
[Sol, ~] = solver.solve(space, NuMap, SigmaMap, V_mag, is_bnd);

% --- 4. 验证阻抗 ---
fprintf('[Step 3] Validating Impedance...\n');
I_phasor = Sol.I;
Z_calc = V_mag / I_phasor; % V / I

% 理论阻抗 Z = R + j*w*L
Z_theo = R_circuit + 1j * omega * L_ref;

fprintf('  - Calculated Z: %.4f + %.4fi Ohm\n', real(Z_calc), imag(Z_calc));
fprintf('  - Theoretical Z: %.4f + %.4fi Ohm\n', real(Z_theo), imag(Z_theo));

err_real = abs(real(Z_calc) - real(Z_theo)) / abs(real(Z_theo));
err_imag = abs(imag(Z_calc) - imag(Z_theo)) / abs(imag(Z_theo));

fprintf('  - Real Part Error: %.2e\n', err_real);
fprintf('  - Imag Part Error: %.2e\n', err_imag);

if err_real < 1e-4 && err_imag < 1e-4
    fprintf('[PASS] Frequency coupled solver matches static inductance model.\n');
else
    fprintf('[FAIL] Impedance mismatch.\n');
end