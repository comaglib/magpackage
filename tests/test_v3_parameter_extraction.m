% test_v3_parameter_extraction.m
clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: Parameter Extraction (Inductance Matrix)        \n');
fprintf('=========================================================\n');

% --- 1. 网格生成 ---
% 生成一个包含两个线圈区域的网格
% 区域 1: 线圈 1 (左)
% 区域 2: 线圈 2 (右)
% 区域 3: 空气背景 (隐式包含，或者是 Mesh 默认区域)
% 这里为了简化，我们在一个大盒子里定义两个分离的区域 Tag
fprintf('[Step 1] Generating Mesh...\n');

% 我们用 create_test_mesh_direct 或者手动生成简单网格
% 这里手动生成一个 3x1x1 的长方体，中间切开
[X, Y, Z] = meshgrid(linspace(0,3,13), linspace(0,1,5), linspace(0,1,5));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh();
mesh.P = DT.Points';
mesh.T = DT.ConnectivityList';
mesh.NumNodes = size(mesh.P, 2); mesh.NumElements = size(mesh.T, 2);

% 定义区域: X < 1.0 -> Region 1, X > 2.0 -> Region 2, 中间 -> Region 3
centroids = (mesh.P(:, mesh.T(1,:)) + mesh.P(:, mesh.T(2,:)) + ...
             mesh.P(:, mesh.T(3,:)) + mesh.P(:, mesh.T(4,:))) / 4;
cx = centroids(1, :);

tags = ones(1, mesh.NumElements) * 3; % Default Air
tags(cx < 1.0) = 1; % Coil 1
tags(cx > 2.0) = 2; % Coil 2
mesh.RegionTags = tags;

mesh.generateEdges();

% --- 2. 自由度 ---
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
if isprop(dofHandler, 'NumGlobalDofs') && dofHandler.NumGlobalDofs==0
    all = dofHandler.DofMaps(space.toString()); dofHandler.NumGlobalDofs = max(all(:));
end

% --- 3. 定义绕组 ---
fprintf('[Step 2] Defining Windings...\n');
% Coil 1: N=100, R=1, Region=1, Dir=Z
w1 = Winding('Primary', 1, 100, 1.0, 1.0*1.0, [0,0,1]);
% Coil 2: N=200, R=2, Region=2, Dir=Z
w2 = Winding('Secondary', 2, 200, 2.0, 1.0*1.0, [0,0,1]);

windings = {w1, w2};

% --- 4. 提取参数 ---
fprintf('[Step 3] Running Parameter Extractor...\n');
assembler = Assembler(mesh, dofHandler);
extractor = ParameterExtractor(assembler);

mu0 = 4*pi*1e-7;
% NuMap 需要覆盖所有区域 ID (1, 2, 3)
% 都是空气
NuMap = [1/mu0; 1/mu0; 1/mu0]; 

[L_mat, R_mat] = extractor.extractMatrices(space, NuMap, windings);

% --- 5. 结果验证 ---
fprintf('\n[Results]\n');
disp('Inductance Matrix L (H):');
disp(L_mat);

disp('Resistance Matrix R (Ohm):');
disp(R_mat);

% 验证互易性 (Reciprocity) L12 == L21
L12 = L_mat(1, 2);
L21 = L_mat(2, 1);
diff_rel = abs(L12 - L21) / max(abs(L12), 1e-9);

fprintf('Reciprocity Error |L12 - L21|/L12: %.2e\n', diff_rel);

% 验证自感正定性
if L_mat(1,1) > 0 && L_mat(2,2) > 0
    fprintf('[PASS] Self inductances are positive.\n');
else
    fprintf('[FAIL] Self inductances invalid.\n');
end

% 验证耦合系数 k = M / sqrt(L1 L2) < 1
k_coupling = L12 / sqrt(L_mat(1,1) * L_mat(2,2));
fprintf('Coupling Coefficient k: %.4f\n', k_coupling);

if k_coupling >= 0 && k_coupling < 1.0
    fprintf('[PASS] Coupling coefficient is physically valid (0 <= k < 1).\n');
else
    fprintf('[WARN] Coupling coefficient out of range (maybe mesh too coarse or tight coupling).\n');
end

if diff_rel < 1e-5
    fprintf('[PASS] Reciprocity verified.\n');
else
    fprintf('[FAIL] Reciprocity check failed.\n');
end