% test_v3_step5.m
clear; clc;
addpath(genpath('src'));

fprintf('=== Step 5: Nonlinear Jacobian Test ===\n');

% 1. 构造单单元测试案例
mesh = Mesh();
mesh.P = [0 1 0 0; 0 0 1 0; 0 0 0 1]; % 标准四面体
mesh.T = [1; 2; 3; 4];
mesh.RegionTags = 1;
mesh.generateEdges();

dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

% 2. 定义非线性材料
% 构造一个简单的 nu = 1 + B^2 的材料
% nu' = 1.0
% B = 0 -> 10 的曲线
B_curve = linspace(0, 10, 100);
nu_vals = 1000 + 100 * B_curve.^2; % 随 B 增加而增加 (非物理，但便于测试导数)
% 注意: 物理材料 nu 通常随 B 增加而增加 (饱和)，或者 B-H 曲线斜率下降
H_curve = nu_vals .* B_curve;

MatLibData(1) = MaterialLib.createNonlinear(B_curve, H_curve);

% 3. 准备数据
packedData = dofHandler.packForKernel(space);
packedData.P = mesh.P;
packedData.T = mesh.T;
packedData.RegionTags = mesh.RegionTags;

% 4. 随机产生一个解 A0
numDofs = 6;
rng(42);
A0 = rand(numDofs, 1);

% 5. 计算解析 Jacobian 和 Residual
fprintf('[Test] Calculating Analytical Jacobian...\n');
config = FemConfig.Default();
[I, J, V, R0] = assemble_jacobian_kernel(packedData, A0, MatLibData, config);
K_tan_ana = full(sparse(I, J, V, numDofs, numDofs));

% 6. 有限差分验证
fprintf('[Test] Calculating Finite Difference Jacobian...\n');
epsilon = 1e-7;
K_tan_fd = zeros(numDofs, numDofs);

for i = 1:numDofs
    A_pert = A0;
    A_pert(i) = A_pert(i) + epsilon;
    
    [~, ~, ~, R_pert] = assemble_jacobian_kernel(packedData, A_pert, MatLibData, config);
    
    % J_col = (R(A+eps) - R(A)) / eps
    col = (R_pert - R0) / epsilon;
    K_tan_fd(:, i) = full(col);
end

% 7. 对比
diff_norm = norm(K_tan_ana - K_tan_fd, 'fro');
rel_err = diff_norm / norm(K_tan_ana, 'fro');

fprintf('  - Analytical vs FD Difference Norm: %.4e\n', diff_norm);
fprintf('  - Relative Error: %.4e\n', rel_err);

if rel_err < 1e-4 % 由于积分误差，FD 不会完全等于解析解
    fprintf('  [PASS] Jacobian is consistent.\n');
else
    fprintf('  [FAIL] Jacobian mismatch!\n');
    disp('Analytical:'); disp(K_tan_ana);
    disp('FD:'); disp(K_tan_fd);
end