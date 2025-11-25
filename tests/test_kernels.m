% tests/test_kernels.m
% 验证生成的有限元内核的基本物理性质

addpath(genpath('src'));

% 1. 定义一个标准单位四面体
P = [0 1 0 0; 
     0 0 1 0; 
     0 0 0 1]; % N1=(0,0,0), N2=(1,0,0)...
nu = 1000;     % 磁阻率
sigma = 5.8e7; % 电导率

% 2. 计算矩阵
K = Ke_curl_curl(P, nu);
M = Me_edge_edge(P, sigma);

% 3. 验证性质
fprintf('Testing Kernel Properties:\n');

% 3.1 零特征值测试 (常数场应产生零能量)
% 棱单元的零空间包含梯度假定。如果 A = grad(V)，则 Curl(A) = 0，故 K*A 应为 0。
% 在单位四面体上，如果 V=x，则 grad(V) = (1,0,0)。
% 需将 grad(V) 投影到棱上: A_i = (grad V) dot (edge_vector)
edge_vecs = P(:,[2 3 4 3 4 4]) - P(:,[1 1 1 2 2 3]); % 6条棱向量
% Edge list: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
% Edge 1 (1->2): (1,0,0)
grad_V = [1; 0; 0];
A_test = edge_vecs' * grad_V; % 投影

energy = A_test' * K * A_test;
if abs(energy) < 1e-10
    fprintf('[PASS] Curl-Curl Matrix Null-Space Check (Energy = %g)\n', energy);
else
    fprintf('[FAIL] Curl-Curl Matrix Null-Space Check (Energy = %g)\n', energy);
end

% 3.2 对称性测试
if norm(K - K', 'fro') < 1e-10
    fprintf('[PASS] Stiffness Matrix Symmetry\n');
else
    fprintf('[FAIL] Stiffness Matrix Symmetry\n');
end

if norm(M - M', 'fro') < 1e-10
    fprintf('[PASS] Mass Matrix Symmetry\n');
else
    fprintf('[FAIL] Mass Matrix Symmetry\n');
end

% 3.3 耦合矩阵一致性
% div(J) = 0 -> 强耦合.