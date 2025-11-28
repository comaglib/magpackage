% test_v3_step3.m
clear; clc;
addpath(genpath('src'));

fprintf('=== Step 3: Source Integration (LinearForm) ===\n');

% 1. 构造一个长方体网格 [0,1]x[0,1]x[0,1]
% 为了积分验证，我们把它切成 6 个四面体
fprintf('[Setup] Creating Unit Cube Mesh...\n');
mesh = Mesh();
mesh.P = [0 1 0 1 0 1 0 1; 
          0 0 1 1 0 0 1 1; 
          0 0 0 0 1 1 1 1]; % 8个顶点
% 手动 Delaunay 或者调用 mesh generator (如果有)
% 这里简单起见，我们定义一个简单的 T
dt = delaunayTriangulation(mesh.P');
mesh.T = dt.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); % 都是区域 1
mesh.NumNodes = 8;
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

% 2. 准备数据
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
packedData = dofHandler.packForKernel(space);
packedData.P = mesh.P;
packedData.T = mesh.T;
packedData.RegionTags = mesh.RegionTags;

% 3. 定义源项
% 在区域 1 施加 J = [0, 0, 1] A/m^2
% 总体积 = 1 m^3
% 总电流矩 (Current Moment) Integral(J) dV = [0, 0, 1]
J_mag = 1.0;
src = CurrentSource(1, [0, 0, J_mag]);

% 预处理 SourceMap
% 假设最大 Region ID 是 1
SourceMap = zeros(1, 3);
SourceMap(src.RegionID, :) = src.J_Vector;

% 4. 组装 F
conf = FemConfig.Default();
F = assemble_source_kernel(packedData, SourceMap, conf);

fprintf('  - Assembled F vector with %d non-zeros.\n', nnz(F));

% 5. 物理验证
% FEM 方程: K * A = F
% 右端项 F 的物理意义:
% F_i = Integral( N_i . J ) dV
% 让我们检验 F 的"总和"是否某种程度上反映了电流方向？
% 不容易直接检验 F 的和。但我们可以检验做功。
% 如果 A 是匀强场 A = [-y/2, x/2, 0] (对应的 B=[0,0,1])? 不对。
% 让我们用一个更简单的检验：
% 如果 A 是常数向量 A_test = [0, 0, 1] (非物理，但用于数学检验)
% Integral( A_test . J ) dV = Integral( 1 * Jz ) dV = Vol * Jz = 1.0 * 1.0 = 1.0
% 离散形式: A_vec' * F 应该等于 1.0
% 这里的 A_vec 是将 A_test 投影到 Edge DoFs 上的系数

fprintf('[Verification] Checking Integral( A_test . J ) dV ...\n');

% 构造测试场 A_test = [0, 0, 1]
% Edge DoF_i = Integral( A_test . t_i ) dl (对于 P1 Nedelec, 近似为 A_test . vec_edge)
% 实际上就是 A_test 点乘 棱向量
A_test = [0; 0; 1];
edges = mesh.Edges;
numEdges = size(edges, 2);
A_vec_test = zeros(numEdges, 1);

for e = 1:numEdges
    n1 = edges(1, e);
    n2 = edges(2, e);
    p1 = mesh.P(:, n1);
    p2 = mesh.P(:, n2);
    edge_vec = p2 - p1;
    
    A_vec_test(e) = dot(A_test, edge_vec);
end

% 计算做功: W = A_vec_test' * F
Work = A_vec_test' * F;

fprintf('  - Calculated Work (A_vec . F): %.4f\n', Work);
fprintf('  - Expected Work (J . A * Vol): %.4f\n', 1.0);

if abs(Work - 1.0) < 1e-4 % 放宽一点公差
    fprintf('  [PASS] Source Integral is physically consistent.\n');
else
    fprintf('  [FAIL] Source Integral mismatch!\n');
end

% 6. 测试反向
SourceMap_Neg = SourceMap * -1;
F_neg = assemble_source_kernel(packedData, SourceMap_Neg, conf);
if norm(F + F_neg) < 1e-10
    fprintf('  [PASS] Linearity verified (F_neg = -F).\n');
else
    fprintf('  [FAIL] Linearity check failed.\n');
end