% test_v3_step2.m
clear; clc;
addpath(genpath('src'));

fprintf('=== Step 2 (Fix & Opt): Testing Kernels with Config ===\n');

% 1. 验证配置模块
conf = FemConfig.Default();
conf.AssemblyChunkSize = 20000; % 测试自定义配置
fprintf('[Config] ChunkSize set to %d\n', conf.AssemblyChunkSize);

% 2. (跳过基础单元测试，直接测试组装)
% ... (之前的 Lagrange/Jacobian 测试可以保留) ...

fprintf('[Test 3] Testing Kernel Assembly (2 Elements)...\n');

% 构造测试网格
mesh = Mesh();
% 两个共面四面体 (金字塔切分)
mesh.P = [0 1 1 0 0; 0 0 1 1 0; 0 0 0 0 1]; 
mesh.T = [1 2; 2 3; 4 4; 5 5]; 
mesh.RegionTags = [1, 1];
mesh.NumNodes = 5;
mesh.NumElements = 2;
mesh.generateEdges();

% 准备数据
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);

packedData = dofHandler.packForKernel(space);
packedData.P = mesh.P;
packedData.T = mesh.T;
packedData.Nu = [1.0; 1.0];

% 显式检查 Signs 类型
if isa(packedData.Signs, 'int8')
    fprintf('  [Check] Signs is int8 as expected. Verifying kernel fix...\n');
end

% 调用内核 (传入 Config)
try
    [I, J, V] = assemble_curl_curl_kernel(packedData, conf);
    K = sparse(I, J, V);
    
    fprintf('  - Assembled Matrix Size: %d x %d\n', size(K,1), size(K,2));
    fprintf('  - Non-zeros: %d\n', nnz(K));
    
    % 检查对称性
    if norm(K - K', 'fro') < 1e-10
        fprintf('  [PASS] Stiffness Matrix is Symmetric.\n');
    else
        fprintf('  [FAIL] Matrix is NOT Symmetric! Err=%.2e\n', norm(K - K', 'fro'));
    end
    
catch ME
    fprintf('  [FAIL] Kernel execution failed: %s\n', ME.message);
    rethrow(ME);
end

% ... (Grad Field Test 保持不变) ...