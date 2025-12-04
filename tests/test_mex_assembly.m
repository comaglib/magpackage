% TEST_MEX_ASSEMBLY 验证 C++/MEX 组装内核的正确性与性能
%
% 此脚本会：
% 1. 检查 MEX 文件是否存在。
% 2. 在单个四面体上进行精确度对比 (Correctness Check)。
% 3. 在大规模虚拟网格上进行速度对比 (Benchmark)。

clc;
fprintf('============================================================\n');
fprintf('           MagPackage MEX Assembly Test Suite               \n');
fprintf('============================================================\n');

%% 0. 环境检查
mex_name = 'assemble_curl_curl_kernel_mex';
if exist(mex_name, 'file') ~= 3
    fprintf(2, '[ERROR] MEX file "%s" not found.\n', mex_name);
    fprintf('Please run make.m first to compile the C++ kernels.\n');
    return;
else
    fprintf('[INFO] Found MEX file: %s.%s\n', mex_name, mexext);
end

% 确保路径包含 src
if ~exist('assemble_curl_curl_kernel', 'file')
    addpath(genpath('../src')); % 假设在 tests 目录下运行
end

%% 1. 正确性验证 (Correctness Test)
fprintf('\n--- [Test 1] Correctness Verification (Single Element) ---\n');

% 构造单单元数据 (标准四面体)
P = [0, 1, 0, 0; 
     0, 0, 1, 0; 
     0, 0, 0, 1]; % [3 x 4]
T = [1; 2; 3; 4]; % [4 x 1]

% 模拟数据
Dofs  = [1; 2; 3; 4; 5; 6];      % [6 x 1]
Signs = [1; -1; 1; -1; 1; -1];   % [6 x 1] 随机符号
Nu    = 500.0;                   % 测试非 1 磁阻率

% 准备积分数据
config = FemConfig.Default();
order = config.DefaultQuadratureOrder;
[q_pts, q_w] = get_quadrature_data('tet', order);
[~, curl_ref] = nedelec_tet_p1(q_pts);

% --- Run MATLAB Kernel ---
PackedData.P = P; 
PackedData.T = T;
PackedData.CellDofs = Dofs; 
PackedData.Signs = Signs;
PackedData.Nu = Nu;

[I_mat, J_mat, V_mat] = assemble_curl_curl_kernel(PackedData, config);
K_mat = sparse(I_mat, J_mat, V_mat, 6, 6);

% --- Run MEX Kernel ---
% 注意：T 在 MATLAB 中通常是 double，传入 C++ 也是 double*
[I_mex, J_mex, V_mex] = assemble_curl_curl_kernel_mex(...
    P, T, Dofs, Signs, Nu, q_w, curl_ref);
K_mex = sparse(I_mex, J_mex, V_mex, 6, 6);

% --- Compare ---
diff_norm = full(max(abs(K_mat(:) - K_mex(:))));
rel_err = diff_norm / max(abs(K_mat(:)));

fprintf('Max Absolute Error: %e\n', diff_norm);
if diff_norm < 1e-13
    fprintf('[PASS] Correctness verification successful.\n');
else
    fprintf(2, '[FAIL] Results mismatch! Error too large.\n');
end

%% 2. 性能基准测试 (Performance Benchmark)
fprintf('\n--- [Test 2] Performance Benchmark (Large Scale) ---\n');

n_rep = 100000; % 单元数量: 10万
fprintf('Generating data for %d elements...\n', n_rep);

% 简单复制数据以构造大数组
big_T = repmat(T, 1, n_rep);     % [4 x N]
big_Dofs = repmat(Dofs, 1, n_rep); % [6 x N] (虚拟DOF ID，仅测试组装速度)
big_Signs = repmat(Signs, 1, n_rep);
big_Nu = repmat(Nu, n_rep, 1);   % [N x 1]

% 更新 PackedData
PackedData.T = big_T;
PackedData.CellDofs = big_Dofs;
PackedData.Signs = big_Signs;
PackedData.Nu = big_Nu;

% 预热 (Warm-up)
assemble_curl_curl_kernel_mex(P, T, Dofs, Signs, Nu, q_w, curl_ref);

% --- Timing MATLAB ---
fprintf('Running MATLAB kernel (Parfor)... ');
tic;
[~, ~, ~] = assemble_curl_curl_kernel(PackedData, config);
t_mat = toc;
fprintf('Done. Time: %.4f s\n', t_mat);

% --- Timing MEX ---
fprintf('Running MEX kernel (OpenMP)...    ');
tic;
[~, ~, ~] = assemble_curl_curl_kernel_mex(...
    P, big_T, big_Dofs, big_Signs, big_Nu, q_w, curl_ref);
t_mex = toc;
fprintf('Done. Time: %.4f s\n', t_mex);

% --- Summary ---
speedup = t_mat / t_mex;
fprintf('\n>>> Speedup: %.2fx <<<\n', speedup);

if speedup > 5
    fprintf('[SUCCESS] MEX is significantly faster.\n');
else
    fprintf('[WARNING] Speedup is lower than expected. Check OpenMP settings.\n');
end

fprintf('\nTest Suite Completed.\n');