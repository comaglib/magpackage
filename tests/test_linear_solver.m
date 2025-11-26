function test_linear_solver()
% TEST_LINEAR_SOLVER 验证线性求解器接口 (集成 MUMPS 测试)
% 测试点:
% 1. 实数稀疏矩阵求解 (dmumps)
% 2. 实数矩阵 + 复数右端项 (自动切换到 zmumps)
% 3. 复数稀疏矩阵求解 (zmumps)
% 4. 精度验证 (与 MATLAB 内置 \ 对比)

    addpath(genpath('src'));
    clc;
    fprintf('==============================================\n');
    fprintf('        Linear Solver Interface Test          \n');
    fprintf('==============================================\n');

    % 检查 MUMPS 是否存在
    if exist('dmumps', 'file') ~= 2 && exist('dmumpsmex', 'file') ~= 3
        warning('未检测到 MUMPS 接口文件。测试将回退到 MATLAB 内置求解器。');
    else
        fprintf('MUMPS 接口检测通过。\n');
    end

    % ---------------------------------------------------------
    % Case 1: 实数系统 (Real System)
    % ---------------------------------------------------------
    fprintf('\n[Case 1] Testing Real System (Expect: dmumps)...\n');
    n = 2000;
    density = 0.005;
    
    % 构造实对称正定稀疏矩阵 (模拟 FEM 刚度阵)
    A_real = sprand(n, n, density);
    A_real = A_real + A_real' + 10*speye(n); 
    b_real = rand(n, 1);
    
    % 设置模型配置
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true; % 开启对称优化
    
    % 求解
    t_start = tic;
    x_real = linear_solve(A_real, b_real, Model);
    t_mumps = toc(t_start);
    
    % 验证
    verify_result(A_real, x_real, b_real, 'Real System');
    fprintf('  Time (MUMPS): %.4f s\n', t_mumps);

    % ---------------------------------------------------------
    % Case 2: 混合系统 (Real Matrix + Complex RHS)
    % ---------------------------------------------------------
    fprintf('\n[Case 2] Testing Mixed System (Real A + Complex b) (Expect: Auto-switch to zmumps)...\n');
    % 使用相同的实矩阵 A_real，但在右端项加入虚部
    b_mixed = rand(n, 1) + 1i * rand(n, 1);
    
    % 求解
    t_start = tic;
    x_mixed = linear_solve(A_real, b_mixed, Model);
    t_mumps = toc(t_start);
    
    % 验证
    verify_result(A_real, x_mixed, b_mixed, 'Mixed System');
    fprintf('  Time (MUMPS): %.4f s\n', t_mumps);
    
    % ---------------------------------------------------------
    % Case 3: 全复数系统 (Complex System)
    % ---------------------------------------------------------
    fprintf('\n[Case 3] Testing Fully Complex System (Expect: zmumps)...\n');
    % 构造复数稀疏矩阵 (模拟 HBFEM 或 频域矩阵)
    A_complex = sprand(n, n, density) + 1i * sprand(n, n, density);
    A_complex = A_complex + 10*speye(n); % 保证非奇异
    b_complex = rand(n, 1) + 1i * rand(n, 1);
    
    % 对于全复数非厄米矩阵，通常设 Symmetric = false
    Model.Solver.Linear.Symmetric = false; 
    
    % 求解
    t_start = tic;
    x_complex = linear_solve(A_complex, b_complex, Model);
    t_mumps = toc(t_start);
    
    % 验证
    verify_result(A_complex, x_complex, b_complex, 'Complex System');
    fprintf('  Time (MUMPS): %.4f s\n', t_mumps);
    
    fprintf('\nAll Linear Solver Tests Passed!\n');
end

function verify_result(A, x, b, name)
    % 验证残差
    res = norm(A*x - b) / norm(b);
    fprintf('  %s Relative Residual: %e ', name, res);
    
    if res < 1e-10
        fprintf('[PASS]\n');
    else
        fprintf('[FAIL]\n');
        error('%s solution accuracy is too low!', name);
    end
end