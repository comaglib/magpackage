% TEST_CUSTOM_TRANSIENT_NONLINEAR_VERIFIED
% 自定义非线性瞬态求解器测试 (带自动验证版)
%
% 功能:
%   1. 构造非线性材料与涡流耦合的测试案例。
%   2. 运行 TransientSolver。
%   3. [新增] 自动检查结果的有效性和范围，输出 PASS/FAIL。

clear; clc;
addpath(genpath('../src')); 

%% 1. 动态生成网格
fprintf('[Test] Generating custom mesh...\n');
lx = 0.2; ly = 0.1; lz = 0.1;
[X, Y, Z] = meshgrid(linspace(0, lx, 11), linspace(0, ly, 5), linspace(0, lz, 5));
P = [X(:), Y(:), Z(:)]'; 
T = delaunay(P(1,:), P(2,:), P(3,:))'; 

mesh = Mesh();
mesh.P = P; mesh.T = T;
mesh.NumNodes = size(P, 2); mesh.NumElements = size(T, 2);

% 定义区域: x < 0.1 (Reg 1, Iron), x >= 0.1 (Reg 2, Coil)
centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4;
mesh.RegionTags = ones(1, mesh.NumElements);
mesh.RegionTags(centers(1,:) >= 0.1) = 2;

mesh.generateEdges(); 

%% 2. 物理场
dofHandler = DofHandler(mesh);
space_A = FunctionSpace('Nedelec', 1);
space_V = FunctionSpace('Lagrange', 1);

%% 3. 定义非线性材料
fprintf('[Test] Creating Material Properties...\n');
MatLibData = containers.Map('KeyType', 'double', 'ValueType', 'any');

% 区域 1: 非线性 (饱和模型)
B_knots = [0, 0.5, 1.0, 1.4, 1.6, 1.8, 2.5, 5.0];
mu_r_init = 2000;
H_knots = [0, 0.5/(mu_r_init*4*pi*1e-7), 1.0/(mu_r_init*4*pi*1e-7), 1000, 5000, 15000, 100000, 500000];
MatLibData(1) = MaterialLib.createNonlinear(B_knots, H_knots);

% 区域 2: 线性
MatLibData(2) = MaterialLib.createLinear(1.0);

% 电导率
SigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
SigmaMap(2) = 5.0e7; 

%% 4. 边界条件
fprintf('[Test] Setting Boundary Conditions...\n');
tol = 1e-6;
is_bnd_node = (abs(P(1,:) - 0) < tol) | (abs(P(1,:) - lx) < tol) | ...
              (abs(P(2,:) - 0) < tol) | (abs(P(2,:) - ly) < tol) | ...
              (abs(P(3,:) - 0) < tol) | (abs(P(3,:) - lz) < tol);
bnd_node_indices = find(is_bnd_node);
edges = mesh.Edges;
is_bnd_edge = ismember(edges(1,:), bnd_node_indices) & ismember(edges(2,:), bnd_node_indices);

dofHandler.distributeDofs(space_A);
global_dofs_A = dofHandler.getGlobalIndices(space_A, find(is_bnd_edge));
fixedDofs_A = false(dofHandler.NumGlobalDofs, 1);
fixedDofs_A(global_dofs_A) = true;
fixedDofs_V = []; % V 场由 Solver 自动处理悬浮电位

%% 5. 激励与时间
freq = 50; Jmax = 1e6;
sourceFunc = @(t) getCustomSource(t, freq, Jmax);
dt = 1 / (freq * 20); 
timeSteps = repmat(dt, 1, 10); 

%% 6. 求解
fprintf('[Test] Running Transient Solver...\n');
assembler = Assembler(mesh, dofHandler);
solver = TransientSolver(assembler);

solver.Dt = dt;
solver.MaxIterations = 25;
solver.Tolerance = 1e-5;
solver.UseLineSearch = true; 
solver.AutoRegularization = true; % 确保开启正则化

try
    tic;
    [results, info] = solver.solve(space_A, space_V, MatLibData, SigmaMap, ...
                                   sourceFunc, timeSteps, fixedDofs_A, fixedDofs_V);
    solveTime = toc;
    fprintf('[Test] Solver completed in %.2f seconds.\n', solveTime);
catch ME
    fprintf('[Test] FAILED: Solver crashed with error: %s\n', ME.message);
    rethrow(ME);
end

%% 7. 结果验证 (Logic Check)
fprintf('\n--------------------------------------\n');
fprintf('           VERIFICATION               \n');
fprintf('--------------------------------------\n');

status = 'PASS';
fail_reason = '';

if isempty(results)
    status = 'FAIL';
    fail_reason = 'No results returned.';
else
    sol_final = results{end};
    final_norm = norm(sol_final);
    fprintf('Final Solution Norm: %.4e\n', final_norm);
    
    % Check 1: Numeric Validity
    if any(isnan(sol_final)) || any(isinf(sol_final))
        status = 'FAIL';
        fail_reason = 'Solution contains NaN or Inf.';
    end
    
    % Check 2: Range Check (Expected ~0.018 from previous valid run)
    % 允许一定的误差范围，例如 [0.001, 1.0]
    expected_range = [0.001, 1.0];
    if strcmp(status, 'PASS')
        if final_norm < expected_range(1)
            status = 'FAIL';
            fail_reason = sprintf('Result too small (Norm < %.3e). Possible trivial solution.', expected_range(1));
        elseif final_norm > expected_range(2)
            status = 'FAIL';
            fail_reason = sprintf('Result too large (Norm > %.3e). Possible divergence.', expected_range(2));
        end
    end
end

if strcmp(status, 'PASS')
    fprintf('[Test] RESULT: [ PASS ]\n');
else
    fprintf('[Test] RESULT: [ FAIL ]\n');
    fprintf('[Test] Reason: %s\n', fail_reason);
    error('Test Failed.');
end
fprintf('--------------------------------------\n');


%% --- 辅助函数 ---
function sMap = getCustomSource(t, freq, Jmax)
    val = Jmax * sin(2*pi*freq*t);
    sMap = containers.Map('KeyType','double','ValueType','any');
    sMap(2) = [0; 0; val];
end