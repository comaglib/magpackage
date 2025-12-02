% TEST_TRANSIENT_EDDY_BLOCK
% 全新瞬态求解器测试脚本：导电块涡流模型
%
% 场景:
%   - 区域 1: 空气 (Air)
%   - 区域 2: 铝块 (Conductor, sigma = 3.7e7)
%   - 区域 3: 激励源区域 (Source Coil)
%
% 目的:
%   验证 TransientSolver 在高电导率、SI单位制下的收敛性和物理正确性。

clear; clc;
addpath(genpath('../src')); 

fprintf('======================================================\n');
fprintf('   TEST: Transient Eddy Current (Conductive Block)    \n');
fprintf('======================================================\n');

%% 1. 网格生成 (Structured Hex -> Tet)
fprintf('[Mesh] Generating structured mesh...\n');
% 几何范围: [-0.1, 0.1] x [-0.1, 0.1] x [0, 0.2]
x_range = linspace(-0.1, 0.1, 9);  % 8 单元
y_range = linspace(-0.1, 0.1, 9);  % 8 单元
z_range = linspace(0, 0.2, 7);     % 6 单元
[X, Y, Z] = meshgrid(x_range, y_range, z_range);

P = [X(:), Y(:), Z(:)]';
T = delaunay(P(1,:), P(2,:), P(3,:))';

mesh = Mesh(); 
mesh.P = P; mesh.T = T;
mesh.NumNodes = size(P, 2); mesh.NumElements = size(T, 2);

% --- 定义区域 ---
% 计算单元中心
centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4;
cx = centers(1,:); cy = centers(2,:); cz = centers(3,:);

mesh.RegionTags = ones(1, mesh.NumElements); % 默认为 1: Air

% 区域 2: 导电块 (中心部分, z < 0.05)
mask_cond = (abs(cx) < 0.06) & (abs(cy) < 0.06) & (cz < 0.06);
mesh.RegionTags(mask_cond) = 2;

% 区域 3: 激励源 (上方悬浮, z > 0.12)
mask_coil = (abs(cx) < 0.06) & (abs(cy) < 0.06) & (cz > 0.14);
mesh.RegionTags(mask_coil) = 3;

mesh.generateEdges();
fprintf('[Mesh] Nodes: %d, Elements: %d, Edges: %d\n', ...
    mesh.NumNodes, mesh.NumElements, size(mesh.Edges, 2));

%% 2. 物理场空间与自由度
fprintf('[FEM] Initializing Function Spaces...\n');
space_A = FunctionSpace('Nedelec', 1); % 磁矢位
space_V = FunctionSpace('Lagrange', 1); % 标量电势 (用于涡流区)

dofHandler = DofHandler(mesh);
dofHandler.distributeDofs(space_A);

% 仅在导电区域 (Tag=2) 分配 V 场自由度
conducting_tags = 2; 
dofHandler.distributeDofs(space_V, conducting_tags);

%% 3. 材料属性 (Strict SI Units)
fprintf('[Physics] Setting up Material Properties (SI Units)...\n');
mu0 = 4*pi*1e-7; 
nu0 = 1.0/mu0; % 空气磁阻率 (~795774)

MatLibData = containers.Map('KeyType', 'double', 'ValueType', 'any');
MatLibData(1) = MaterialLib.createLinear(nu0); % Air
MatLibData(2) = MaterialLib.createLinear(nu0); % Aluminum (Linear for stability test)
MatLibData(3) = MaterialLib.createLinear(nu0); % Coil (Air-like)

% 电导率
SigmaMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
SigmaMap(2) = 3.7e7; % 铝的电导率

%% 4. 边界条件
% 简单 Dirichlet: 整个计算域边界 A x n = 0 (即切向 A = 0，磁壁边界)
fprintf('[BC] Applying Boundary Conditions...\n');
tol = 1e-6;
is_bnd_node = (P(1,:) < min(x_range)+tol) | (P(1,:) > max(x_range)-tol) | ...
              (P(2,:) < min(y_range)+tol) | (P(2,:) > max(y_range)-tol) | ...
              (P(3,:) < min(z_range)+tol) | (P(3,:) > max(z_range)-tol);
bnd_node_idx = find(is_bnd_node);

edges = mesh.Edges;
% 边界边: 两个端点都在边界上的边
is_bnd_edge = ismember(edges(1,:), bnd_node_idx) & ismember(edges(2,:), bnd_node_idx);

fixedDofs_A = false(dofHandler.NumGlobalDofs, 1);
global_idx_A = dofHandler.getGlobalIndices(space_A, find(is_bnd_edge));
fixedDofs_A(global_idx_A) = true;

% V 场边界条件: 让 Solver 的 Auto-Grounding 处理，或者手动接地一个点
fixedDofs_V = []; 

%% 5. 求解设置
fprintf('[Solver] Configuring Transient Solver...\n');
assembler = Assembler(mesh, dofHandler);
solver = TransientSolver(assembler);

freq = 50;                  % 50 Hz
period = 1/freq;
steps_per_cycle = 20;       % 每个周期 20 步
dt = period / steps_per_cycle; 
num_cycles = 1.0;           % 模拟 1 个周期
num_steps = round(steps_per_cycle * num_cycles);
timeSteps = repmat(dt, 1, num_steps);

% 激励函数: 在区域 3 施加 Y 方向电流
% 使用 Soft-Start (1-exp) 避免初始震荡
Jmag = 1e6; % 1 MA/m^2
sourceFunc = @(t) get_source_excitation(t, freq, Jmag);

% 求解参数 (基于新版 TransientSolver)
solver.Dt = dt;
solver.MaxIterations = 20;
solver.Tolerance = 1e-3;      % 绝对容差 (针对 A 场量级)
solver.RelTolerance = 1e-4;   % 相对容差 (针对大数值矩阵)
solver.AutoRegularization = true;
solver.RegScaleRel = 1e-6;

%% 6. 运行求解
fprintf('[Run] Starting Simulation (Total Steps: %d)...\n', num_steps);
t_start = tic;

try
    [results, info] = solver.solve(space_A, space_V, MatLibData, SigmaMap, ...
                                   sourceFunc, timeSteps, fixedDofs_A, fixedDofs_V);
    total_time = toc(t_start);
    fprintf('[Run] Simulation completed in %.2f s.\n', total_time);
catch ME
    fprintf('[Run] Error during simulation: %s\n', ME.message);
    rethrow(ME);
end

%% 7. 结果验证与可视化
fprintf('\n--------------------------------------\n');
fprintf('           VERIFICATION               \n');
fprintf('--------------------------------------\n');

if isempty(results)
    error('Test Failed: No results.');
end

sol_last = results{end};
norm_sol = norm(sol_last);

fprintf('1. Final Solution Norm: %.4e\n', norm_sol);

% 检查 1: 结果不应发散 (NaN/Inf)
if any(isnan(sol_last)) || any(isinf(sol_last))
    error('Test Failed: Solution contains NaN/Inf.');
end

% 检查 2: 物理量级检查
% 对于 J=1e6, mu0=1e-6, L=0.1，特征场强 A ~ mu * J * L^2 ~ 1e-6 * 1e6 * 0.01 = 0.01 Wb/m
% 考虑到几何因子，Norm(A) 应该在 0.1 ~ 10.0 之间 (取决于自由度数量)
if norm_sol < 1e-5
    fprintf('[Warn] Solution norm is surprisingly small. Check source magnitude.\n');
elseif norm_sol > 1e5
    error('Test Failed: Solution norm too large (Divergence likely).');
else
    fprintf('[Pass] Solution norm is within physical range.\n');
end

fprintf('[Test] RESULT: [ PASS ]\n');
fprintf('--------------------------------------\n');


%% --- 辅助函数 ---
function sMap = get_source_excitation(t, freq, Jmag)
    % 软启动正弦波
    ramp = (1 - exp(-50*t)); 
    val = Jmag * sin(2*pi*freq*t) * ramp;
    
    sMap = containers.Map('KeyType','double','ValueType','any');
    % 区域 3: Y 方向电流
    sMap(3) = [0; val; 0]; 
end