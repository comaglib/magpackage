% TEST_TRANSIENT_COUPLED_SYSTEM_FIXED
% 测试场路耦合求解器 (Stranded Coil + Massive Plate)
%
% 修复: 
%   修正 Winding 类实例化参数不足的错误。
%   Winding(name, id, turns, resistance, area, direction)

clear; clc;
addpath(genpath('../src')); 

%% 1. 网格生成
fprintf('[Test] Generating mesh...\n');
lx=0.2; ly=0.1; lz=0.1;
[X, Y, Z] = meshgrid(linspace(0,lx,8), linspace(0,ly,6), linspace(0,lz,6));
P = [X(:), Y(:), Z(:)]';
T = delaunay(P(1,:), P(2,:), P(3,:))';
mesh = Mesh(); mesh.P=P; mesh.T=T; 
mesh.NumNodes=size(P,2); mesh.NumElements=size(T,2);

% 定义区域
centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4;
mesh.RegionTags = ones(1, mesh.NumElements); % 1: Air
mesh.RegionTags(centers(1,:) > 0.05 & centers(1,:) < 0.1) = 2; % 2: Coil
mesh.RegionTags(centers(1,:) > 0.12 & centers(1,:) < 0.18) = 3; % 3: Plate

mesh.generateEdges();
dofHandler = DofHandler(mesh);
space_A = FunctionSpace('Nedelec', 1);
space_V = FunctionSpace('Lagrange', 1);

%% 2. 材料与物理参数
fprintf('[Test] Setting up physics...\n');
MatLib = containers.Map('KeyType','double','ValueType','any');
MatLib(1) = MaterialLib.createLinear(1.0); % Air
MatLib(2) = MaterialLib.createLinear(1.0); % Coil
MatLib(3) = MaterialLib.createLinear(1.0); % Plate

SigmaMap = containers.Map('KeyType','double','ValueType','double');
SigmaMap(3) = 3.5e7; % Plate only

% 线圈定义
coilDir = [0; 1; 0]; % Y 方向
% [FIX] 补全 Winding 构造函数参数: Name, ID, Turns, Resistance, Area, Direction
% Resistance 设为 0，因为在 circuit 结构体中已经定义了总电阻
windingObj = Winding('ExcitationCoil', 2, 100, 0, 0.005, coilDir); 

% 电路参数
circuit.R = 1.0; % Ohm
circuit.L = 1e-3; % Henry
freq = 50;
circuit.V_source_func = @(t) 10 * sin(2*pi*freq*t);

%% 3. 边界条件
% 简单处理：仅固定边界棱边 A=0
fprintf('[Test] Distributing DoFs...\n');
dofHandler.distributeDofs(space_A);

% 找到边界棱边
bndNodes = find(P(1,:)<1e-6 | P(1,:)>lx-1e-6 | ...
                P(2,:)<1e-6 | P(2,:)>ly-1e-6 | ...
                P(3,:)<1e-6 | P(3,:)>lz-1e-6);
edges = mesh.Edges;
is_bnd_edge = ismember(edges(1,:), bndNodes) & ismember(edges(2,:), bndNodes);
fixedDofs_A = false(dofHandler.NumGlobalDofs, 1);
global_idx = dofHandler.getGlobalIndices(space_A, find(is_bnd_edge));
fixedDofs_A(global_idx) = true;

%% 4. 求解
fprintf('[Test] Initializing Solver...\n');
assembler = Assembler(mesh, dofHandler);
solver = TransientCoupledSolver(assembler);
solver.Dt = 1e-3;
solver.MaxIterations = 20;
solver.UseLineSearch = true;
solver.AutoRegularization = true; % 确保开启正则化以稳定空气域

timeSteps = repmat(1e-3, 1, 10); % 10 steps
init_State = [];

fprintf('[Test] Running Coupled Solver...\n');
try
    [results, info] = solver.solve(space_A, space_V, MatLib, SigmaMap, ...
                                   circuit, windingObj, ...
                                   timeSteps, fixedDofs_A, [], init_State);
                               
    % 验证结果
    if isempty(results)
        error('Solver returned empty results.');
    end
    
    sol_end = results{end};
    I_end = sol_end(end); % 电流是最后一个自由度
    fprintf('\n[Result] Final Current I = %.4f A\n', I_end);
    fprintf('[Result] Final Time = %.4f s\n', info.FinalTime);
    
    % 检查数值有效性
    if isnan(I_end) || isinf(I_end)
        fprintf('[Test] FAIL: Current is NaN/Inf.\n');
    elseif abs(I_end) < 1e-9
        fprintf('[Test] WARN: Current is near zero (Check V_source or Connection).\n');
    else
        fprintf('[Test] PASS: Current generated successfully.\n');
    end
    
catch ME
    fprintf('[Test] FAIL: Solver crashed: %s\n', ME.message);
    % rethrow(ME); % Optional: rethrow to see stack trace
end