% TEST_FREQUENCY_COUPLED
% 测试场路耦合非线性时谐求解器 (FrequencyCoupledSolver)
%
% 场景:
%   1. 电压源驱动线圈 (Stranded, Region 2)
%   2. 下方铝板 (Massive Conductor, Region 3) 产生涡流
%   3. 验证电流相位滞后及涡流产生

clear; clc;
addpath(genpath('../src')); 

%% 1. 网格生成 (Plate + Coil)
fprintf('[Test] Generating Mesh...\n');
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
mesh.RegionTags(centers(1,:) > 0.12 & centers(1,:) < 0.18) = 3; % 3: Plate (Conductor)

mesh.generateEdges();
dofHandler = DofHandler(mesh);
space_A = FunctionSpace('Nedelec', 1);
space_V = FunctionSpace('Lagrange', 1);

%% 2. 物理参数
fprintf('[Test] Setting up Physics...\n');
% 构造非线性材料 (模拟铁芯饱和，虽然这里是线圈，但为了测试 Solver 的非线性逻辑)
B_curve = linspace(0, 2, 10);
H_curve = B_curve / (1000 * 4*pi*1e-7); % mu_r = 1000
H_curve(end) = H_curve(end) * 5; % 饱和
MatLib = containers.Map('KeyType','double','ValueType','any');
MatLib(1) = MaterialLib.createLinear(1.0); 
MatLib(2) = MaterialLib.createNonlinear(B_curve, H_curve); % Region 2 Nonlinear
MatLib(3) = MaterialLib.createLinear(1.0);

SigmaMap = containers.Map('KeyType','double','ValueType','double');
SigmaMap(3) = 3.5e7; % Plate 导电

% 线圈与电路
coilDir = [0; 1; 0];
windingObj = Winding('TestCoil', 2, 100, 0, 0.005, coilDir);

freq = 50;
circuit.R = 0.1; 
circuit.L = 0; % 纯电阻电路，感抗由 FEM 提供
circuit.V_source = 10.0 + 1j*0; % 10V, 0 deg

%% 3. 求解器配置
fprintf('[Test] Initializing Solver...\n');
assembler = Assembler(mesh, dofHandler);
solver = FrequencyCoupledSolver(assembler, freq);

% 强制边界条件 (A=0 at boundary)
dofHandler.distributeDofs(space_A);
bndNodes = find(P(1,:)<1e-6 | P(1,:)>lx-1e-6);
edges = mesh.Edges;
is_bnd = ismember(edges(1,:), bndNodes) & ismember(edges(2,:), bndNodes);
fixedDofs_A = false(dofHandler.NumGlobalDofs, 1);
fixedDofs_A(dofHandler.getGlobalIndices(space_A, find(is_bnd))) = true;

%% 4. 运行
fprintf('[Test] Running Coupled Frequency Solve...\n');
try
    [A_sol, V_sol, I_sol, info] = solver.solve(space_A, space_V, MatLib, SigmaMap, ...
                                               circuit, windingObj, ...
                                               fixedDofs_A, []);
    
    fprintf('\n---------------------------------\n');
    fprintf('           RESULTS               \n');
    fprintf('---------------------------------\n');
    fprintf('Final Current: %.4f + j%.4f A\n', real(I_sol), imag(I_sol));
    fprintf('Current Mag  : %.4f A\n', abs(I_sol));
    fprintf('Phase Angle  : %.2f deg\n', angle(I_sol)*180/pi);
    
    % 验证逻辑
    if abs(I_sol) < 1e-6
        error('Current is zero.');
    end
    
    % 简单的物理验证:
    % 感性负载，电流应滞后于电压 (Phase < 0)
    phi = angle(I_sol)*180/pi;
    if phi > 0.1
        warning('Current is leading Voltage? Check formulation signs.');
    else
        fprintf('[Check] Phase lag observed (Expected for inductive load).\n');
    end
    
    if ~isempty(V_sol)
        fprintf('[Check] Eddy currents computed in Region 3.\n');
        fprintf('        Max V_sol (phys): %.4e V\n', max(abs(V_sol)));
    else
        error('V_sol is empty (Eddy currents not computed).');
    end
    
    fprintf('[Test] PASS\n');
    
catch ME
    fprintf('[Test] FAIL: %s\n', ME.message);
    rethrow(ME);
end