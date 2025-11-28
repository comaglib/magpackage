% test_v3_step1.m
clear; clc;
addpath(genpath('src'));

% 1. 加载网格 (OOP 方式)
% 假设当前目录下有 example_iron_coil.msh，或者你可以造一个简单的 mesh
mesh = Mesh.load('tests/example_iron_coil.msh'); 
mesh.stats();

% 2. 定义空间
V_h = FunctionSpace('Lagrange', 1); % 标量电位
Q_h = FunctionSpace('Nedelec', 1);  % 磁矢位

% 3. 初始化 DofHandler
dofHandler = DofHandler(mesh);

% 4. 分发自由度
nodesMap = dofHandler.distributeDofs(V_h);
edgesMap = dofHandler.distributeDofs(Q_h);

% 5. 验证 Nedelec 映射是否触发了 Edge 生成
mesh.stats(); % 现在应该显示 Edges 数量了

% 6. 检查数据剥离 (Data Stripping)
packed_Q = dofHandler.packForKernel(Q_h);

fprintf('\n[Success] Step 1 Refactoring Complete.\n');
fprintf('Packed Nedelec Data contains Signs: %d\n', isfield(packed_Q, 'Signs'));