% test_hbfem_solve.m
% 测试纯 HBFEM 求解器 (非耦合)
% 谐波配置: DC + 基波 + 3次 + 5次
% 修复: 修正 info 结构体字段访问错误 (FinalResidual -> Residuals(end))

clear; clc;
addpath(genpath('src'));

fprintf('=========================================================\n');
fprintf('   Test: HBFEM Field Solver (Harmonics: 0, 1, 3, 5)\n');
fprintf('=========================================================\n');

%% 1. 网格生成 (Mesh Generation)
% 生成一个 1m x 0.2m x 0.2m 的简易长方体
[X, Y, Z] = meshgrid(linspace(0,1,6), linspace(0,0.2,3), linspace(0,0.2,3));
DT = delaunayTriangulation(X(:), Y(:), Z(:));
mesh = Mesh(); 
mesh.P = DT.Points'; 
mesh.T = DT.ConnectivityList';
mesh.RegionTags = ones(1, size(mesh.T, 2)); % 单一区域
mesh.NumNodes = size(mesh.P, 2); 
mesh.NumElements = size(mesh.T, 2);
mesh.generateEdges();

fprintf('[Mesh] Nodes: %d, Elements: %d, Edges: %d\n', mesh.NumNodes, mesh.NumElements, size(mesh.Edges, 2));

%% 2. 物理空间 (Function Space)
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
numDofs = dofHandler.NumGlobalDofs;

%% 3. AFT 与 谐波设置 (Harmonics)
% 要求配置: [0, 1, 3, 5]
harmonics = [0, 1, 3, 5];
base_freq = 50; 
% 自动计算所需时间步数 (Nyquist: > 2*max(H) + 1 => 11)
% 为了 FFT 效率取 16 或 32
num_time_steps = 32; 

% 注意: 根据 AFT.m 定义，参数顺序为 (harmonics, n_time_steps, base_freq)
aft = AFT(harmonics, num_time_steps, base_freq);
fprintf('[AFT] Freq: %g Hz, Steps: %d, Harmonics: %s\n', base_freq, num_time_steps, mat2str(harmonics));

%% 4. 材料定义 (Manual Spline Material)
% 构造非线性材料数据
mat_iron = struct();
mat_iron.Type = 'Nonlinear';
mat_iron.Name = 'TestIron';
mat_iron.nu0 = 1 / (1000 * 4*pi*1e-7); % Initial mu_r = 1000
mat_iron.MaxBSq = 20.0; 

% 生成 B^2 - Nu 曲线
b_sq_pts = linspace(0, 5, 20); 
mu_r_eff = 1000 ./ (1 + 2.0 * b_sq_pts) + 10; % 饱和模型
nu_pts = (1 ./ (mu_r_eff * 4*pi*1e-7))';

mat_iron.SplineNu = spline(b_sq_pts, nu_pts);
% 维度修正
if size(mat_iron.SplineNu.coefs, 1) ~= (length(mat_iron.SplineNu.breaks)-1)
    mat_iron.SplineNu.coefs = mat_iron.SplineNu.coefs'; 
end

MatLibData = containers.Map('KeyType', 'double', 'ValueType', 'any');
MatLibData(1) = mat_iron;

%% 5. 组装源项 (Source Assembly)
assembler = Assembler(mesh, dofHandler);

% 定义空间上的电流密度 J = 1e6 A/m^2 (沿 Z 轴)
J_vec = [0; 0; 1e6];
SourceMap = containers.Map({1}, {J_vec});

% 计算空间源向量 F_spatial (NumDofs x 1)
F_spatial = assembler.assembleSource(space, SourceMap);

% 构造 HBFEM 右端项矩阵 F_total (NumDofs x NumHarmonics)
F_total = zeros(numDofs, aft.NumHarmonics);

% 将源施加在基波分量上 (Harmonic index 2 corresponds to h=1)
% Harmonics array: [0, 1, 3, 5] -> Indices: 1, 2, 3, 4
idx_fund = find(harmonics == 1); 
F_total(:, idx_fund) = F_spatial;

fprintf('[Source] Applied Jz=1e6 on Harmonic k=%d (Order 1)\n', idx_fund);

%% 6. 求解器配置与运行 (Solver)
solver = HBFEMSolver(assembler, aft);
solver.Tolerance = 1e-4;
solver.MaxIter = 20;
solver.LinearSolver.Method = 'MUMPS'; 

% 边界条件: 外表面 A=0
fixedDofs = BoundaryCondition.findOuterBoundaryDofs(mesh, dofHandler, space);

fprintf('[Solver] Starting Newton Iterations...\n');
x_init = zeros(numDofs, aft.NumHarmonics);

try
    % 调用求解
    % 接口: (space, matLib, RHS_Matrix, fixedDofs, x_init)
    [x_sol, info] = solver.solve(space, MatLibData, F_total, fixedDofs, x_init);
    
    fprintf('   > Converged: %d\n', info.Converged);
    fprintf('   > Iterations: %d\n', info.Iterations);
    
    % [修复点] 使用 Residuals(end) 获取最终残差
    final_res = info.Residuals(end);
    fprintf('   > Final Residual: %.2e\n', final_res);
    
    %% 7. 结果检查
    if info.Converged
        % 提取基波分量
        A_fund = x_sol(:, idx_fund);
        max_A = max(abs(A_fund));
        fprintf('[Result] Max Fundamental A: %.4e\n', max_A);
        
        % 提取 3 次谐波分量 (检查是否产生非线性谐波)
        idx_3rd = find(harmonics == 3);
        A_3rd = x_sol(:, idx_3rd);
        max_A3 = max(abs(A_3rd));
        fprintf('[Result] Max 3rd Harmonic A: %.4e (Ratio to Fund: %.2f%%)\n', ...
            max_A3, (max_A3/max_A)*100);
        
        % 简单的可视化
        viz = Visualizer(PostProcessor(assembler));
        figure(1); clf;
        viz.plotFieldMagnitude(A_fund, space, 'B');
        title('Fundamental B Field');
        colorbar;
        drawnow;
    end
    
catch ME
    fprintf('[ERROR] %s\n', ME.message);
    % disp(ME.stack(1));
end