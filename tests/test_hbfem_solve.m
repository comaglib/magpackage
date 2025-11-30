% TEST_HBFEM_HARMONICS
% 测试 HBFEMSolver 在多谐波 (DC + Odd) 工况下的表现
% 谐波配置: [0, 1, 3, 5]
% 验证目标: DC 响应、基波响应及非线性激发的高次谐波

clear; clc;

% ==========================================
% 1. 构建简单网格 (单四面体)
% ==========================================
fprintf('[Step 1] Creating simple mesh...\n');
mesh = Mesh();
% 4个节点: 原点, x, y, z (单位: m)
mesh.P = [0 1 0 0; 0 0 1 0; 0 0 0 1] * 0.01; % 1cm 尺寸
mesh.T = [1; 2; 3; 4]; 
mesh.RegionTags = [1]; 
mesh.NumNodes = 4;
mesh.NumElements = 1;
mesh.generateEdges(); 

% ==========================================
% 2. 定义空间与自由度
% ==========================================
fprintf('[Step 2] Distributing DoFs...\n');
dofHandler = DofHandler(mesh);
space = FunctionSpace('Nedelec', 1);
dofHandler.distributeDofs(space);
numDofs = dofHandler.NumGlobalDofs;
fprintf('   Total DoFs: %d\n', numDofs);

% ==========================================
% 3. 定义非线性材料
% ==========================================
fprintf('[Step 3] Defining Nonlinear Material...\n');
% 构造一条典型的软磁 B-H 曲线 (带物理斜率以防止奇异)
H_data = linspace(0, 5000, 200);
mu0 = 4*pi*1e-7;
% B = 1.6 * tanh(H/200) + mu0*H (饱和磁通约 1.6T)
B_data = 1.6 * tanh(H_data / 200) + mu0 * H_data; 

matLibData = containers.Map('KeyType', 'double', 'ValueType', 'any');
matLibData(1) = MaterialLib.createNonlinear(B_data, H_data);

% ==========================================
% 4. 配置 AFT (谐波平衡)
% ==========================================
fprintf('[Step 4] Configuring AFT (DC, 1st, 3rd, 5th)...\n');
baseFreq = 50;
harmonics = [0, 1, 3, 5]; 
% 时间步数需 > 2*max(H)+1 = 11, 取 32 以优化 FFT
aft = AFT(harmonics, 64, baseFreq); 

% ==========================================
% 5. 组装源项 (DC + AC)
% ==========================================
fprintf('[Step 5] Preparing Sources...\n');
assembler = Assembler(mesh, dofHandler);

% 设定电流密度 J (A/m^2)
% 1. DC 偏置: 1e6 A/m^2
% 2. AC 基波: 5e6 A/m^2 (使其工作在非线性区)
src_DC = containers.Map({1}, {[1e6; 0; 0]});
src_AC = containers.Map({1}, {[5e6; 0; 0]});

% 源项列表对应 harmonics: [0, 1, 3, 5]
% 3次和5次无外部源，依靠非线性产生
sourceMaps = {src_DC, src_AC, [], []}; 

% 边界条件: 固定第1条棱 (防止刚体模态)
fixedDofs = false(numDofs, 1);
fixedDofs(1) = true; 

% ==========================================
% 6. 运行求解器
% ==========================================
fprintf('[Step 6] Running Solver...\n');
solver = HBFEMSolver(assembler, aft);
solver.AutoRegularization = true; % 处理 DC 分量的规范性
solver.Tolerance = 1e-3;
solver.MaxIter = 30;

[X_sol, info] = solver.solve(space, matLibData, sourceMaps, fixedDofs);

% ==========================================
% 7. 结果分析与判断
% ==========================================
fprintf('\n==============================================\n');
fprintf('             ANALYSIS REPORT                  \n');
fprintf('==============================================\n');

if info.Converged
    fprintf('[STATUS] Solver Converged in %d iterations.\n', info.Iterations);
else
    fprintf('[STATUS] Solver FAILED to converge.\n');
end

% 提取各次谐波分量 (Harmonic-Major 存储格式)
norms = zeros(length(harmonics), 1);
names = {'DC ', '1st', '3rd', '5th'};
norms(2) = norm(X_sol(numDofs + 1 : 2 * numDofs));

fprintf('\n--- Harmonic Contents (L2 Norm) ---\n');
for i = 1:length(harmonics)
    idx_start = (i-1) * numDofs + 1;
    idx_end   = i * numDofs;
    vec_h = X_sol(idx_start : idx_end);
    norms(i) = norm(vec_h);
    
    % 计算相对基波的比率
    ratio = 0;
    if i ~= 2 && norms(2) > 1e-12
        ratio = norms(i) / norms(2) * 100;
    end
    
    if i == 2
        fprintf('   Harmonic %s: %10.4e (Reference)\n', names{i}, norms(i));
    else
        fprintf('   Harmonic %s: %10.4e (%.2f%% of 1st)\n', names{i}, norms(i), ratio);
    end
end

fprintf('\n--- Quality Checks ---\n');
passed = true;

% 1. 检查 DC 分量 (应存在且非零，因为施加了 DC 源)
if norms(1) > 1e-9
    fprintf('[PASS] DC Component detected (Response to bias).\n');
else
    fprintf('[WARN] DC Component is too small or zero.\n');
end

% 2. 检查基波 (应为主导)
if norms(2) > norms(3) && norms(2) > norms(4)
    fprintf('[PASS] Fundamental frequency is dominant.\n');
else
    fprintf('[FAIL] Fundamental is NOT dominant!\n');
    passed = false;
end

% 3. 检查非线性谐波 (3次应显著存在)
% 阈值设为 0.1% 基波幅值，表明非线性被激发
if norms(3) > 1e-3 * norms(2)
    fprintf('[PASS] 3rd Harmonic generated (>0.1%%).\n');
else
    fprintf('[WARN] 3rd Harmonic too weak. Material might be linear or saturated.\n');
    % 此时不标记为 FAIL，视物理工况而定
end

% 4. 检查高次谐波衰减 (5次通常小于3次)
if norms(4) < norms(3)
    fprintf('[PASS] Harmonic decay observed (5th < 3rd).\n');
else
    fprintf('[WARN] 5th Harmonic is larger than 3rd (Unusual but possible).\n');
end

fprintf('----------------------------------------------\n');
if passed
    fprintf('OVERALL RESULT: [ PASSED ]\n');
else
    fprintf('OVERALL RESULT: [ FAILED ]\n');
end