function test_hbfem_voltage()
% TEST_HBFEM_VOLTAGE 电压驱动 HBFEM 验证脚本 (注释增强版)
%
% =========================================================================
% 测试目的 (Objective):
%   验证谐波平衡有限元求解器 (HBFEM) 在电压驱动模式下处理强非线性场路耦合问题的能力。
%
% 物理场景 (Physical Scenario):
%   一个口字型 (UI-Core) 闭合铁芯电感器，线圈绕在左柱上。
%   当施加正弦电压源时，由于铁芯材料 (Steel) 的非线性 B-H 特性，
%   随着磁通密度进入饱和区，电感急剧下降，导致励磁电流出现尖峰畸变 (奇次谐波)。
%
% 预期结果 (Expected Result):
%   1. 电流波形应呈现典型的"尖峰"状 (Peaky)，波峰因数 CF > 1.414。
%   2. 频谱分析应显示显著的 3次、5次谐波分量。
%   3. 求解过程应稳定收敛。
% =========================================================================

    % 确保求解器在路径中
    if ~exist('solve_hbfem_voltage', 'file'), addpath(genpath('src')); end

    fprintf('==============================================\n');
    fprintf('   HBFEM Voltage Driven Test (Closed Core)    \n');
    fprintf('==============================================\n');
    
    %% 1. 网格生成 (Mesh Generation)
    % 使用特制的鲁棒生成器，确保铁芯与空气在交界面共用节点。
    % 只要节点连通，磁通就能顺利通过铁芯，避免出现非物理的气隙磁阻。
    msh_file = 'test_closed_core.msh';
    create_closed_core_mesh_robust(msh_file); 
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    %% 2. 材料定义 (Material Definition)
    % 材料 1: 空气 (线性)
    Model.Materials.Lib(1).Name = 'Air'; 
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1; 
    Model.Materials.Lib(1).Sigma = 0; 
    
    % 材料 2: 软磁钢 (非线性)
    Model.Materials.Lib(2).Name = 'Robust-Steel'; 
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Mu_r = 5000; % 初始相对磁导率 (线性区)
    Model.Materials.Lib(2).Sigma = 0;   % 忽略涡流损耗，专注磁滞非线性
    
    % 定义 B-H 曲线 (典型的电工钢特性)
    % H: 磁场强度 (A/m)
    % B: 磁通密度 (T)
    % 膝点约为 1.5T (H=400 A/m)，之后进入深度饱和
    H_data = [0, 50, 100, 200, 400, 800, 1600, 5000, 10000, 50000]; 
    B_data = [0, 0.5, 0.9, 1.3, 1.5, 1.65, 1.75, 1.85, 1.95, 2.1];
    Model.Materials.Lib(2).BH_Curve.H = H_data; 
    Model.Materials.Lib(2).BH_Curve.B = B_data;
    
    % 预处理材料 (生成样条插值)
    Model = preprocess_materials(Model);
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    %% 3. 线圈定义 (Coil Setup)
    % 绕组位置: 铁芯左柱中心 [-0.2, 0, 0]
    CoilParams.Center = [-0.2, 0, 0]; 
    CoilParams.Normal = [0, 0, 1];    % 轴向沿 Z 轴
    CoilParams.Length = 0.4;          % 线圈高度
    CoilParams.Radius = 0.15;         % 半径 (需包围截面 0.1x0.1 的铁芯)
    CoilParams.N_seg = 36;            % 圆弧离散段数
    CoilParams.Turns = 100;           % 匝数 (关键参数，决定安匝数)
    CoilParams.Current = 1.0;         % 几何参考电流 (求解器会自动缩放)
    
    Coil = create_racetrack_coil(CoilParams);
    
    %% 4. 电路与求解参数 (Circuit & Solver Settings)
    Circuit.R = 0.5;         % 线圈直流电阻 (Ohm)
    Circuit.L_leak = 1e-5;   % 漏感 (H)，用于数值稳定性
    
    % 电压设定: 300V
    % 估算: V = 4.44 * f * N * A * B_max
    % 300 = 4.44 * 50 * 100 * 0.01 * B_max  => B_max ≈ 1.35 T
    % 考虑到非线性，峰值电流会激增，实际 B 场会更高，足以进入饱和区。
    Circuit.Voltage = [300; 0; 0]; 
    
    % 谐波平衡设置
    HBFEMParams.Frequency = 50;             % 基频 50Hz
    HBFEMParams.Harmonics = [1, 3, 5, 7, 9]; % 求解奇次谐波
    HBFEMParams.TimeSteps = 32;             % FFT 采样点数
    
    % 边界条件: 自动识别最外层表面并设为 Dirichlet 0
    if isfield(Raw, 'FaceTags'), fe = identify_boundary_edges(Model.Mesh, Raw, 100); else, fe = []; end
    Model.Runtime.FixedEdges = fe;
    
    %% 5. 调用求解器 (Execute Solver)
    fprintf('  [Test] Calling Solver (V=%.0f V)...\n', Circuit.Voltage(1));
    [Sol, Info] = solve_hbfem_voltage(Model, Coil, Circuit, HBFEMParams);
    
    %% 6. 结果分析 (Result Analysis)
    % 提取基波和三次谐波幅值
    I_fund = abs(Sol.Current(1)); 
    I_3rd  = abs(Sol.Current(2));
    
    fprintf('\n----------------------------------------------\n');
    fprintf('   Results Summary\n');
    fprintf('----------------------------------------------\n');
    fprintf('  Current (Fund): %.2f A\n', I_fund);
    fprintf('  Current (3rd) : %.2f A (Distortion: %.2f%%)\n', I_3rd, (I_3rd/I_fund)*100);
    fprintf('  Iterations    : %d\n', Info.Iterations);
    
    %% 7. 绘图与验证 (Visualization)
    % 初始化 AFT 模块进行时域重构
    AFT = aft_module(); 
    Info_Recon = AFT.prepare_indices(HBFEMParams.Harmonics, HBFEMParams.TimeSteps);
    
    % 构造电压行向量 (确保维度正确 1 x NumHarm)
    NumHarm = length(HBFEMParams.Harmonics);
    V_harm_vec = zeros(1, NumHarm);
    V_len = min(length(Circuit.Voltage), NumHarm);
    V_harm_vec(1:V_len) = Circuit.Voltage(1:V_len).'; 
    
    % IFFT: 频域 -> 时域
    V_t = AFT.freq2time(V_harm_vec, Info_Recon);
    I_t = AFT.freq2time(Sol.Current.', Info_Recon); 
    
    % 生成时间轴
    t_vec = linspace(0, 1/50, HBFEMParams.TimeSteps+1); t_vec(end) = [];
    
    % 图表 1: 时域波形 (V-I)
    figure('Name', 'Closed Core HBFEM Result');
    subplot(2,1,1); 
    yyaxis left; plot(t_vec*1000, V_t, 'b-', 'LineWidth', 1.5); ylabel('Voltage (V)');
    yyaxis right; plot(t_vec*1000, I_t, 'r-', 'LineWidth', 2); ylabel('Current (A)');
    xlabel('Time (ms)'); title('V-I Response (Peaky Current Indicates Saturation)'); 
    legend('Voltage', 'Current'); grid on;
    
    % 图表 2: 频域频谱 (Harmonics)
    subplot(2,1,2);
    bar(HBFEMParams.Harmonics, abs(Sol.Current), 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Harmonic Order'); ylabel('Current Magnitude (A)'); 
    title('Current Harmonics Spectrum'); grid on;
    
    % 计算波峰因数 (Crest Factor)
    % 正弦波 CF = 1.414。饱和电流 CF 通常 > 1.5。
    I_max = max(abs(I_t)); 
    CF_I = I_max / rms(I_t);
    fprintf('  Current Crest Factor: %.2f (Sine=1.41)\n', CF_I);
    
    % 自动判定测试通过
    if I_3rd/I_fund > 0.02 || CF_I > 1.5
        fprintf('[PASS] Significant distortion detected (Saturation achieved).\n');
    else
        fprintf('[WARN] Distortion low. Check magnetic circuit connectivity.\n');
    end
    
    % 清理临时文件
    if exist(msh_file,'file'), try delete(msh_file); catch; end; end
end

% -------------------------------------------------------------------------
% 辅助函数: 网格生成器 (Robust Connected Mesh)
% -------------------------------------------------------------------------
function create_closed_core_mesh_robust(filename)
% CREATE_CLOSED_CORE_MESH_ROBUST 生成闭合磁路网格
% 
% 技术细节:
%   使用统一坐标轴刻度 (Ticks) 生成点云，确保空气和铁芯在交界面处
%   拥有完全相同的节点坐标。这保证了 Delaunay 剖分后的网格是物理连通的。
%
% 结构:
%   UI型铁芯，尺寸约 0.6m x 0.8m，截面 0.1m x 0.1m。

    fprintf('  [MeshGen] Generating Connected Mesh (Shared Nodes)...\n');
    
    % 定义关键坐标刻度 (单位: m)
    x_ticks = [-0.6, -0.3, -0.1, 0.1, 0.3, 0.6];
    y_ticks = [-0.4, -0.05, 0.05, 0.4];
    z_ticks = [-0.8, -0.4, -0.3, 0.3, 0.4, 0.8];
    
    % 细分刻度以控制网格密度 (SubDiv=2 产生较粗网格，适合快速测试)
    x_grid = refine_ticks(x_ticks, 2); 
    y_grid = refine_ticks(y_ticks, 2);
    z_grid = refine_ticks(z_ticks, 2);
    
    [X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);
    P = [X(:), Y(:), Z(:)];
    nNodes = size(P, 1);
    fprintf('    - Total Nodes: %d (Optimized)\n', nNodes);
    
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList;
    
    % 标记材料区域 (根据重心坐标)
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    tol = 1e-5;
    in_outer = abs(Centers(:,1)) <= 0.3+tol & abs(Centers(:,2)) <= 0.05+tol & abs(Centers(:,3)) <= 0.4+tol;
    in_hole  = abs(Centers(:,1)) < 0.1-tol  & abs(Centers(:,3)) < 0.3-tol; 
    in_iron = in_outer & ~in_hole;
    
    ElemTags = ones(size(T, 1), 1); % Default Air
    ElemTags(in_iron) = 2;          % Iron
    
    TR = triangulation(T, P);
    FB = freeBoundary(TR);
    nFaces = size(FB, 1);
    
    % 写入 Gmsh 2.2 格式
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    fprintf(fid, '%d %.6f %.6f %.6f\n', [(1:nNodes); P']);
    fprintf(fid, '$EndNodes\n');
    
    nElems = size(T, 1);
    nTotal = nFaces + nElems;
    idx_air = find(ElemTags==1);
    idx_iron = find(ElemTags==2);
    
    % 动态计算块数量 (Surface + Air? + Iron?)
    nBlocks = 1 + ~isempty(idx_air) + ~isempty(idx_iron);
    fprintf(fid, '$Elements\n%d %d 1 %d\n', nBlocks, nTotal, nTotal);
    
    % 写入表面单元 (Tag 100)
    fprintf(fid, '2 100 2 %d\n', nFaces);
    fprintf(fid, '%d %d %d %d\n', [(1:nFaces); FB']);
    id_start = nFaces + 1;
    
    % 写入空气单元 (Tag 1)
    if ~isempty(idx_air)
        n = length(idx_air);
        fprintf(fid, '3 1 4 %d\n', n);
        fprintf(fid, '%d %d %d %d %d\n', [(id_start : id_start+n-1); T(idx_air, :)']);
        id_start = id_start + n;
    end
    
    % 写入铁芯单元 (Tag 2)
    if ~isempty(idx_iron)
        n = length(idx_iron);
        fprintf(fid, '3 2 4 %d\n', n);
        fprintf(fid, '%d %d %d %d %d\n', [(id_start : id_start+n-1); T(idx_iron, :)']);
    end
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end

function grid_out = refine_ticks(ticks, sub_divs)
    grid_out = [];
    for k = 1:length(ticks)-1
        pts = linspace(ticks(k), ticks(k+1), sub_divs+1);
        grid_out = [grid_out, pts(1:end-1)]; %#ok<AGROW>
    end
    grid_out = [grid_out, ticks(end)];
    grid_out = unique(grid_out);
end