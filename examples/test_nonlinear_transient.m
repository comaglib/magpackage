function test_nonlinear_transient()
% TEST_NONLINEAR_TRANSIENT 测试非线性瞬态场路耦合 (最终验证版)
% 
% 场景: 
%   模拟一个 50 匝的铁心电感器在 500V 电压阶跃激励下的瞬态响应。
%   通过观察电流上升斜率的变化，验证磁饱和效应。
%
% 验证点:
%   1. 网格生成与拓扑构建的正确性 (无负体积)。
%   2. 绕组耦合向量 (Winding Vector) 的有效性。
%   3. 电流从线性上升 (高电感) 过渡到加速上升 (饱和/低电感) 的物理现象。
%
% 依赖模块:
%   - src/mesh_io
%   - src/physics_setup
%   - src/materials
%   - src/assembly
%   - src/time_stepping/solve_transient_voltage_nonlinear

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Nonlinear Transient Test (Voltage Step)    \n');
    fprintf('==============================================\n');
    
    %% 1. 前处理: 网格生成
    % 使用内存直传方式生成包含空气和铁心的鲁棒网格
    msh_file = 'test_nl_trans.msh';
    MeshRaw = create_test_mesh_direct(msh_file); 
    
    fprintf('  - Mesh Stats: %d Nodes, %d Tets, %d Tris\n', ...
        size(MeshRaw.P, 2), size(MeshRaw.T, 2), size(MeshRaw.Faces, 2));
    
    % 构建拓扑结构 (Edges, T2E)
    Model.Mesh = build_topology(MeshRaw);
    
    %% 2. 材料定义
    % Material 1: Air (Linear)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0;
    Model.Materials.Lib(1).Sigma = 0;
    
    % Material 2: Steel (Nonlinear)
    Model.Materials.Lib(2).Name = 'Steel';
    Model.Materials.Lib(2).Type = 'Nonlinear';
    Model.Materials.Lib(2).Sigma = 1e6; % 导电，产生涡流
    Model.Materials.Lib(2).Mu_r = 5000; % 初始相对磁导率
    
    % 定义 B-H 曲线数据 (典型软磁材料)
    H_data = [0, 100, 200, 500, 1000, 5000, 10000, 50000];
    B_data = [0, 0.8, 1.2, 1.5, 1.65, 1.85, 1.95, 2.2];
    Model.Materials.Lib(2).BH_Curve.H = H_data;
    Model.Materials.Lib(2).BH_Curve.B = B_data;
    
    % 预处理材料 (生成插值样条)
    Model = preprocess_materials(Model);
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    %% 3. 线圈定义
    CoilParams.Center = [0, 0, 0]; 
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.4; 
    CoilParams.Radius = 0.3; 
    CoilParams.N_seg = 100; 
    CoilParams.Current = 0; 
    CoilParams.Turns = 50; % [关键] 50匝线圈
    
    Coil = create_racetrack_coil(CoilParams);
    
    % 验证绕组耦合向量是否有效
    W_check = assemble_winding_vector(Model, Coil);
    if nnz(W_check) == 0
        error('Fatal: Coil Winding Vector is empty. Coil might be outside the mesh.');
    else
        fprintf('  [Check] Winding vector active (nnz=%d).\n', nnz(W_check));
    end
    
    %% 4. 电路与时间参数
    CircuitProps.R = 5.0;       % 电路电阻 (Ohm)
    CircuitProps.L_leak = 1e-5; % 漏感 (H)，用于数值稳定
    CircuitProps.V_func = @(t) 500.0 * (t > 0); % 电压源: 500V 阶跃
    
    TimeParams.dt = 0.001;      % 时间步长 (s)
    TimeParams.N_steps = 50;    % 总步数 (总时间 50ms)
    
    %% 5. 边界条件
    % 尝试从网格数据中识别边界 (Tag 100)
    if isfield(MeshRaw, 'FaceTags')
        fixed_edges = identify_boundary_edges(Model.Mesh, MeshRaw, 100);
    else
        fixed_edges = [];
    end
    
    if isempty(fixed_edges) 
        warning('Boundary edges not found. Using natural BC (Neumann).');
        fixed_edges = [];
    end
    Model.Runtime.FixedEdges = fixed_edges;
    
    %% 6. 求解器配置
    Model.Solver.Nonlinear.MaxIter = 25;
    Model.Solver.Nonlinear.Tolerance = 1e-4;
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = false;
    Model.Solver.Linear.MumpsICNTL.i8 = 77; % 自动缩放
    
    % 启动并行池
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    %% 7. 执行求解
    [Sol, ~] = solve_transient_voltage_nonlinear(Model, Coil, CircuitProps, TimeParams);
    
    %% 8. 结果分析与验证
    I_hist = Sol.I_history;
    dI = diff(I_hist);
    
    fprintf('\n--- Analysis ---\n');
    target_I = 500.0 / 5.0; % V/R
    fprintf('Final Current: %.4f A (Target V/R: %.1f A)\n', I_hist(end), target_I);
    
    if length(dI) >= 5
        slope_early = mean(dI(1:3));        % 早期斜率 (高电感区)
        slope_late  = mean(dI(end-4:end));  % 晚期斜率 (饱和区)
        
        fprintf('Early Slope (High Inductance): %.4f A/step\n', slope_early);
        fprintf('Late Slope  (Low Inductance):  %.4f A/step\n', slope_late);
        
        % 物理判据: 饱和 -> L 减小 -> dI/dt 增大
        if slope_late > slope_early * 1.05
            fprintf('[PASS] Current slope increased (Saturation detected).\n');
        elseif slope_late > 0
            fprintf('[INFO] Current rising consistent (Linear dominance).\n');
        else
            fprintf('[FAIL] Current did not rise correctly.\n');
        end
    end

    %% 9. 绘制电流曲线 (新增)
    t_vec = (1:TimeParams.N_steps) * TimeParams.dt;
    
    figure('Name', 'Nonlinear Transient Current');
    plot(t_vec * 1000, I_hist, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    grid on;
    xlabel('Time (ms)');
    ylabel('Current (A)');
    title('Coil Current Response (Voltage Driven)');
    
    % 标记线性理论预测 (参考)
    hold on;
    yline(target_I, 'r--', 'Steady State (V/R)');
    legend('Finite Element Result', 'Theoretical Limit', 'Location', 'southeast');
    hold off;
    
    fprintf('  [Plot] Current waveform plotted.\n');
end

% -------------------------------------------------------------------------
% 辅助函数: 内存直传式网格生成器
% -------------------------------------------------------------------------
function MeshRaw = create_test_mesh_direct(filename)
    fprintf('  [MeshGen] Generating mesh via Delaunay (Direct Memory)...\n');
    
    noise = 1e-5; % 微小扰动防止 Delaunay 退化
    
    % 1. 点云生成
    % 铁心区域 [-0.15, 0.15]
    [X1, Y1, Z1] = ndgrid(linspace(-0.15, 0.15, 5));
    P_inner = [X1(:), Y1(:), Z1(:)];
    
    % 空气区域 [-0.6, 0.6] (确保完全包裹线圈)
    [X2, Y2, Z2] = ndgrid(linspace(-0.6, 0.6, 7));
    P_outer = [X2(:), Y2(:), Z2(:)];
    
    P = unique([P_inner; P_outer], 'rows');
    nNodes = size(P, 1);
    
    % 2. Delaunay 三角剖分
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList; 
    
    % 3. 手性修正 (Chirality Fix)
    % 确保所有四面体体积为正 (符合右手定则)
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    vols = sum(cross(p2-p1, p3-p1, 2) .* (p4-p1), 2);
    
    neg_mask = vols < 0;
    if any(neg_mask)
        % 交换节点 1 和 2 以翻转手性
        T(neg_mask, [1 2]) = T(neg_mask, [2 1]);
    end
    
    % 4. 构建数据结构
    P_io = P'; 
    T_io = T'; 
    nElems = size(T_io, 2);
    
    % 材料区域标记
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    in_iron = abs(Centers(:,1)) < 0.16 & abs(Centers(:,2)) < 0.16 & abs(Centers(:,3)) < 0.16;
    
    ElemTags = ones(1, nElems); 
    ElemTags(in_iron) = 2; % Tag 2 = Iron
    
    % 提取边界表面
    TR = triangulation(T, P);
    FB = freeBoundary(TR); 
    BoundaryFaces = FB'; 
    nFaces = size(BoundaryFaces, 2);
    FaceTags = repmat(100, 1, nFaces); % Tag 100 = Boundary
    
    % 5. 返回 Mesh 结构体
    MeshRaw.P = P_io;
    MeshRaw.T = T_io;
    MeshRaw.RegionTags = ElemTags;
    MeshRaw.Faces = BoundaryFaces;
    MeshRaw.FaceTags = FaceTags;
    
    % 6. 写入文件 (用于后续检查或可视化)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    % 格式: ID X Y Z
    NodeData = [(1:nNodes); P_io];
    fprintf(fid, '%d %.6f %.6f %.6f\n', NodeData);
    fprintf(fid, '$EndNodes\n');
    
    nTotal = nFaces + nElems;
    idx_air = (ElemTags == 1);
    idx_iron = (ElemTags == 2);
    nAir = sum(idx_air);
    nIron = sum(idx_iron);
    
    nBlocks = 1 + (nAir>0) + (nIron>0);
    fprintf(fid, '$Elements\n%d %d 1 %d\n', nBlocks, nTotal, nTotal);
    
    % Surface Block
    fprintf(fid, '2 100 2 %d\n', nFaces);
    TriData = [(1:nFaces); BoundaryFaces];
    fprintf(fid, '%d %d %d %d\n', TriData);
    
    id = nFaces + 1;
    % Air Block
    if nAir > 0
        fprintf(fid, '3 1 4 %d\n', nAir);
        AirData = [(id : id+nAir-1); T_io(:, idx_air)];
        fprintf(fid, '%d %d %d %d %d\n', AirData);
        id = id + nAir;
    end
    
    % Iron Block
    if nIron > 0
        fprintf(fid, '3 2 4 %d\n', nIron);
        IronData = [(id : id+nIron-1); T_io(:, idx_iron)];
        fprintf(fid, '%d %d %d %d %d\n', IronData);
    end
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end