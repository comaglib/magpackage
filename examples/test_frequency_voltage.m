function test_frequency_voltage()
% TEST_FREQUENCY_VOLTAGE 测试频域电压驱动求解器 (一致性验证版)
% 
% 验证策略:
%   1. 电压驱动 (Method B): 给定 U -> 求解得到 I_B 和 A_B
%   2. 电流驱动 (Method A): 给定 I = I_B -> 求解得到 A_A
%   3. 比较 A_A 和 A_B。如果一致，证明场路耦合矩阵构造正确。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Field-Circuit Coupling Test (Consistency)  \n');
    fprintf('==============================================\n');
    
    % 1. 准备模型
    msh_file = 'test_voltage_delaunay.msh';
    create_test_mesh(msh_file); 
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 材料 (空气 + 铁心)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0; Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Iron';
    Model.Materials.Lib(2).Type = 'Linear';
    Model.Materials.Lib(2).Mu_r = 1000.0; Model.Materials.Lib(2).Sigma = 0; 
    
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.4; CoilParams.Radius = 0.3; 
    CoilParams.N_seg = 72;
    CoilParams.Current = 0; 
    Coil = create_racetrack_coil(CoilParams);
    
    % 边界
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    if isempty(fixed_edges), warning('Boundary not found, using natural BC.'); end
    Model.Runtime.FixedEdges = fixed_edges;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 测试参数
    freq = 50;
    FreqParams.Frequency = freq;
    
    % ---------------------------------------------------------
    % Step 1: 电压驱动求解 (Voltage Driven)
    % ---------------------------------------------------------
    fprintf('\n--- Step 1: Voltage Driven Solve (U=10V) ---\n');
    
    U_target = 10.0;
    CircuitProps.R = 1e-3; % 给一点微小电阻防止奇异 (如果纯感性)
    CircuitProps.L_leak = 0;
    CircuitProps.U = U_target;
    
    [Sol_Voltage, ~] = solve_frequency_voltage(Model, Coil, CircuitProps, FreqParams);
    
    I_calculated = Sol_Voltage.I;
    A_voltage = Sol_Voltage.A;
    
    fprintf('  -> Calculated Current I = %.4f + j%.4f A\n', real(I_calculated), imag(I_calculated));
    
    % ---------------------------------------------------------
    % Step 2: 电流驱动验证 (Current Driven)
    % ---------------------------------------------------------
    fprintf('\n--- Step 2: Current Driven Verification (I=I_calc) ---\n');
    
    Coil_Check = Coil;
    Coil_Check.I = I_calculated; % 将上一步算出的电流代入
    
    [Sol_Current, ~] = solve_frequency_linear(Model, Coil_Check, FreqParams);
    
    A_current = Sol_Current.A;
    
    % ---------------------------------------------------------
    % Step 3: 一致性对比
    % ---------------------------------------------------------
    fprintf('\n--- Consistency Check ---\n');
    
    % 比较两个解向量 A
    diff_norm = norm(A_voltage - A_current);
    ref_norm = norm(A_voltage);
    rel_err = diff_norm / max(ref_norm, 1e-10);
    
    fprintf('|A_voltage| = %.4e\n', ref_norm);
    fprintf('|A_current| = %.4e\n', norm(A_current));
    fprintf('Relative Error: %.4e\n', rel_err);
    
    if rel_err < 1e-5
        fprintf('[PASS] Voltage and Current solvers are mathematically consistent.\n');
        fprintf('       (Coupling Matrix W and Source Vector Q are correct transpose pairs)\n');
    else
        error('[FAIL] Solvers are inconsistent. Check coupling matrices.');
    end
end

function create_test_mesh(filename)
    % 生成包含铁心和空气的鲁棒网格
    % 范围: [-0.5, 0.5]
    
    % 点云
    [X1, Y1, Z1] = ndgrid(linspace(-0.15, 0.15, 5));
    P_inner = [X1(:), Y1(:), Z1(:)];
    [X2, Y2, Z2] = ndgrid(linspace(-0.5, 0.5, 6));
    P_outer = [X2(:), Y2(:), Z2(:)];
    P = unique([P_inner; P_outer], 'rows');
    nNodes = size(P, 1);
    
    % Delaunay
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList; 
    
    % Chirality Fix
    p1 = P(T(:,1), :); p2 = P(T(:,2), :); p3 = P(T(:,3), :); p4 = P(T(:,4), :);
    vols = sum(cross(p2-p1, p3-p1, 2) .* (p4-p1), 2);
    if any(vols < 0), T(vols < 0, [1 2]) = T(vols < 0, [2 1]); end
    
    P_io = P'; T_io = T'; nElems = size(T_io, 2);
    
    % Tags
    Centers = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:) + P(T(:,4),:)) / 4;
    in_iron = abs(Centers(:,1)) < 0.16 & abs(Centers(:,2)) < 0.16 & abs(Centers(:,3)) < 0.16;
    ElemTags = ones(1, nElems);
    ElemTags(in_iron) = 2;
    
    % Boundary
    TR = triangulation(T, P);
    FB = freeBoundary(TR); BoundaryFaces = FB'; nFaces = size(BoundaryFaces, 2);
    
    % Write
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 %d 1 %d\n2 1 0 %d\n', nNodes, nNodes, nNodes);
    fprintf(fid, '%d %.6f %.6f %.6f\n', [(1:nNodes); P_io]);
    fprintf(fid, '$EndNodes\n');
    
    nTotal = nFaces + nElems;
    idx_air = (ElemTags == 1); idx_iron = (ElemTags == 2);
    nAir = sum(idx_air); nIron = sum(idx_iron);
    
    fprintf(fid, '$Elements\n%d %d 1 %d\n', 3, nTotal, nTotal);
    fprintf(fid, '2 100 2 %d\n', nFaces);
    fprintf(fid, '%d %d %d %d\n', [(1:nFaces); BoundaryFaces]);
    id = nFaces + 1;
    
    if nAir > 0
        fprintf(fid, '3 1 4 %d\n', nAir);
        fprintf(fid, '%d %d %d %d %d\n', [(id:id+nAir-1); T_io(:, idx_air)]);
        id = id + nAir;
    end
    if nIron > 0
        fprintf(fid, '3 2 4 %d\n', nIron);
        fprintf(fid, '%d %d %d %d %d\n', [(id:id+nIron-1); T_io(:, idx_iron)]);
    end
    
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end