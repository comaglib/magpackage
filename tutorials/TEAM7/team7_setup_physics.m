function [Model, Coil] = team7_setup_physics(ModelIn)
% TEAM7_SETUP_PHYSICS [修复版 - 奇异性消除]
% 
% 修正:
% 1. [关键] 给空气赋予微小电导率 (1.0 S/m)，防止 V 场方程奇异。
% 2. 保持 Coil 积分点设置 (4x4)。

    % 路径检查
    if ~exist('analyze_coil_geometry', 'file'), addpath(genpath('../../src')); end
    
    fprintf('==============================================\n');
    fprintf('   Step 2: Physics Setup (Regularized V)      \n');
    fprintf('==============================================\n');

    %% 1. 网格处理
    if nargin > 0 && isfield(ModelIn, 'Mesh')
        Model = ModelIn;
        MeshRaw = Model.MeshRaw; 
        if ~isfield(MeshRaw, 'Faces')
             mesh_file = 'team7.mphtxt';
             MeshRaw = read_mphtxt(mesh_file, 'mm');
        end
        fprintf('  [Mesh] 使用外部传入的网格数据。\n');
    else
        mesh_file = 'team7.mphtxt';
        if ~isfile(mesh_file)
            curr_path = fileparts(mfilename('fullpath'));
            mesh_file = fullfile(curr_path, 'team7.mphtxt');
        end
        MeshRaw = read_mphtxt(mesh_file, 'mm'); 
        Model.Mesh = build_topology(MeshRaw);
        Model.MeshRaw = MeshRaw; 
    end
    
    max_coord = max(abs(Model.Mesh.P(:)));
    fprintf('  [Unit] 坐标范围检查: Max = %.3f m\n', max_coord);

    %% 2. 材料定义 (修复奇异性)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0; 
    
    % [CRITICAL FIX] 空气电导率不能为 0，否则 V 场方程奇异。
    % 设置一个较小的值 (相对于铝的 3e7，1.0 是可以忽略的，但能稳定矩阵)
    Model.Materials.Lib(1).Sigma = 1.0; 
    
    Model.Materials.Lib(2).Name = 'Aluminum';
    Model.Materials.Lib(2).Type = 'Linear';
    Model.Materials.Lib(2).Mu_r = 1.0;
    Model.Materials.Lib(2).Sigma = 3.526e7; 
    
    mat_map = ones(1, max(Model.Mesh.RegionTags)); 
    if ismember(2, Model.Mesh.RegionTags), mat_map(2) = 2; end
    Model.Materials.ActiveMap = mat_map(Model.Mesh.RegionTags);
    
    %% 3. 线圈定义
    fprintf('  [Coil] 分析线圈几何...\n');
    TotalAmpTurns = 2742; 
    CoilParams = analyze_coil_geometry(Model.Mesh, [3, 4, 5, 6]); 
    CoilParams.Current = TotalAmpTurns; 

    % 使用 4x4 积分点
    IntegrationOrder = 4; 
    fprintf('  [Coil] 生成源场积分点 (Gauss Order: %dx%d)...\n', IntegrationOrder, IntegrationOrder);
    
    Coil = create_thick_rounded_rectangle_coil(CoilParams, IntegrationOrder);
    
    num_integration_points = size(Coil.P1, 2); 
    current_weight = TotalAmpTurns / num_integration_points;
    Coil.I = repmat(current_weight, size(Coil.I, 1), num_integration_points);
    
    fprintf('         总积分段数: %d, 单段权重电流: %.4f A\n', num_integration_points, current_weight);

    %% 4. 边界条件
    fprintf('  [BC] 正在检测外边界 (Tag Analysis)...\n');
    outer_tags = detect_outer_boundary_tags(MeshRaw);
    
    if isempty(outer_tags)
        warning('未检测到外边界 Tag，使用默认空边界。');
        Model.Runtime.FixedEdges = [];
    else
        fprintf('       识别到外边界 Tags: [%s]\n', num2str(outer_tags));
        Model.Runtime.FixedEdges = identify_boundary_edges(Model.Mesh, MeshRaw, outer_tags);
        fprintf('       已施加 Dirichlet BC (A=0) 于 %d 条棱。\n', length(Model.Runtime.FixedEdges));
    end
    
    fprintf('  [Step 2] 设置完成。\n');
end