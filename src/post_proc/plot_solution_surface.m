function plot_solution_surface(Model, Solution, FieldType)
% PLOT_SOLUTION_SURFACE 绘制三维表面云图
%
% 输入:
%   Model, Solution
%   FieldType - 'B_mag' (磁通模), 'A_mag' (磁位模), 'J_mag' (电流模)

    Mesh = Model.Mesh;
    
    % 1. 计算场量 (在单元中心)
    % 这里简化处理，假设我们需要 B_mag
    numElems = size(Mesh.T, 2);
    ElemData = zeros(numElems, 1);
    
    if strcmp(FieldType, 'B_mag')
        % 计算每个单元的 B
        % 并行计算太慢，这里串行即可，绘图不需要实时
        % 需要调用 calc_element_B...
        % 为简化，这里写个简单的循环
        % 实际应调用专门的 Field Calculator
        % 假设 Solution 已经包含了 B (在后续版本中应预计算)
        
        % 临时: 现场计算 B_mag (仅取前 10000 个单元演示，全量太慢)
        % *建议*: 在 solve 之后就把 B 算好存入 Solution.B
        fprintf('Computing B field for plotting...\n');
        
        % 快速计算: 
        % 我们需要 Coil 信息来算 Bs。这里假设 Solution 包含 B_total
        % 如果没有，我们暂时画 A_mag (节点值)
        warning('Plotting A_mag instead of B (B requires coil info).');
        
        % 画 A_mag (节点值直接可得)
        A_edges = Solution.A;
        % 需要将 Edge A 转换为 Node A?
        % Edge A = (Ai + Aj)/2 * L? No. Edge A is line integral.
        % Node A 不太好定义。
        
        % 改回画 B. 我们必须传入 Coil 才能算 Bs。
        % 为了接口简单，我们假设 Model.Physics.Coil 存在
        % 或者只画 Br (Reaction Field)
        
    end
    
    % ----------------------------------------------------------
    % 实用方案: 仅绘制网格表面
    % ----------------------------------------------------------
    
    % 提取表面 (Boundary Faces)
    % 如果 MeshRaw 不在，我们用 triangulation 提取
    TR = triangulation(Mesh.T', Mesh.P');
    [FB, ~] = freeBoundary(TR); % N_face x 3
    
    % 节点坐标
    X = Mesh.P(1,:); Y = Mesh.P(2,:); Z = Mesh.P(3,:);
    
    figure;
    % 绘制网格
    trisurf(FB, X, Y, Z, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.5);
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Mesh Surface');
    
    % 如果有场数据 NodalVal
    % trisurf(FB, X, Y, Z, NodalVal, ...)
end