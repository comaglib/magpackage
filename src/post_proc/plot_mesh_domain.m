function plot_mesh_domain(Mesh, TargetRegions, Options)
% PLOT_MESH_DOMAIN 绘制指定域的网格 (仅显示外表面)
%
% 功能:
%   提取指定 RegionTags 对应的四面体，计算其外包络面 (Hull)，并绘制。
%   这使得在查看复杂装配体内部结构（如线圈、铁心）时非常清晰。
%
% 输入:
%   Mesh          - 网格结构体 (含 .P, .T, .RegionTags)
%   TargetRegions - 需要绘制的区域 ID 列表 (向量)
%   Options       - (可选) 绘图选项结构体
%       .FaceColor - 颜色 [r g b] 或 'r', 'b' 等
%       .FaceAlpha - 透明度 (0-1)
%       .EdgeColor - 边框颜色 ('none' 或 'k')
%       .Title     - 标题
%       .NewFigure - 是否新建图形窗口 (true/false)

    if nargin < 3, Options = struct(); end
    
    % 默认选项
    if ~isfield(Options, 'FaceColor'), Options.FaceColor = [0.8 0.8 1.0]; end
    % [修改] 默认设置为半透明，便于观察内部重叠
    if ~isfield(Options, 'FaceAlpha'), Options.FaceAlpha = 0.6; end 
    if ~isfield(Options, 'EdgeColor'), Options.EdgeColor = 'k'; end
    if ~isfield(Options, 'NewFigure'), Options.NewFigure = true; end
    
    P = Mesh.P;
    T = Mesh.T;
    Tags = Mesh.RegionTags;
    
    % 1. 筛选单元
    % 找出属于 TargetRegions 的所有单元索引
    mask = ismember(Tags, TargetRegions);
    T_sub = T(:, mask);
    
    if isempty(T_sub)
        warning('未找到属于 Region [%s] 的单元。', num2str(TargetRegions));
        return;
    end
    
    numElems = size(T_sub, 2);
    fprintf('正在绘制域 [%s]: %d 个四面体...\n', num2str(TargetRegions), numElems);
    
    % 2. 提取所有面 (4 * Ne)
    % 四面体面定义 (Local ID)
    f1 = T_sub([1, 2, 3], :);
    f2 = T_sub([1, 3, 4], :);
    f3 = T_sub([1, 4, 2], :);
    f4 = T_sub([2, 4, 3], :);
    
    AllFaces = [f1, f2, f3, f4]; % 3 x (4*Ne)
    
    % 3. 核心算法: 提取外表面 (只出现一次的面)
    % 转置为 (4*Ne) x 3 以利用 unique('rows')
    AllFacesT = sort(AllFaces, 1)'; 
    
    % 使用累积计数法找出唯一面
    [C, ~, ic] = unique(AllFacesT, 'rows');
    counts = accumarray(ic, 1);
    
    % 提取只出现 1 次的面索引 (即外边界)
    BoundaryFaces = C(counts == 1, :);
    
    fprintf('  - 提取到 %d 个外表面三角形。\n', size(BoundaryFaces, 1));
    
    % 4. 绘图
    if Options.NewFigure
        figure;
    else
        hold on;
    end
    
    % 绘制表面
    trisurf(BoundaryFaces, P(1,:), P(2,:), P(3,:), ...
        'FaceColor', Options.FaceColor, ...
        'FaceAlpha', Options.FaceAlpha, ...
        'EdgeColor', Options.EdgeColor, ...
        'EdgeAlpha', 0.3); % 边框稍微淡一点，避免太杂乱
        
    % [修改] 标准工程视图设置
    axis equal; 
    axis tight; 
    axis vis3d;
    grid off;
    view(3);
    camlight; % 默认灯光，提供基本的立体感，但不使用 gouraud 插值，保持网格硬朗感
    
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    
    if isfield(Options, 'Title')
        title(Options.Title);
    else
        title(['Mesh Domain: ' num2str(TargetRegions)]);
    end
    
    hold off;
end