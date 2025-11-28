function Coil = create_rounded_rectangle_coil(CoilParams)
% CREATE_ROUNDED_RECTANGLE_COIL 生成圆角矩形/跑道形线圈几何 (Biot-Savart源)
%
% 输入:
%   CoilParams - 结构体，包含:
%       .Center    : [x; y; z] 线圈几何中心 (米)
%       .Lx_mean   : X方向平均跨度 (中心到中心 + 直径) (米)
%                    即线圈红色路径的BoundingBox X宽度
%       .Ly_mean   : Y方向平均跨度 (米)
%       .R_mean    : 圆角平均半径 (米)
%       .Current   : 电流 (安培)
%
% 输出:
%   Coil       - 结构体，包含 .P1, .P2, .I, .Turns

    % 提取参数
    cx = CoilParams.Center(1); 
    cy = CoilParams.Center(2); 
    cz = CoilParams.Center(3);
    
    Lx = CoilParams.Lx_mean;
    Ly = CoilParams.Ly_mean;
    R  = CoilParams.R_mean;
    
    % 几何校验: 半径不能超过短边的一半
    max_R = min(Lx, Ly) / 2;
    if R > max_R
        warning('create_rounded_rectangle_coil:RadiusTooLarge', ...
            'R_mean (%.2f mm) 大于短边的一半，已自动限制为 %.2f mm 以形成跑道形。', ...
            R*1000, max_R*1000);
        R = max_R;
    end
    
    % 计算直线段的一半长度 (Half Length of Straight section)
    % 几何关系: Total_Span = 2*R + 2*dx  =>  dx = Span/2 - R
    dx = Lx/2 - R;
    dy = Ly/2 - R;
    
    % 防止浮点微小负数
    dx = max(0, dx);
    dy = max(0, dy);
    
    Nodes = []; 
    N_arc = 12; % 90度圆弧的分段数
    
    % --- 逆时针生成节点 ---
    
    % 1. 右侧直线 (Right Edge): (cx+Lx/2, cy-dy) -> (cx+Lx/2, cy+dy)
    %    起始点位于圆角结束处
    Nodes = [Nodes, [cx + dx + R; cy - dy; cz]]; 
    Nodes = [Nodes, [cx + dx + R; cy + dy; cz]];
    
    % 2. 右上角 (Top-Right): Center (cx+dx, cy+dy)
    t = linspace(0, pi/2, N_arc+1); t(1)=[];
    for i=1:length(t)
        Nodes = [Nodes, [cx + dx + R*cos(t(i)); cy + dy + R*sin(t(i)); cz]];
    end
    
    % 3. 上侧直线 (Top Edge)
    Nodes = [Nodes, [cx - dx; cy + dy + R; cz]];
    
    % 4. 左上角 (Top-Left): Center (cx-dx, cy+dy)
    t = linspace(pi/2, pi, N_arc+1); t(1)=[];
    for i=1:length(t)
        Nodes = [Nodes, [cx - dx + R*cos(t(i)); cy + dy + R*sin(t(i)); cz]];
    end
    
    % 5. 左侧直线 (Left Edge)
    Nodes = [Nodes, [cx - dx - R; cy - dy; cz]];
    
    % 6. 左下角 (Bottom-Left): Center (cx-dx, cy-dy)
    t = linspace(pi, 3*pi/2, N_arc+1); t(1)=[];
    for i=1:length(t)
        Nodes = [Nodes, [cx - dx + R*cos(t(i)); cy - dy + R*sin(t(i)); cz]];
    end
    
    % 7. 下侧直线 (Bottom Edge)
    Nodes = [Nodes, [cx + dx; cy - dy - R; cz]];
    
    % 8. 右下角 (Bottom-Right): Center (cx+dx, cy-dy)
    t = linspace(3*pi/2, 2*pi, N_arc+1); t(1)=[];
    for i=1:length(t)
        Nodes = [Nodes, [cx + dx + R*cos(t(i)); cy - dy + R*sin(t(i)); cz]];
    end
    
    % 构建输出
    Coil.P1 = Nodes;
    Coil.P2 = [Nodes(:,2:end), Nodes(:,1)]; % 闭合回路
    Coil.I = repmat(CoilParams.Current, 1, size(Nodes, 2));
    Coil.Turns = 1; 
end