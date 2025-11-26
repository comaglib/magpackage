function CoilSegments = create_racetrack_coil(Params)
% CREATE_RACETRACK_COIL 生成跑道型线圈的离散线段数据 (修正版)
% 修复了起点和终点重合导致的零长度线段问题

    L = Params.Length;
    R = Params.Radius;
    N_arc = Params.N_seg;
    I_val = Params.Current;
    
    Nodes = [];
    
    % --- 直线 1 (底部) ---
    % 起点: [-L/2, -R] -> 终点: [L/2, -R]
    Nodes(:, end+1) = [-L/2; -R; 0];
    Nodes(:, end+1) = [ L/2; -R; 0];
    
    % --- 圆弧 1 (右侧) ---
    % 从 -pi/2 到 pi/2
    theta1 = linspace(-pi/2, pi/2, N_arc+1);
    theta1 = theta1(2:end); % 去掉起点(已在直线中)
    for t = theta1
        Nodes(:, end+1) = [L/2 + R*cos(t); R*sin(t); 0]; %#ok<AGROW>
    end
    
    % --- 直线 2 (顶部) ---
    % 此时我们在 [L/2, R]. 下一个点是 [-L/2, R]
    Nodes(:, end+1) = [-L/2; R; 0]; 
    
    % --- 圆弧 2 (左侧) ---
    % 从 pi/2 到 3pi/2
    % 3pi/2 对应的坐标是 [-L/2, -R]，这正是 Nodes(:,1)
    % [关键修正] 使用 2:end-1，去掉起点(已在直线中) 和 终点(即 Nodes(1))
    theta2 = linspace(pi/2, 3*pi/2, N_arc+1);
    theta2 = theta2(2:end-1); 
    
    for t = theta2
        Nodes(:, end+1) = [-L/2 + R*cos(t); R*sin(t); 0]; %#ok<AGROW>
    end
    
    % 闭合线圈
    % 此时 Nodes(end) 是左侧圆弧倒数第二个点
    % Nodes(1) 是左侧圆弧的终点(也是整个线圈起点)
    % P2 连接它们，自然形成了最后一段圆弧线段
    
    num_nodes = size(Nodes, 2);
    P1 = Nodes;
    P2 = [Nodes(:, 2:end), Nodes(:, 1)]; % 循环连接
    
    % 2. 坐标变换: 局部 -> 全局
    if norm(Params.Normal) < 1e-12
         w_axis = [0;0;1];
    else
         w_axis = Params.Normal(:) / norm(Params.Normal);
    end
    
    aux = [0; 0; 1];
    if abs(dot(w_axis, aux)) > 0.99
        aux = [0; 1; 0];
    end
    u_axis = cross(aux, w_axis); u_axis = u_axis / norm(u_axis);
    v_axis = cross(w_axis, u_axis);
    
    RotMat = [u_axis, v_axis, w_axis]; 
    Center = Params.Center(:);
    
    P1_global = RotMat * P1 + Center;
    P2_global = RotMat * P2 + Center;
    
    CoilSegments.P1 = P1_global;
    CoilSegments.P2 = P2_global;
    CoilSegments.I = repmat(I_val, 1, num_nodes);
end