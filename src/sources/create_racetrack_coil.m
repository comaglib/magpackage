function CoilSegments = create_racetrack_coil(Params)
% CREATE_RACETRACK_COIL 生成跑道型线圈的离散线段数据 (支持匝数)
% 
% 输入 Params:
%   .Turns - 匝数 (默认1)
%   ... (其他同前)

    % 默认匝数
    if isfield(Params, 'Turns')
        Turns = Params.Turns;
    else
        Turns = 1.0;
    end

    L = Params.Length;
    R = Params.Radius;
    N_arc = Params.N_seg;
    I_val = Params.Current;
    
    Nodes = [];
    
    % --- 直线 1 ---
    Nodes(:, end+1) = [-L/2; -R; 0];
    Nodes(:, end+1) = [ L/2; -R; 0];
    
    % --- 圆弧 1 ---
    theta1 = linspace(-pi/2, pi/2, N_arc+1);
    theta1 = theta1(2:end); 
    for t = theta1
        Nodes(:, end+1) = [L/2 + R*cos(t); R*sin(t); 0]; %#ok<AGROW>
    end
    
    % --- 直线 2 ---
    Nodes(:, end+1) = [-L/2; R; 0]; 
    
    % --- 圆弧 2 ---
    theta2 = linspace(pi/2, 3*pi/2, N_arc+1);
    theta2 = theta2(2:end-1); 
    for t = theta2
        Nodes(:, end+1) = [-L/2 + R*cos(t); R*sin(t); 0]; %#ok<AGROW>
    end
    
    num_nodes = size(Nodes, 2);
    P1 = Nodes;
    P2 = [Nodes(:, 2:end), Nodes(:, 1)]; 
    
    % 坐标变换
    if norm(Params.Normal) < 1e-12, w_axis = [0;0;1]; else, w_axis = Params.Normal(:) / norm(Params.Normal); end
    aux = [0; 0; 1]; if abs(dot(w_axis, aux)) > 0.99, aux = [0; 1; 0]; end
    u_axis = cross(aux, w_axis); u_axis = u_axis / norm(u_axis);
    v_axis = cross(w_axis, u_axis);
    RotMat = [u_axis, v_axis, w_axis]; 
    Center = Params.Center(:);
    
    P1_global = RotMat * P1 + Center;
    P2_global = RotMat * P2 + Center;
    
    CoilSegments.P1 = P1_global;
    CoilSegments.P2 = P2_global;
    CoilSegments.I = repmat(I_val, 1, num_nodes);
    
    % [新增] 保存匝数
    CoilSegments.Turns = Turns;
end