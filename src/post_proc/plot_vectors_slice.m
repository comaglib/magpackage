function plot_vectors_slice(Model, Solution, Coil, Plane)
% PLOT_VECTORS_SLICE 在指定平面绘制 B 矢量
% Plane: 'xy', 'yz', 'zx' (过原点)

    fprintf('绘制切片矢量图 (%s)...\n', Plane);
    
    Mesh = Model.Mesh;
    P = Mesh.P; T = Mesh.T;
    
    % 1. 筛选切片附近的单元
    Centers = (P(:,T(1,:)) + P(:,T(2,:)) + P(:,T(3,:)) + P(:,T(4,:)))/4;
    
    tol = 0.02; % 切片厚度
    switch Plane
        case 'xy', mask = abs(Centers(3,:)) < tol; proj_idx=[1,2];
        case 'yz', mask = abs(Centers(1,:)) < tol; proj_idx=[2,3];
        case 'zx', mask = abs(Centers(2,:)) < tol; proj_idx=[3,1];
    end
    
    target_elems = find(mask);
    
    % 2. 计算 B
    % 使用 compute_biot_savart_B_serial 和 calc_element_B
    % 串行计算即可
    
    numT = length(target_elems);
    if numT == 0, warning('No elements on slice.'); return; end
    
    Pts = Centers(:, target_elems);
    B_vecs = zeros(3, numT);
    
    % Calc Bs
    Bs = compute_biot_savart_B_serial(Coil, Pts);
    
    % Calc Br
    for k = 1:numT
        e = target_elems(k);
        nodes = T(:,e);
        pts_e = P(:,nodes);
        edges = Model.Mesh.T2E(:,e);
        s = Model.Mesh.T2E_Sign(:,e);
        a_v = Solution.A(edges) .* double(s);
        
        Br = calc_element_B(pts_e, a_v);
        B_vecs(:, k) = Br + Bs(:, k);
    end
    
    % 3. Quiver Plot
    figure;
    quiver3(Pts(1,:), Pts(2,:), Pts(3,:), ...
            B_vecs(1,:), B_vecs(2,:), B_vecs(3,:));
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['B Field Vectors on ' Plane ' Plane']);
end