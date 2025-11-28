function [val, grad] = lagrange_tet_p1(q_pts)
% LAGRANGE_TET_P1 计算 Lagrange P1 基函数及其参考梯度
% 输入:
%   q_pts - [3 x N_q] 参考坐标 (xi, eta, zeta)
% 输出:
%   val   - [4 x N_q] 基函数值
%   grad  - [4 x 3 x N_q] 参考梯度 [dN/dxi, dN/deta, dN/dzeta]

    xi = q_pts(1, :);
    eta = q_pts(2, :);
    zeta = q_pts(3, :);
    
    % 节点顺序: 1(0,0,0), 2(1,0,0), 3(0,1,0), 4(0,0,1)
    % N1 = 1 - xi - eta - zeta
    % N2 = xi
    % N3 = eta
    % N4 = zeta
    
    n_pts = size(q_pts, 2);
    
    % 1. 函数值
    val = zeros(4, n_pts);
    val(1,:) = 1 - xi - eta - zeta;
    val(2,:) = xi;
    val(3,:) = eta;
    val(4,:) = zeta;
    
    % 2. 参考梯度 (常量)
    if nargout > 1
        grad = zeros(4, 3, n_pts);
        % N1 gradients: [-1, -1, -1]
        grad(1, 1, :) = -1; grad(1, 2, :) = -1; grad(1, 3, :) = -1;
        % N2 gradients: [1, 0, 0]
        grad(2, 1, :) = 1;
        % N3 gradients: [0, 1, 0]
        grad(3, 2, :) = 1;
        % N4 gradients: [0, 0, 1]
        grad(4, 3, :) = 1;
    end
end