function [J, detJ, iJt] = compute_jacobian_tet(P_elem)
% COMPUTE_JACOBIAN_TET 计算四面体单元的雅可比矩阵
% 输入:
%   P_elem - [3 x 4] 节点坐标矩阵 [x1 x2 x3 x4; y1...]
% 输出:
%   J      - [3 x 3] 雅可比矩阵 dx/dxi
%   detJ   - 标量，行列式
%   iJt    - [3 x 3] inv(J).' (即 J 的逆转置矩阵，用于梯度变换)
%
% 映射公式: x = x1 + (x2-x1)xi + (x3-x1)eta + (x4-x1)zeta

    % 向量化计算边向量
    x1 = P_elem(:, 1);
    x2 = P_elem(:, 2);
    x3 = P_elem(:, 3);
    x4 = P_elem(:, 4);
    
    d1 = x2 - x1; % dx/dxi
    d2 = x3 - x1; % dx/deta
    d3 = x4 - x1; % dx/dzeta
    
    J = [d1, d2, d3];
    
    % 手动计算 3x3 行列式 (比 det() 快)
    detJ = d1(1)*(d2(2)*d3(3) - d2(3)*d3(2)) - ...
           d1(2)*(d2(1)*d3(3) - d2(3)*d3(1)) + ...
           d1(3)*(d2(1)*d3(2) - d2(2)*d3(1));
           
    if nargout > 2
        % 计算 inv(J).' = inv(J^T)
        % M = inv(A) => M * det(A) = adj(A)
        % adj(A)^T = [cross(d2,d3), cross(d3,d1), cross(d1,d2)]
        
        c1 = cross(d2, d3);
        c2 = cross(d3, d1);
        c3 = cross(d1, d2);
        
        % iJt 的列就是 J^{-T} 的列，也就是梯度算子在参考系的基
        iJt = [c1, c2, c3] / detJ;
    end
end