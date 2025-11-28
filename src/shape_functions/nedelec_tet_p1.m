function [val, curl_val] = nedelec_tet_p1(q_pts)
% NEDELEC_TET_P1 计算 Nedelec (Edge) P1 基函数及其参考旋度
% 参考: Whitney Edge Elements
% Edge定义需与 Mesh 类一致: 
% 1:1-2, 2:1-3, 3:1-4, 4:2-3, 5:2-4, 6:3-4
%
% 输入:
%   q_pts - [3 x N_q] 参考坐标
% 输出:
%   val      - [6 x 3 x N_q] 参考矢量值 W_ref
%   curl_val - [6 x 3 x N_q] 参考旋度值 curl(W_ref)

    % 首先获取 Lagrange 梯度 (用于构造 Whitney 形式)
    % W_ij = lambda_i * grad(lambda_j) - lambda_j * grad(lambda_i)
    [L_val, L_grad] = lagrange_tet_p1(q_pts);
    
    n_pts = size(q_pts, 2);
    edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    
    val = zeros(6, 3, n_pts);
    curl_val = zeros(6, 3, n_pts);
    
    % 参考梯度是常数，可以提取出来简化计算
    % G1=[-1,-1,-1], G2=[1,0,0], G3=[0,1,0], G4=[0,0,1]
    % 预计算 curl(W_ref) = 2 * cross(grad_i, grad_j)
    % 由于是常数梯度，参考旋度在整个单元内是常数
    
    for e = 1:6
        i = edge_defs(e, 1);
        j = edge_defs(e, 2);
        
        % --- 计算 W_ref ---
        % W = Li * G_j - Lj * G_i
        % 这是一个 [3 x N_q] 的向量
        
        G_i = L_grad(i, :, 1); % [1x3]
        G_j = L_grad(j, :, 1);
        
        % 向量化计算
        % val(e, c, p) = Li(p)*Gj(c) - Lj(p)*Gi(c)
        term1 = bsxfun(@times, reshape(G_j, 1, 3), reshape(L_val(i,:), n_pts, 1));
        term2 = bsxfun(@times, reshape(G_i, 1, 3), reshape(L_val(j,:), n_pts, 1));
        
        val(e, :, :) = permute(term1 - term2, [3, 2, 1]);
        
        % --- 计算 Curl(W_ref) ---
        % Curl = 2 * cross(grad_i, grad_j)
        c_vec = 2 * cross(G_i, G_j); % [1x3] constant
        
        curl_val(e, :, :) = repmat(reshape(c_vec, 1, 3, 1), 1, 1, n_pts);
    end
end