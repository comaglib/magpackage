classdef SpectralTimeUtils
    % SPECTRALTIMEUTILS 谱时间离散工具类 (v2.0 - Fixed Node Ordering)
    % 
    % 修复日志:
    %   1. [Critical] GLL 节点生成后强制升序排序 (sort)。
    %      修复了之前生成降序节点导致 Prediction 步出现负时间步长 (dt < 0) 
    %      以及积分矩阵 Q 构建错误的问题。
    
    methods (Static)
        function [nodes, weights] = gll(P)
            % GLL 计算 P 阶 (P+1 个点) 的 GLL 节点和权重
            % 节点范围: [-1, 1]
            
            N = P; 
            % 初始猜测 (Chebyshev 节点) [1, ..., -1]
            nodes = cos(pi * (0:N) / N)';
            
            % 牛顿迭代
            P_N = zeros(N+1, 1);
            maxIter = 100; tol = 1e-12;
            
            for k = 1:length(nodes)
                x = nodes(k);
                for iter = 1:maxIter
                    [L, dL] = SpectralTimeUtils.legendre_poly(x, N);
                    val = (1-x^2) * dL;
                    if abs(val) < tol, break; end
                    
                    d2L = (2*x*dL - N*(N+1)*L) / (1-x^2);
                    d_val = -2*x*dL + (1-x^2)*d2L;
                    x = x - val / d_val;
                end
                nodes(k) = x;
            end
            
            % [CRITICAL FIX] 强制升序排列 [-1, ..., 1]
            % 确保 dt = nodes(i) - nodes(i-1) > 0
            nodes = sort(nodes, 'ascend');
            
            % 计算权重 (基于排序后的节点)
            weights = zeros(N+1, 1);
            for k = 1:length(nodes)
                [L, ~] = SpectralTimeUtils.legendre_poly(nodes(k), N);
                weights(k) = 2 / (N * (N+1) * L^2);
            end
        end
        
        function Q = integration_matrix(nodes)
            % INTEGRATION_MATRIX 构建谱积分矩阵 Q
            % Q_ij = Integral_{-1}^{tau_i} L_j(tau) dtau
            
            N = length(nodes) - 1;
            Q = zeros(N+1, N+1);
            
            % 使用高阶高斯积分
            [q_pts, q_w] = SpectralTimeUtils.gl_quad(N + 2); 
            
            for i = 1:N+1
                b = nodes(i);
                a = -1; 
                % 既然节点已排序，第一点一定是 -1，积分区间为 0，跳过
                if abs(b - a) < 1e-14
                    continue; 
                end
                
                % 映射积分区间 [-1, 1] -> [a, b]
                mapped_pts = 0.5 * ((b-a)*q_pts + (a+b));
                mapped_w   = 0.5 * (b-a) * q_w;
                
                for j = 1:N+1
                    L_vals = SpectralTimeUtils.lagrange_basis(mapped_pts, nodes, j);
                    Q(i, j) = sum(mapped_w .* L_vals);
                end
            end
        end
        
        function [S, Lambda, S_inv, pair_map] = decompose_for_sdc(Q, dt)
            % DECOMPOSE_FOR_SDC 执行特征分解并识别共轭对
            
            [S, D] = eig(Q);
            lam = diag(D);
            S_inv = inv(S);
            
            M = length(lam);
            pair_map = repmat(struct('Type', '', 'ConjIdx', 0), M, 1);
            
            tol = 1e-10;
            visited = false(M, 1);
            
            for i = 1:M
                if visited(i), continue; end
                
                val = lam(i);
                if abs(imag(val)) < tol
                    pair_map(i).Type = 'Real';
                    pair_map(i).ConjIdx = i;
                    visited(i) = true;
                else
                    dist = abs(lam - conj(val));
                    [~, idx] = min(dist);
                    
                    if idx == i
                        pair_map(i).Type = 'Real'; 
                    else
                        if imag(val) > 0
                            idx_A = i; idx_B = idx;
                        else
                            idx_A = idx; idx_B = i;
                        end
                        pair_map(idx_A).Type = 'ComplexA'; 
                        pair_map(idx_A).ConjIdx = idx_B;
                        pair_map(idx_B).Type = 'ComplexB'; 
                        pair_map(idx_B).ConjIdx = idx_A;
                        visited(idx_A) = true; visited(idx_B) = true;
                    end
                end
            end
            Lambda = D;
        end
        
        % --- 辅助函数 ---
        function [L, dL] = legendre_poly(x, N)
            if N == 0
                L = 1; dL = 0;
            elseif N == 1
                L = x; dL = 1;
            else
                L_prev2 = 1; L_prev1 = x;
                dL_prev2 = 0; dL_prev1 = 1;
                for k = 2:N
                    L = ( (2*k-1)*x*L_prev1 - (k-1)*L_prev2 ) / k;
                    dL = dL_prev2 + (2*k-1)*L_prev1;
                    L_prev2 = L_prev1; L_prev1 = L;
                    dL_prev2 = dL_prev1; dL_prev1 = dL;
                end
            end
        end
        
        function vals = lagrange_basis(eval_pts, nodes, j)
            N = length(nodes);
            vals = ones(size(eval_pts));
            xj = nodes(j);
            for k = 1:N
                if k == j, continue; end
                xk = nodes(k);
                vals = vals .* (eval_pts - xk) / (xj - xk);
            end
        end
        
        function [x, w] = gl_quad(N)
            beta = .5 ./ sqrt(1-(2*(1:N-1)).^(-2));
            T = diag(beta,1) + diag(beta,-1);
            [V, D] = eig(T);
            x = diag(D); 
            [x, i] = sort(x);
            w = 2 * V(1,i).^2;
            w = w(:);
        end
    end
end