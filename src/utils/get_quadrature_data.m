function [q_pts, q_w] = get_quadrature_data(elemType, order)
% GET_QUADRATURE_DATA 获取单元的高斯积分点和权重
% 输入:
%   elemType - 'tet', 'tri', 'hex' 等
%   order    - 积分精度阶数 (近似)
% 输出:
%   q_pts    - [3 x N_q] 参考坐标 (xi, eta, zeta)
%   q_w      - [1 x N_q] 权重

    if strcmpi(elemType, 'tet')
        % 四面体积分规则 (Keast, 1986)
        if order <= 1
            % 1点积分 (a=0.25, b=0.25, c=0.25, w=1/6)
            % 注意：标准四面体体积为 1/6，这里的权重通常归一化到体积 1 或 1/6
            % 我们约定：权重之和 = 1/6 (参考体积)
            q_pts = [0.25; 0.25; 0.25];
            q_w   = 1.0/6.0;
        else
            % 4点积分 (精度阶数 2)
            a = 0.58541020; b = 0.13819660;
            q_pts = [a b b b;
                     b a b b;
                     b b a b]; % [xi; eta; zeta]
            % 最后一个坐标由 1-xi-eta-zeta 隐式定义用于形状函数，
            % 但这里我们直接给出前三个独立坐标
            
            w = 1.0/24.0; % 4个点的权重相同，和为 1/6
            q_w = repmat(w, 1, 4);
        end
    else
        error('Unsupported element type: %s', elemType);
    end
end