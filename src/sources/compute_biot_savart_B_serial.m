function B_field = compute_biot_savart_B_serial(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_B_SERIAL 计算源磁感应强度 B (Serial)
%
% 输入:
%   CoilSegments: 结构体, 包含 P1, P2 (3xN), I (1xN or scalar)
%   TargetPoints: 3xM 目标点
% 输出:
%   B_field: 3xM 磁感应强度 (Tesla)

    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    S1 = CoilSegments.P1; 
    S2 = CoilSegments.P2; 
    I_raw = CoilSegments.I;
    
    nSegs = size(S1, 2);
    nPoints = size(TargetPoints, 2);
    
    % 鲁棒性处理: 确保 I 是 1xN
    if isscalar(I_raw)
        I_vec = repmat(I_raw, 1, nSegs);
    elseif size(I_raw, 1) > 1 && size(I_raw, 2) == nSegs
        % 如果 I 是 TimeSteps x Segments，取第一行
        I_vec = I_raw(1, :);
    else
        I_vec = I_raw;
    end
    
    % 预计算线段向量
    SegVec = S2 - S1; 
    L_sq = sum(SegVec.^2, 1);
    
    % 过滤极短线段
    valid_mask = L_sq > 1e-20;
    S1 = S1(:, valid_mask);
    S2 = S2(:, valid_mask);
    I_vec = I_vec(valid_mask);
    SegVec = SegVec(:, valid_mask);
    L_seg = sqrt(L_sq(valid_mask));
    UnitDir = bsxfun(@rdivide, SegVec, L_seg);
    
    nSegsValid = size(S1, 2);
    B_field = zeros(3, nPoints);
    
    for i = 1:nSegsValid
        p1 = S1(:, i);
        % p2 = S2(:, i);
        u = UnitDir(:, i);
        curr = I_vec(i);
        
        % 向量化计算所有目标点对当前线段的贡献
        % R1 = P - P1
        R1_vec = bsxfun(@minus, TargetPoints, p1);
        R1_sq = sum(R1_vec.^2, 1);
        R1 = sqrt(R1_sq);
        
        % 投影 s = R1 . u
        s = sum(bsxfun(@times, R1_vec, u), 1);
        
        % 垂直距离平方 d^2 = R1^2 - s^2
        d_sq = R1_sq - s.^2;
        d_sq(d_sq < 1e-20) = 1e-20; % 保护
        
        % 积分因子
        % Integral = (s / (d^2 * R1)) - ((s - L) / (d^2 * R2))
        % 这里 R2^2 = d^2 + (s-L)^2
        L = L_seg(i);
        s_minus_L = s - L;
        R2 = sqrt(d_sq + s_minus_L.^2);
        
        term1 = s ./ R1;
        term2 = s_minus_L ./ R2;
        factor = (term1 - term2) ./ d_sq;
        
        % 叉乘方向 (u x R1)
        % u = [ux; uy; uz], R1 = [rx; ry; rz]
        % cross_x = uy*rz - uz*ry
        cross_x = u(2)*R1_vec(3,:) - u(3)*R1_vec(2,:);
        cross_y = u(3)*R1_vec(1,:) - u(1)*R1_vec(3,:);
        cross_z = u(1)*R1_vec(2,:) - u(2)*R1_vec(1,:);
        
        % 累加 B
        scalar_k = (ConstFactor * curr) * factor;
        B_field(1,:) = B_field(1,:) + scalar_k .* cross_x;
        B_field(2,:) = B_field(2,:) + scalar_k .* cross_y;
        B_field(3,:) = B_field(3,:) + scalar_k .* cross_z;
    end
end