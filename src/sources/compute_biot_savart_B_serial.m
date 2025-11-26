function B_field = compute_biot_savart_B_serial(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_B_SERIAL 串行版 (用于 parfor 内部调用)
% 
% 修复: 增加对标量 CoilSegments.I 的支持，防止索引越界。

    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    % 1. 预处理线圈
    S1 = CoilSegments.P1; 
    S2 = CoilSegments.P2;
    I_vec = CoilSegments.I;
    
    nSegsTotal = size(S1, 2);
    
    % [关键修复] 兼容标量电流输入
    if isscalar(I_vec)
        I_vec = repmat(I_vec, 1, nSegsTotal);
    end
    
    SegVec = S2 - S1;
    L_sq = sum(SegVec.^2, 1);
    valid_idx = L_sq > 1e-24;
    
    S1 = S1(:, valid_idx);
    S2 = S2(:, valid_idx);
    I_vec = I_vec(:, valid_idx);
    
    nSegs = size(S1, 2);
    nPoints = size(TargetPoints, 2);
    
    B_field = zeros(3, nPoints);
    
    % 2. 循环线段
    for i = 1:nSegs
        p1 = S1(:, i);
        p2 = S2(:, i);
        curr = I_vec(i);
        
        R1 = TargetPoints - p1; 
        R2 = TargetPoints - p2; 
        
        Cx = R1(2,:).*R2(3,:) - R1(3,:).*R2(2,:);
        Cy = R1(3,:).*R2(1,:) - R1(1,:).*R2(3,:);
        Cz = R1(1,:).*R2(2,:) - R1(2,:).*R2(1,:);
        
        cross_sq = Cx.^2 + Cy.^2 + Cz.^2;
        cross_sq(cross_sq < 1e-20) = 1e-20; 
        
        d1 = sqrt(sum(R1.^2, 1));
        d2 = sqrt(sum(R2.^2, 1));
        
        u = p2 - p1;
        dot1 = u(1)*R1(1,:) + u(2)*R1(2,:) + u(3)*R1(3,:);
        dot2 = u(1)*R2(1,:) + u(2)*R2(2,:) + u(3)*R2(3,:);
        
        scalar = dot1 ./ d1 - dot2 ./ d2;
        coeff = (ConstFactor * curr) * scalar ./ cross_sq;
        
        B_field(1,:) = B_field(1,:) + Cx .* coeff;
        B_field(2,:) = B_field(2,:) + Cy .* coeff;
        B_field(3,:) = B_field(3,:) + Cz .* coeff;
    end
end