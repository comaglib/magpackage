function B_field = compute_biot_savart_B(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_B 计算线圈在目标点产生的磁通密度 B (HPC修正版 v2)
% 
% 修正: 彻底消除 TargetPoints 在 parfor 中的广播开销

    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    % 1. 预处理线圈
    S1_raw = CoilSegments.P1; 
    S2_raw = CoilSegments.P2;
    I_raw = CoilSegments.I;
    
    SegVec_raw = S2_raw - S1_raw;
    L_sq = sum(SegVec_raw.^2, 1);
    valid_idx = L_sq > 1e-24;
    
    S1 = S1_raw(:, valid_idx);
    S2 = S2_raw(:, valid_idx);
    I_vec = I_raw(:, valid_idx);
    
    nSegs = size(S1, 2);
    nPoints = size(TargetPoints, 2);
    
    % 2. 切片 (Pre-chunking)
    chunkSize = 10000;
    numChunks = ceil(nPoints / chunkSize);
    TargetBlocks = cell(numChunks, 1);
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, nPoints);
        % 将切片存入 Cell，切断与原矩阵的联系
        TargetBlocks{k} = TargetPoints(:, idx_s:idx_e);
    end
    
    B_cell = cell(numChunks, 1);
    
    % 3. 广播常量
    C_S1 = parallel.pool.Constant(S1);
    C_S2 = parallel.pool.Constant(S2);
    C_I = parallel.pool.Constant(I_vec);
    
    parfor k = 1:numChunks
        % 仅访问 Cell 元素，MATLAB 能够识别
        P_targets = TargetBlocks{k};
        nLocal = size(P_targets, 2);
        B_local = zeros(3, nLocal);
        
        loc_S1 = C_S1.Value;
        loc_S2 = C_S2.Value;
        loc_I = C_I.Value;
        
        for i = 1:nSegs
            p1 = loc_S1(:, i);
            p2 = loc_S2(:, i);
            curr = loc_I(i);
            
            R1 = P_targets - p1; 
            R2 = P_targets - p2; 
            
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
            
            B_local(1,:) = B_local(1,:) + Cx .* coeff;
            B_local(2,:) = B_local(2,:) + Cy .* coeff;
            B_local(3,:) = B_local(3,:) + Cz .* coeff;
        end
        B_cell{k} = B_local;
    end
    
    B_field = [B_cell{:}];
end