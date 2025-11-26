function A_field = compute_biot_savart_A(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_A 计算线圈在目标点产生的磁矢量位 A (HPC优化版)
%
% 修复: 增加对标量 CoilSegments.I 的支持。

    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    nSegsTotal = size(CoilSegments.P1, 2);
    nPoints = size(TargetPoints, 2);
    
    S1_raw = CoilSegments.P1; 
    S2_raw = CoilSegments.P2;
    I_raw = CoilSegments.I;
    
    % [关键修复] 兼容标量电流
    if isscalar(I_raw)
        I_raw = repmat(I_raw, 1, nSegsTotal);
    end
    
    SegVec_raw = S2_raw - S1_raw;
    L_sq = sum(SegVec_raw.^2, 1);
    
    % 找出长度极短的无效线段
    valid_mask = L_sq > 1e-24;
    
    % 压缩有效数据
    S1 = S1_raw(:, valid_mask);
    S2 = S2_raw(:, valid_mask);
    I_vec = I_raw(:, valid_mask);
    
    % 重新计算方向
    SegVec = S2 - S1;
    L_seg = sqrt(sum(SegVec.^2, 1));
    UnitDir = SegVec ./ L_seg;
    
    nSegs = size(S1, 2);
    
    % Pre-chunking
    chunkSize = 10000;
    numChunks = ceil(nPoints / chunkSize);
    TargetBlocks = cell(numChunks, 1);
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, nPoints);
        TargetBlocks{k} = TargetPoints(:, idx_s:idx_e);
    end
    
    A_cell = cell(numChunks, 1);
    
    C_S1 = parallel.pool.Constant(S1);
    C_S2 = parallel.pool.Constant(S2);
    C_UnitDir = parallel.pool.Constant(UnitDir);
    C_I = parallel.pool.Constant(I_vec);
    
    parfor k = 1:numChunks
        P_targets = TargetBlocks{k};
        nLocal = size(P_targets, 2);
        A_local = zeros(3, nLocal);
        
        loc_S1 = C_S1.Value;
        loc_S2 = C_S2.Value;
        loc_Dir = C_UnitDir.Value;
        loc_I = C_I.Value;
        
        for i = 1:nSegs
            p1 = loc_S1(:, i);
            p2 = loc_S2(:, i);
            curr = loc_I(i);
            dir = loc_Dir(:, i);
            
            R1_vec = P_targets - p1; 
            R2_vec = P_targets - p2; 
            
            R1 = sqrt(sum(R1_vec.^2, 1));
            R2 = sqrt(sum(R2_vec.^2, 1));
            
            L = norm(p2 - p1);
            
            numer = R1 + R2 + L;
            denom = R1 + R2 - L;
            denom(denom < 1e-12) = 1e-12; 
            
            log_term = log(numer ./ denom);
            
            val = (ConstFactor * curr) * log_term; 
            A_local = A_local + dir * val; 
        end
        A_cell{k} = A_local;
    end
    
    A_field = [A_cell{:}];
end