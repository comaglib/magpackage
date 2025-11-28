function A_field = compute_biot_savart_A(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_A [鲁棒版] 计算磁矢量位
% 
% 改进: 对数项分母保护，防止数值爆炸。

    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    S1 = CoilSegments.P1; S2 = CoilSegments.P2; I_raw = CoilSegments.I;
    nSegsTotal = size(S1, 2);
    if isscalar(I_raw), I_raw = repmat(I_raw, 1, nSegsTotal); end
    
    SegVec = S2 - S1;
    L_sq = sum(SegVec.^2, 1);
    valid_mask = L_sq > 1e-20;
    
    S1 = S1(:, valid_mask); S2 = S2(:, valid_mask); I_vec = I_raw(:, valid_mask);
    SegVec = S2 - S1; L_seg = sqrt(sum(SegVec.^2, 1)); UnitDir = bsxfun(@rdivide, SegVec, L_seg);
    
    nSegs = size(S1, 2); nPoints = size(TargetPoints, 2);
    
    % 分块并行
    chunkSize = 5000; numChunks = ceil(nPoints / chunkSize);
    TargetBlocks = cell(numChunks, 1);
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1; idx_e = min(k*chunkSize, nPoints);
        TargetBlocks{k} = TargetPoints(:, idx_s:idx_e);
    end
    
    A_cell = cell(numChunks, 1);
    C_S1 = parallel.pool.Constant(S1); C_S2 = parallel.pool.Constant(S2);
    C_L = parallel.pool.Constant(L_seg); C_Dir = parallel.pool.Constant(UnitDir);
    C_I = parallel.pool.Constant(I_vec);
    
    parfor k = 1:numChunks
        P_targets = TargetBlocks{k};
        nLocal = size(P_targets, 2);
        A_local = zeros(3, nLocal);
        
        loc_S1 = C_S1.Value; loc_S2 = C_S2.Value; loc_L = C_L.Value; 
        loc_Dir = C_Dir.Value; loc_I = C_I.Value;
        
        EPS_DENOM = 1e-14;
        
        for i = 1:nSegs
            p1 = loc_S1(:, i); p2 = loc_S2(:, i); L = loc_L(i); curr = loc_I(i); dir = loc_Dir(:, i);
            
            R1_vec = bsxfun(@minus, P_targets, p1);
            R2_vec = bsxfun(@minus, P_targets, p2);
            R1 = sqrt(sum(R1_vec.^2, 1)); R2 = sqrt(sum(R2_vec.^2, 1));
            
            numer = R1 + R2 + L;
            denom = R1 + R2 - L;
            denom(denom < EPS_DENOM) = EPS_DENOM; % 保护
            
            val = (ConstFactor * curr) * log(numer ./ denom);
            A_local = A_local + bsxfun(@times, dir, val);
        end
        A_cell{k} = A_local;
    end
    A_field = [A_cell{:}];
end