function A_field = compute_biot_savart_A(CoilSegments, TargetPoints)
% COMPUTE_BIOT_SAVART_A 计算线圈在目标点产生的磁矢量位 A (HPC优化版)
%
% 优化:
% 1. 对 TargetPoints 进行预分块 (Pre-chunking)，消除广播通信开销。
% 2. 保留了对零长度线段的 NaN 防护。

    % fprintf('正在计算源磁场 (Biot-Savart A)...\n');
    % t_start = tic;
    
    mu0 = 4*pi*1e-7;
    ConstFactor = mu0 / (4*pi);
    
    nSegs = size(CoilSegments.P1, 2);
    nPoints = size(TargetPoints, 2);
    
    S1 = CoilSegments.P1; 
    S2 = CoilSegments.P2;
    I_vec = CoilSegments.I;
    
    % 预计算线段矢量与长度
    SegVec = S2 - S1; 
    L_seg = sqrt(sum(SegVec.^2, 1));
    
    % [NaN防护] 找出长度极短的无效线段
    valid_mask = L_seg > 1e-12;
    
    UnitDir = zeros(3, nSegs);
    UnitDir(:, valid_mask) = SegVec(:, valid_mask) ./ L_seg(valid_mask);
    
    % --------------------------------------------------------
    % 并行策略优化: 数据预分块 (Pre-chunking)
    % --------------------------------------------------------
    chunkSize = 10000; % 任务包大小
    numChunks = ceil(nPoints / chunkSize);
    
    % [核心优化] 将大矩阵 TargetPoints 切分为 Cell 数组
    % 这避免了在 parfor 中对 TargetPoints 进行动态切片导致的整阵广播
    TargetBlocks = cell(numChunks, 1);
    for k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, nPoints);
        TargetBlocks{k} = TargetPoints(:, idx_start:idx_end);
    end
    
    % 结果容器
    A_cell = cell(numChunks, 1);
    
    % 广播变量封装 (线圈数据较小且需全量访问，使用 Constant 优化)
    C_S1 = parallel.pool.Constant(S1);
    C_S2 = parallel.pool.Constant(S2);
    C_UnitDir = parallel.pool.Constant(UnitDir);
    C_I = parallel.pool.Constant(I_vec);
    
    parfor k = 1:numChunks
        % 1. 获取当前任务包 (无通信开销，仅传输必要的块)
        P_targets = TargetBlocks{k};
        nLocal = size(P_targets, 2);
        
        A_local = zeros(3, nLocal);
        
        % 获取线圈常数
        loc_S1 = C_S1.Value;
        loc_S2 = C_S2.Value;
        loc_Dir = C_UnitDir.Value;
        loc_I = C_I.Value;
        
        % 2. 线圈线段循环
        for i = 1:nSegs
            % 跳过无效线段 (方向全0)
            dir = loc_Dir(:, i);
            if all(dir == 0), continue; end
            
            p1 = loc_S1(:, i);
            p2 = loc_S2(:, i);
            curr = loc_I(i);
            
            % 矢量化计算 R1, R2
            R1_vec = P_targets - p1; 
            R2_vec = P_targets - p2; 
            
            R1 = sqrt(sum(R1_vec.^2, 1));
            R2 = sqrt(sum(R2_vec.^2, 1));
            
            L = norm(p2 - p1);
            
            % 解析公式 (防止除零)
            numer = R1 + R2 + L;
            denom = R1 + R2 - L;
            denom(denom < 1e-12) = 1e-12; 
            
            log_term = log(numer ./ denom);
            
            val = (ConstFactor * curr) * log_term; 
            A_local = A_local + dir * val; 
        end
        
        A_cell{k} = A_local;
    end
    
    % 3. 合并结果
    A_field = [A_cell{:}];
    
    % t_end = toc(t_start);
    % fprintf('  - 计算完成. 耗时: %.4f 秒 (Pts: %d, Segs: %d)\n', t_end, nPoints, nSegs);
end