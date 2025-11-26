function b = assemble_rhs_reduced(Model, CoilSegments)
% ASSEMBLE_RHS_REDUCED 组装 A-V 表述下的右端载荷向量 (HPC修正版 v2)
% 修复: 消除 parfor 中的所有广播变量警告
% 物理修正: 使用 (nu0 - nu) * Bs 作为源项，正确模拟磁化效应。
% HPC修正:  使用预分块 (Pre-chunking) 消除 parfor 广播开销。

    fprintf('正在组装右端项 (RHS Reduced)...\n');
    t_start = tic;
    
    Mesh = Model.Mesh;
    numEdges = size(Mesh.Edges, 2);
    numElems = size(Mesh.T, 2);
    
    % 1. 筛选材料
    mu0 = 4*pi*1e-7;
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    Nu_vec = zeros(1, numElems);
    IsMagnetic = false(1, numElems);
    
    for i = 1:numElems
        mat_id = MatMap(i);
        if mat_id > length(MatLib), continue; end
        
        mu_r = MatLib(mat_id).Mu_r;
        if isempty(mu_r), mu_r = 1.0; end % 防御性编程
        
        if abs(mu_r - 1.0) > 1e-4
            IsMagnetic(i) = true;
        end
        Nu_vec(i) = 1 / (mu_r * mu0);
    end
    
    target_elems = find(IsMagnetic);
    numTarget = length(target_elems);
    fprintf('  - 磁性区域单元数: %d / %d\n', numTarget, numElems);
    
    if numTarget == 0
        b = sparse(numEdges, 1);
        return;
    end
    
    % 2. 严格的数据打包 (Chunking)
    chunkSize = 2000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    % 提取原始数组 (避免在循环中引用 Mesh 结构体)
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        
        % 获取当前块的全局单元索引
        g_idxs = target_elems(s_idx:e_idx);
        
        % 只拷贝需要的数据片段到 Cell
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.T2E = T2E_raw(:, g_idxs);
        ElementBlocks{k}.Signs = Signs_raw(:, g_idxs);
        ElementBlocks{k}.Nu = Nu_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    b_cell = cell(numChunks, 1);
    
    % 广播常量
    C_P = parallel.pool.Constant(Mesh.P);
    C_Coil = parallel.pool.Constant(CoilSegments);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k}; % 获取本地副本
        
        count = Block.Count;
        local_T = Block.T;
        local_T2E = Block.T2E;
        local_Signs = Block.Signs;
        local_Nu = Block.Nu;
        
        local_P_all = C_P.Value;
        local_Coil = C_Coil.Value;
        
        idx_store = zeros(6 * count, 1);
        val_store = zeros(6 * count, 1);
        ptr = 0;
        
        % 批量计算中心点
        Centers = zeros(3, count);
        for e = 1:count
            nodes = local_T(:, e);
            pts = local_P_all(:, nodes);
            Centers(:, e) = mean(pts, 2); 
        end
        
        Bs_all = compute_biot_savart_B(local_Coil, Centers);
        
        for e = 1:count
            nodes = local_T(:, e);
            pts = local_P_all(:, nodes);
            nu_val = local_Nu(e);
            Bs = Bs_all(:, e); 
            
            % 几何计算 (Curl N)
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            J_mat = [v21, v31, v41];
            detJ = det(J_mat);
            invJ_T = inv(J_mat)';
            
            G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
            G_phy = invJ_T * G_ref;
            
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            Curl_N = zeros(3, 6);
            for i = 1:6
                na = edge_pairs(i,1); nb = edge_pairs(i,2);
                Curl_N(:, i) = 2 * cross(G_phy(:, na), G_phy(:, nb));
            end
            
            Vol = abs(detJ) / 6.0;
            
            % RHS = Vol * (Curl N)' * (nu0 - nu) * Bs
            M_contrast = (1/(4*pi*1e-7) - nu_val) * Bs;
            local_b = Vol * (Curl_N' * M_contrast);
            
            local_b = local_b .* local_Signs(:, e);
            
            range = ptr+1 : ptr+6;
            idx_store(range) = local_T2E(:, e);
            val_store(range) = local_b;
            ptr = ptr + 6;
        end
        
        b_cell{k} = sparse(idx_store, 1, val_store, numEdges, 1);
    end
    
    b = sum([b_cell{:}], 2); 
    
    t_end = toc(t_start);
    fprintf('右端项组装完成。耗时: %.4f 秒\n', t_end);
end