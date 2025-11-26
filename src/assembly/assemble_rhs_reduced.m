function b = assemble_rhs_reduced(Model, CoilSegments, Nu_Override)
% ASSEMBLE_RHS_REDUCED 组装 A-V 表述下的右端载荷向量 (HPC重构版 - 修复)
% 
% 修正:
% 1. 修复了 parfor 循环内部变量名拼写错误 (local_T2E -> loc_T2E)。
% 2. 保持了 HPC 优化结构。

    % fprintf('正在组装右端项 (RHS Reduced)...\n');
    % t_start = tic;
    
    Mesh = Model.Mesh;
    numEdges = size(Mesh.Edges, 2);
    numElems = size(Mesh.T, 2);
    
    mu0 = 4*pi*1e-7;
    nu0 = 1/mu0;
    
    % 1. 确定 Nu 向量
    if nargin > 2 && ~isempty(Nu_Override)
        Nu_vec = Nu_Override;
    else
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        Nu_vec = zeros(numElems, 1);
        
        % 简单的循环 (MatLib 很小)
        LibNu = zeros(length(MatLib), 1);
        for i = 1:length(MatLib)
            m_r = MatLib(i).Mu_r;
            if isempty(m_r), m_r=1.0; end
            LibNu(i) = 1/(m_r * mu0);
        end
        Nu_vec = LibNu(MatMap);
    end
    
    % 2. 筛选有效单元 (nu != nu0)
    % 允许一点浮点误差
    tol = nu0 * 1e-4;
    target_elems = find(abs(Nu_vec - nu0) > tol);
    
    numTarget = length(target_elems);
    % fprintf('  - 磁化源积分区域: %d / %d 单元\n', numTarget, numElems);
    
    if numTarget == 0
        b = sparse(numEdges, 1);
        return;
    end
    
    % 3. 数据打包
    chunkSize = 2000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        g_idxs = target_elems(s_idx:e_idx);
        
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.T2E = T2E_raw(:, g_idxs);
        ElementBlocks{k}.Signs = Signs_raw(:, g_idxs);
        ElementBlocks{k}.Nu = Nu_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    b_cell = cell(numChunks, 1);
    
    C_P = parallel.pool.Constant(Mesh.P);
    C_Coil = parallel.pool.Constant(CoilSegments);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        % 统一变量名为 loc_*
        loc_T = Block.T; 
        loc_T2E = Block.T2E; 
        loc_Signs = Block.Signs;
        loc_Nu = Block.Nu;
        
        loc_P_all = C_P.Value;
        loc_Coil = C_Coil.Value;
        
        idx_store = zeros(6 * count, 1);
        val_store = zeros(6 * count, 1);
        ptr = 0;
        
        % 批量计算中心点
        Centers = zeros(3, count);
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            Centers(:, e) = mean(pts, 2); 
        end
        
        % 调用串行版 Biot-Savart
        Bs_all = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            nu_val = loc_Nu(e);
            Bs = Bs_all(:, e); 
            
            % 几何计算
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            J_mat = [v21, v31, v41]; detJ = det(J_mat); invJ_T = inv(J_mat)';
            G_phy = invJ_T * [-1 1 0 0; -1 0 1 0; -1 0 0 1];
            
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            Curl_N = zeros(3, 6);
            for i = 1:6
                na = edge_pairs(i,1); nb = edge_pairs(i,2);
                Curl_N(:, i) = 2 * cross(G_phy(:, na), G_phy(:, nb));
            end
            Vol = abs(detJ) / 6.0;
            
            % RHS = Vol * (Curl N)' * (nu0 - nu) * Bs
            M_contrast = (nu0 - nu_val) * Bs;
            local_b = Vol * (Curl_N' * M_contrast);
            
            local_b = local_b .* loc_Signs(:, e);
            
            range = ptr+1 : ptr+6;
            idx_store(range) = loc_T2E(:, e);
            val_store(range) = local_b;
            ptr = ptr + 6;
        end
        
        b_cell{k} = sparse(idx_store, 1, val_store, numEdges, 1);
    end
    
    b = sum([b_cell{:}], 2); 
    
    % t_end = toc(t_start);
    % fprintf('右端项组装完成。耗时: %.4f 秒\n', t_end);
end