function G = assemble_scalar_laplacian(Model, DofData)
% ASSEMBLE_SCALAR_LAPLACIAN 组装标量电导矩阵 (支持压缩 DoF)
%
% 输入: Model, DofData
% 输出: G (TotalDoFs x TotalDoFs) 
%       注意: 为了方便叠加，我们直接输出与系统同维度的稀疏矩阵
%       或者输出 (N_active x N_active) 的小矩阵，由调用者拼接。
%       *为了清晰，这里输出 (N_active x N_active)* 的局部矩阵。

    % fprintf('正在组装标量电导矩阵 (G_sigma)...\n');
    
    Mesh = Model.Mesh;
    % T_raw = Mesh.T; % Unused directly
    
    % 1. 识别导电单元 (同前)
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    numElems = size(Mesh.T, 2);
    Sigma_vec = zeros(1, numElems);
    IsConductive = false(1, numElems);
    
    for i = 1:numElems
        mat_id = MatMap(i);
        if mat_id > length(MatLib), continue; end
        sig = MatLib(mat_id).Sigma;
        if sig > 1e-12
            Sigma_vec(i) = sig;
            IsConductive(i) = true;
        end
    end
    
    target_elems = find(IsConductive);
    numTarget = length(target_elems);
    
    numActiveV = DofData.NumActiveNodes;
    if numTarget == 0
        G = sparse(numActiveV, numActiveV);
        return;
    end
    
    % 2. 预分块
    chunkSize = 5000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw = Mesh.T;
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, numTarget);
        g_idxs = target_elems(idx_s:idx_e);
        
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.Sigma = Sigma_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    % 广播
    C_P = parallel.pool.Constant(Mesh.P);
    
    % Node2DoF 映射 (本地化 V 索引: GlobalNodeID -> 1..numActive)
    % DofData.Node2DoF 存的是 Global DoF (e.g. numEdges+1...)
    % 我们需要将其转换为 1-based 的局部 V 索引
    GlobalDoFMap = DofData.Node2DoF;
    V_offset = DofData.NumEdges;
    
    % LocalVMap: NodeID -> 1..numActive
    % 只有 active nodes 有值
    LocalVMap = zeros(size(GlobalDoFMap));
    mask = GlobalDoFMap > 0;
    LocalVMap(mask) = GlobalDoFMap(mask) - V_offset;
    
    C_VMap = parallel.pool.Constant(LocalVMap);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        local_T = Block.T;
        local_Sigma = Block.Sigma;
        local_P_all = C_P.Value;
        local_VMap = C_VMap.Value;
        
        i_block = zeros(16 * count, 1);
        j_block = zeros(16 * count, 1);
        v_block = zeros(16 * count, 1);
        ptr = 0;
        
        for e = 1:count
            nodes = local_T(:, e);
            
            % 获取 1-based V 索引
            v_indices = local_VMap(nodes);
            
            % 安全检查 (虽然理论上导电单元的节点必然是 Active)
            if any(v_indices == 0), continue; end
            
            p_elem = local_P_all(:, nodes);
            sig_val = local_Sigma(e);
            
            Ke = Ke_grad_grad(p_elem, sig_val);
            
            R_idx = repmat(v_indices, 1, 4);
            C_idx = repmat(v_indices', 4, 1);
            
            range = ptr+1 : ptr+16;
            i_block(range) = R_idx(:);
            j_block(range) = C_idx(:);
            v_block(range) = Ke(:);
            ptr = ptr + 16;
        end
        
        % 截断
        if ptr < 16*count
            i_block = i_block(1:ptr);
            j_block = j_block(1:ptr);
            v_block = v_block(1:ptr);
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    G = sparse(I, J, V, numActiveV, numActiveV);
end