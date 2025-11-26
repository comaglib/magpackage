function C = assemble_coupling_matrix(Model, mode, DofData)
% ASSEMBLE_COUPLING_MATRIX 并行组装 (支持压缩 DoF)
% 
% 输入: Model, mode, DofData
%       mode - 'Gauge' (默认): 系数 = 1.0, 用于 div A = 0 约束
%              'Physical'    : 系数 = Sigma, 用于 J = sigma (-dA/dt - grad V)
% 输出: C (NumEdges x NumActiveNodes)

    if nargin < 2, mode = 'Gauge'; end
    % 如果 mode='Gauge'，通常不需要压缩 V，或者全部 V 都 Active
    % 这里主要针对 'Physical' 模式
    
    fprintf('正在组装耦合矩阵 (Mode: %s)...\n', mode);
    
    Mesh = Model.Mesh;
    C_P = parallel.pool.Constant(Mesh.P);
    
    numElems = size(Mesh.T, 2);
    numEdges = size(Mesh.Edges, 2);
    
    if nargin < 3
        % 兼容旧调用: 假设全活跃
        numActiveV = size(Mesh.P, 2);
        LocalVMap = 1:numActiveV;
    else
        numActiveV = DofData.NumActiveNodes;
        % 构建 NodeID -> 1..NumActive
        GlobalDoFMap = DofData.Node2DoF;
        V_offset = DofData.NumEdges;
        LocalVMap = zeros(size(GlobalDoFMap));
        mask = GlobalDoFMap > 0;
        LocalVMap(mask) = GlobalDoFMap(mask) - V_offset;
    end
    
    % 1. 准备系数
    Coeff_vec = zeros(1, numElems);
    if strcmpi(mode, 'Physical')
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        for i = 1:numElems
            mat_id = MatMap(i);
            Coeff_vec(i) = MatLib(mat_id).Sigma;
        end
        target_indices = find(Coeff_vec > 1e-12);
    else
        Coeff_vec(:) = 1.0;
        target_indices = 1:numElems;
    end
    
    numTarget = length(target_indices);
    if numTarget == 0
        C = sparse(numEdges, numActiveV);
        return;
    end

    % 2. 预分块
    chunkSize = 5000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    T_raw = Mesh.T; T2E_raw = Mesh.T2E; Signs_raw = double(Mesh.T2E_Sign);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        g_idxs = target_indices(s_idx:e_idx);
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.T2E = T2E_raw(:, g_idxs);
        ElementBlocks{k}.Signs = Signs_raw(:, g_idxs);
        ElementBlocks{k}.Coeff = Coeff_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    C_VMap = parallel.pool.Constant(LocalVMap);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        local_T = Block.T; local_T2E = Block.T2E; local_Signs = Block.Signs;
        local_Coeff = Block.Coeff;
        local_P_all = C_P.Value;
        local_VMap = C_VMap.Value;
        
        i_block = zeros(24 * count, 1);
        j_block = zeros(24 * count, 1);
        v_block = zeros(24 * count, 1);
        ptr = 0;
        
        for e = 1:count
            nodes = local_T(:, e);
            v_indices = local_VMap(nodes);
            if any(v_indices == 0), continue; end % 仅组装全导电单元
            
            val = local_Coeff(e);
            p_elem = local_P_all(:, nodes);
            Ce = Ce_edge_grad(p_elem, val);
            s = local_Signs(:, e);
            Ce_corr = Ce .* s;
            
            row_indices = local_T2E(:, e);
            
            R_idx = repmat(row_indices, 1, 4);
            C_idx = repmat(v_indices', 6, 1);
            
            range = ptr+1 : ptr+24;
            i_block(range) = R_idx(:);
            j_block(range) = C_idx(:);
            v_block(range) = Ce_corr(:);
            ptr = ptr + 24;
        end
        
        if ptr < 24*count
            i_block = i_block(1:ptr); j_block = j_block(1:ptr); v_block = v_block(1:ptr);
        end
        
        I_cell{k} = i_block; J_cell{k} = j_block; V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:}); J = vertcat(J_cell{:}); V = vertcat(V_cell{:});
    C = sparse(I, J, V, numEdges, numActiveV);
end