function C = assemble_coupling_matrix(Model)
% ASSEMBLE_COUPLING_MATRIX 并行组装 棱-节点 耦合矩阵
% 用于 Lagrange Multiplier 规范条件或 A-V 耦合
%
% 输入: Model
% 输出: C (Sparse Matrix, N_edges x N_nodes)
%       C_ij = int( N_edge_i * grad(Phi_node_j) )

    fprintf('正在组装耦合矩阵 (Edge-Node Coupling)...\n');
    t_start = tic;

    Mesh = Model.Mesh;
    % P 需要广播
    C_P = parallel.pool.Constant(Mesh.P);
    
    T = Mesh.T;
    T2E = Mesh.T2E;
    Signs = double(Mesh.T2E_Sign);
    
    numElems = size(T, 2);
    numEdges = size(Mesh.Edges, 2);
    numNodes = size(Mesh.P, 2);
    
    % 分块策略
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    for k = 1:numChunks
        idx = (k-1)*chunkSize + 1 : min(k*chunkSize, numElems);
        ElementBlocks{k}.T = T(:, idx);
        ElementBlocks{k}.T2E = T2E(:, idx);
        ElementBlocks{k}.Signs = Signs(:, idx);
        ElementBlocks{k}.Count = length(idx);
    end
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        local_P_all = C_P.Value;
        Block = ElementBlocks{k};
        
        count = Block.Count;
        local_T = Block.T;
        local_T2E = Block.T2E;
        local_Signs = Block.Signs;
        
        % 每个单元 6 edges * 4 nodes = 24 entries
        i_block = zeros(24 * count, 1);
        j_block = zeros(24 * count, 1);
        v_block = zeros(24 * count, 1);
        
        idx_ptr = 0;
        
        % 系数 (耦合矩阵不需要材料参数，或者是 Sigma=1 的情况)
        % Kernel Ce_edge_grad(P, sigma) 计算 int(N * sigma * grad L)
        % 我们这里设 sigma=1.0 纯几何积分
        coeff = 1.0;
        
        for e = 1:count
            nodes = local_T(:, e);
            p_elem = local_P_all(:, nodes);
            
            % 调用内核
            % Ce_local: 6x4 矩阵 (Rows: Edges, Cols: Nodes)
            Ce = Ce_edge_grad(p_elem, coeff);
            
            % 符号修正
            % 仅修正行 (Edge basis 有方向)，列 (Node basis) 无方向
            s = local_Signs(:, e); % 6x1
            % s * ones(1,4) -> 将符号作用于每一列
            Ce_corrected = Ce .* s; 
            
            % 索引
            row_indices = local_T2E(:, e); % 6x1 (Edge Global IDs)
            col_indices = nodes;           % 4x1 (Node Global IDs)
            
            % 构造网格
            % R_idx (6x4): 每一列都是 row_indices
            % C_idx (6x4): 每一行都是 col_indices'
            R_idx = repmat(row_indices, 1, 4);
            C_idx = repmat(col_indices', 6, 1);
            
            range = idx_ptr+1 : idx_ptr+24;
            i_block(range) = R_idx(:);
            j_block(range) = C_idx(:);
            v_block(range) = Ce_corrected(:);
            
            idx_ptr = idx_ptr + 24;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    C = sparse(I, J, V, numEdges, numNodes);
    
    t_end = toc(t_start);
    fprintf('耦合矩阵组装完成。耗时: %.4f 秒, 维度: %dx%d\n', ...
        t_end, size(C,1), size(C,2));
end