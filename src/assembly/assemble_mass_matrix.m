function M = assemble_mass_matrix(Model, region_type)
% ASSEMBLE_MASS_MATRIX 并行组装磁场质量矩阵 (Edge-Edge)
% M_ij = int( N_i * alpha * N_j )
%
% 输入: 
%   Model
%   region_type - (可选) 'All' 或 'ConductorOnly'. 默认 'All'.
%                 如果是 'ConductorOnly'，则只在导体区积分电导率。
%                 如果是 'All' (且未指定 sigma)，通常用于 Regularization，设 alpha=1。

    if nargin < 2, region_type = 'All'; end

    fprintf('正在组装质量矩阵 (%s)...\n', region_type);
    t_start = tic;

    Mesh = Model.Mesh;
    C_P = parallel.pool.Constant(Mesh.P);
    
    T = Mesh.T;
    T2E = Mesh.T2E;
    Signs = double(Mesh.T2E_Sign);
    
    numElems = size(T, 2);
    numEdges = size(Mesh.Edges, 2);
    
    % 准备系数 alpha
    % Case A: 真实的电导率矩阵 (sigma)
    % Case B: 正则化单位矩阵 (alpha = 1)
    
    Alpha_element = zeros(1, numElems);
    
    if strcmpi(region_type, 'ConductorOnly')
        % 物理模式: 仅在导电区有值
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        for i = 1:numElems
            mat_id = MatMap(i);
            Alpha_element(i) = MatLib(mat_id).Sigma;
        end
    else
        % 正则化模式: 全域为 1.0 (或者很小的值，由外部控制乘子)
        Alpha_element(:) = 1.0; 
    end
    
    % 如果全是0 (例如全是空气且计算涡流)，直接返回零矩阵
    if all(Alpha_element == 0)
        M = sparse(numEdges, numEdges);
        fprintf('质量矩阵全为零 (无导体).\n');
        return;
    end
    
    % 分块 (同上)
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    for k = 1:numChunks
        idx = (k-1)*chunkSize + 1 : min(k*chunkSize, numElems);
        ElementBlocks{k}.T = T(:, idx);
        ElementBlocks{k}.T2E = T2E(:, idx);
        ElementBlocks{k}.Signs = Signs(:, idx);
        ElementBlocks{k}.Alpha = Alpha_element(idx);
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
        local_Alpha = Block.Alpha;
        
        i_block = zeros(36 * count, 1);
        j_block = zeros(36 * count, 1);
        v_block = zeros(36 * count, 1);
        idx_ptr = 0;
        
        for e = 1:count
            val = local_Alpha(e);
            if val == 0, idx_ptr = idx_ptr + 36; continue; end % 跳过零值以保持稀疏性
            
            nodes = local_T(:, e);
            p_elem = local_P_all(:, nodes);
            
            % 调用内核 Me_edge_edge(P, sigma)
            Me = Me_edge_edge(p_elem, val);
            
            % 符号修正 (Rows and Cols)
            s = local_Signs(:, e);
            Me_corrected = Me .* (s * s');
            
            row_indices = local_T2E(:, e);
            R_idx = repmat(row_indices, 1, 6);
            C_idx = repmat(row_indices', 6, 1);
            
            range = idx_ptr+1 : idx_ptr+36;
            i_block(range) = R_idx(:);
            j_block(range) = C_idx(:);
            v_block(range) = Me_corrected(:);
            
            idx_ptr = idx_ptr + 36;
        end
        
        % 截断未使用的空间 (如果是导体模式，很多单元是0)
        % find(v_block ~= 0) 会比较慢，这里简单处理
        mask = v_block ~= 0;
        I_cell{k} = i_block(mask);
        J_cell{k} = j_block(mask);
        V_cell{k} = v_block(mask);
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    M = sparse(I, J, V, numEdges, numEdges);
    
    t_end = toc(t_start);
    fprintf('质量矩阵组装完成。耗时: %.4f 秒, 非零元: %d\n', t_end, nnz(M));
end