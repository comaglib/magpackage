function M = assemble_mass_matrix(Model, region_type)
% ASSEMBLE_MASS_MATRIX 并行组装磁场质量矩阵 (HPC加速版)
% M_ij = int( N_i * alpha * N_j )
%
% 输入: 
%   Model
%   region_type - (可选) 'All' 或 'ConductorOnly'. 默认 'All'.
%                 如果是 'ConductorOnly'，则只在导体区积分电导率。
%                 如果是 'All' (且未指定 sigma)，通常用于 Regularization，设 alpha=1。
% 
% 优化:
% 1. Early Filtering: 仅对系数非零的单元进行循环和分块，大幅提升 ConductorOnly 模式速度。
% 2. Strict Chunking: 消除所有广播变量。
    
    if nargin < 2, region_type = 'All'; end

    % fprintf('正在组装质量矩阵 (%s)...\n', region_type);
    t_start = tic;

    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    numEdges = size(Mesh.Edges, 2);
    
    % 1. 准备系数向量 (Alpha)
    Alpha_vec = zeros(1, numElems);
    
    if strcmpi(region_type, 'ConductorOnly')
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        % 预读取材料电导率
        SigmaLib = zeros(length(MatLib), 1);
        for i = 1:length(MatLib)
            SigmaLib(i) = MatLib(i).Sigma;
        end
        % 向量化映射
        Alpha_vec = SigmaLib(MatMap)'; 
    else
        % Regularization or All
        Alpha_vec(:) = 1.0; 
    end
    
    % 2. 提前筛选有效单元 (核心优化)
    % 仅处理 alpha > 0 的单元
    target_elems = find(Alpha_vec > 1e-12);
    numTarget = length(target_elems);
    
    if numTarget == 0
        M = sparse(numEdges, numEdges);
        % fprintf('  - [Info] 质量矩阵为空.\n');
        return;
    end
    
    % 3. 数据预分块
    chunkSize = 5000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    % 提取原始大数组
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        
        % 获取全局单元索引
        g_idxs = target_elems(s_idx:e_idx);
        
        % 仅切片有效数据
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.T2E = T2E_raw(:, g_idxs);
        ElementBlocks{k}.Signs = Signs_raw(:, g_idxs);
        ElementBlocks{k}.Alpha = Alpha_vec(g_idxs); % 已经是紧凑的非零值
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    % 广播 P
    C_P = parallel.pool.Constant(Mesh.P);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T;
        loc_T2E = Block.T2E;
        loc_Signs = Block.Signs;
        loc_Alpha = Block.Alpha;
        
        loc_P_all = C_P.Value;
        
        i_block = zeros(36 * count, 1);
        j_block = zeros(36 * count, 1);
        v_block = zeros(36 * count, 1);
        ptr = 0;
        
        for e = 1:count
            % 无需再判断 if val == 0，因为传入的都是非零单元
            val = loc_Alpha(e);
            nodes = loc_T(:, e);
            p_elem = loc_P_all(:, nodes);
            
            % Kernel: Me_edge_edge(P, val)
            Me = Me_edge_edge(p_elem, val);
            
            s = loc_Signs(:, e);
            Me_corr = Me .* (s * s');
            
            rows = loc_T2E(:, e);
            R_idx = repmat(rows, 1, 6);
            C_idx = repmat(rows', 6, 1);
            
            range = ptr+1 : ptr+36;
            i_block(range) = R_idx;
            j_block(range) = C_idx;
            v_block(range) = Me_corr;
            ptr = ptr + 36;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    % Reduce
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    M = sparse(I, J, V, numEdges, numEdges);
    
    t_end = toc(t_start);
    fprintf('质量矩阵组装完成 (%s)。耗时: %.4f 秒, 非零元: %d\n', ...
        region_type, t_end, nnz(M));
end