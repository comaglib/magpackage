function K = assemble_magnetic_stiffness(Model)
% ASSEMBLE_MAGNETIC_STIFFNESS 并行组装磁场刚度矩阵 (HPC优化版)
% 
% 优化点:
% 1. 使用 Cell 数组预分块，消除 T, T2E, Signs 的广播开销。
% 2. 使用 parallel.pool.Constant 封装 P，优化节点坐标的广播效率。

    fprintf('正在组装刚度矩阵 (并行模式 + HPC优化)...\n');
    t_start = tic;

    % 1. 提取基础数据
    Mesh = Model.Mesh;
    % P 将被封装，此处先不提取
    T = Mesh.T;           
    T2E = Mesh.T2E;       
    Signs = double(Mesh.T2E_Sign); 
    
    numElems = size(T, 2);
    numEdges = size(Mesh.Edges, 2);
    
    % 2. 准备材料参数 (计算每个单元的磁阻率 Nu)
    if isfield(Model.Materials, 'ActiveMap')
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        % 预计算 Nu 向量 (SoA)
        % 假设线性，非线性将在迭代循环中处理
        Nu_element = zeros(1, numElems);
        for i = 1:numElems
            mat_id = MatMap(i);
            Nu_element(i) = 1 / (MatLib(mat_id).Mu_r * 4*pi*1e-7);
        end
    else
        % 默认空气
        nu_air = 1/(4*pi*1e-7);
        Nu_element = repmat(nu_air, 1, numElems);
    end

    % 3. 数据预分块 (Pre-chunking)
    % 将所有基于单元的数据切分放入 Cell 数组
    % 这样 parfor 就会识别为 Slice 变量，而非 Broadcast
    
    chunkSize = 5000; % 每个任务包的大小
    numChunks = ceil(numElems / chunkSize);
    
    ElementBlocks = cell(numChunks, 1);
    
    fprintf('  - 数据预分块 (Chunks: %d)...\n', numChunks);
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numElems);
        indices = s_idx:e_idx;
        
        % 将切片存入结构体
        ElementBlocks{k}.T     = T(:, indices);
        ElementBlocks{k}.T2E   = T2E(:, indices);
        ElementBlocks{k}.Signs = Signs(:, indices);
        ElementBlocks{k}.Nu    = Nu_element(indices);
        ElementBlocks{k}.Count = length(indices);
    end
    
    % 4. 优化节点坐标 P 的广播
    % P 是随机访问的，无法切片。使用 Constant 显式告知 Pool 这是一个常数广播
    % 这样数据只传输一次，且在 Worker 端缓存
    C_P = parallel.pool.Constant(Mesh.P);

    % 5. 开启并行循环
    % 预分配结果容器
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);

    parfor k = 1:numChunks
        % 5.1 获取 Worker 本地数据
        % 访问 Constant 变量的值
        local_P_all = C_P.Value; 
        
        % 获取当前块的数据 (无通信开销，因为是切片输入)
        Block = ElementBlocks{k};
        
        local_T     = Block.T;
        local_T2E   = Block.T2E;
        local_Signs = Block.Signs;
        local_Nu    = Block.Nu;
        count       = Block.Count;
        
        % 5.2 预分配计算缓存
        i_block = zeros(36 * count, 1);
        j_block = zeros(36 * count, 1);
        v_block = zeros(36 * count, 1);
        
        idx_ptr = 0;
        
        % 5.3 密集计算内核
        for e = 1:count
            % 获取当前单元的节点坐标
            nodes = local_T(:, e);
            p_elem = local_P_all(:, nodes); % 随机访问本地缓存的 P
            
            % 物理参数
            nu_val = local_Nu(e);
            
            % Kernel 计算
            Ke = Ke_curl_curl(p_elem, nu_val);
            
            % 符号修正 (利用 Sign 向量外积)
            s = local_Signs(:, e); 
            Ke_corrected = Ke .* (s * s');
            
            % 填充索引
            row_indices = local_T2E(:, e);
            
            % 快速生成索引对 (比 meshgrid 更快的方法)
            % Col 优先: 11, 21, 31...
            R_idx = repmat(row_indices, 1, 6); % 6x6
            C_idx = repmat(row_indices', 6, 1); 
            
            % 存入 Buffer
            range = idx_ptr+1 : idx_ptr+36;
            i_block(range) = R_idx; % 自动按列展开
            j_block(range) = C_idx;
            v_block(range) = Ke_corrected;
            
            idx_ptr = idx_ptr + 36;
        end
        
        % 存入 Cell
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    % 6. Reduce 阶段
    fprintf('  - 正在合并并行数据...\n');
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    fprintf('  - 构建稀疏矩阵...\n');
    K = sparse(I, J, V, numEdges, numEdges);
    
    t_end = toc(t_start);
    fprintf('刚度矩阵组装完成。耗时: %.4f 秒, 非零元: %d\n', t_end, nnz(K));
end