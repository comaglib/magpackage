function As_vec = project_source_A_on_edges(Model, Coil)
% PROJECT_SOURCE_A_ON_EDGES 将源磁矢量位投影到网格棱边上
% 
% 输入: Model, Coil (单位电流线圈)
% 输出: As_vec (N_edges x 1) 投影值

    fprintf('正在投影源磁矢量位 (Source A Projection)...\n');
    t_start = tic;
    
    Mesh = Model.Mesh;
    P = Mesh.P;
    Edges = Mesh.Edges;
    numEdges = size(Edges, 2);
    
    % 1. 预分块处理 (Chunking)
    chunkSize = 10000;
    numChunks = ceil(numEdges / chunkSize);
    
    As_cell = cell(numChunks, 1);
    
    % 广播常量
    C_P = parallel.pool.Constant(P);
    C_Edges = parallel.pool.Constant(Edges);
    C_Coil = parallel.pool.Constant(Coil);
    
    parfor k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, numEdges);
        indices = idx_s:idx_e;
        count = length(indices);
        
        local_P = C_P.Value;
        local_Edges = C_Edges.Value(:, indices);
        local_Coil = C_Coil.Value;
        
        % 获取棱的端点
        p1 = local_P(:, local_Edges(1, :));
        p2 = local_P(:, local_Edges(2, :));
        
        % 计算棱中点 (Midpoint Rule 积分)
        mid_pts = 0.5 * (p1 + p2);
        
        % 计算棱矢量 dL
        edge_vecs = p2 - p1;
        
        % 计算中点处的 A (3 x count)
        A_vals = compute_biot_savart_A(local_Coil, mid_pts);
        
        % 投影: A dot dL
        % dot(A, B, 1) 按列求点积
        proj_vals = dot(A_vals, edge_vecs, 1)';
        
        As_cell{k} = proj_vals;
    end
    
    As_vec = vertcat(As_cell{:});
    
    t_end = toc(t_start);
    fprintf('源场投影完成。耗时: %.4f 秒\n', t_end);
end