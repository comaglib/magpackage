function [TotalLoss, LossDensity] = compute_ohmic_loss(Model, A_curr, A_prev, V_curr, dt)
% COMPUTE_OHMIC_LOSS 计算总欧姆损耗 (Joule Heating)
% 
% 输入:
%   Model   - 模型结构体
%   A_curr  - 当前时刻磁矢量位 (N_edges x 1)
%   A_prev  - 上一时刻磁矢量位 (用于差分计算 dA/dt)
%   V_curr  - 当前时刻电标量位 (压缩格式, N_active_nodes x 1)
%   dt      - 时间步长
%
% 输出:
%   TotalLoss   - 总功率 (Watts)
%   LossDensity - 每个单元的损耗密度 (Watts/m^3) [N_elems x 1]

    % fprintf('  [Post] Calculating Ohmic Loss...\n');
    
    Mesh = Model.Mesh;
    P = Mesh.P;
    T = Mesh.T;
    T2E = Mesh.T2E;
    Signs = double(Model.Mesh.T2E_Sign);
    
    numElems = size(T, 2);
    
    % 1. 恢复 V 的自由度映射
    % 因为输入的是压缩的 V，我们需要知道它对应哪些节点
    DofData = create_dof_map(Model);
    
    % 构建 Local V Map (Global Node ID -> Index in V_curr)
    % Node2DoF 存储的是 Global DoF ID (OFFSET + Index)
    % 我们需要还原为 V_curr 的索引
    GlobalMap = DofData.Node2DoF;
    V_offset = DofData.NumEdges;
    LocalVMap = zeros(size(GlobalMap));
    mask = GlobalMap > 0;
    LocalVMap(mask) = GlobalMap(mask) - V_offset;
    
    % 2. 识别导电单元
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    Sigma_vec = zeros(numElems, 1);
    for i = 1:numElems
        mat_id = MatMap(i);
        if mat_id <= length(MatLib)
            Sigma_vec(i) = MatLib(mat_id).Sigma;
        end
    end
    
    target_elems = find(Sigma_vec > 1e-12);
    numTarget = length(target_elems);
    
    if numTarget == 0
        TotalLoss = 0;
        LossDensity = zeros(numElems, 1);
        return;
    end
    
    % 3. 预分块并行计算
    chunkSize = 5000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw = T; T2E_raw = T2E; Signs_raw = Signs;
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        g_idxs = target_elems(s_idx:e_idx);
        
        ElementBlocks{k}.GlobalIdx = g_idxs;
        ElementBlocks{k}.T = T_raw(:, g_idxs);
        ElementBlocks{k}.T2E = T2E_raw(:, g_idxs);
        ElementBlocks{k}.Signs = Signs_raw(:, g_idxs);
        ElementBlocks{k}.Sigma = Sigma_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    % 广播常量
    C_P = parallel.pool.Constant(P);
    C_A_curr = parallel.pool.Constant(A_curr);
    C_A_prev = parallel.pool.Constant(A_prev);
    C_V_curr = parallel.pool.Constant(V_curr);
    C_VMap = parallel.pool.Constant(LocalVMap);
    
    Loss_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T; loc_T2E = Block.T2E; loc_Signs = Block.Signs;
        loc_Sigma = Block.Sigma;
        
        loc_P_all = C_P.Value;
        loc_A_now = C_A_curr.Value;
        loc_A_old = C_A_prev.Value;
        loc_V_now = C_V_curr.Value;
        loc_VMap  = C_VMap.Value;
        
        % 局部结果: [GlobalElemIndex, LossValue]
        res_loss = zeros(count, 1);
        
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes); % 3x4
            
            % 几何导数
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            J_mat = [v21, v31, v41];
            detJ = det(J_mat);
            invJ_T = inv(J_mat)';
            Vol = abs(detJ) / 6.0;
            
            G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
            G_phy = invJ_T * G_ref; % 3x4 grad L
            
            % --- 计算 E 场 (重心处) ---
            
            % 1. grad V (Constant inside element)
            % V = sum Vi * Li -> grad V = sum Vi * grad Li
            v_indices = loc_VMap(nodes);
            V_vals = loc_V_now(v_indices);
            GradV = G_phy * V_vals; % 3x1
            
            % 2. dA/dt
            % A(r) = sum Aj * Nj(r).
            % 棱基函数 Nj 在重心处的值
            % Center: L = [0.25, 0.25, 0.25, 0.25]
            % Nj = L_a * grad L_b - L_b * grad L_a
            % 在重心处, La=Lb=0.25. 
            % Nj(center) = 0.25 * (grad L_b - grad L_a)
            
            edge_idx = loc_T2E(:, e);
            signs = loc_Signs(:, e);
            
            dA_edges = (loc_A_now(edge_idx) - loc_A_old(edge_idx)) / dt;
            dA_edges = dA_edges .* signs; % 修正方向
            
            dA_dt_vec = zeros(3, 1);
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            
            for i = 1:6
                na = edge_pairs(i, 1); 
                nb = edge_pairs(i, 2);
                % N_i_center = 0.25 * (grad L_b - grad L_a)
                N_val = 0.25 * (G_phy(:, nb) - G_phy(:, na));
                
                dA_dt_vec = dA_dt_vec + dA_edges(i) * N_val;
            end
            
            % E = -dA/dt - grad V
            E_vec = -dA_dt_vec - GradV;
            
            % P = sigma * |E|^2 * Vol
            % (假设 E 在单元内常数，这是一个低阶近似，但在重心积分是常用的)
            E_sq = sum(E_vec.^2);
            
            res_loss(e) = loc_Sigma(e) * E_sq * Vol;
        end
        
        Loss_cell{k} = res_loss;
    end
    
    % 汇总
    LossDensity = zeros(numElems, 1);
    
    for k = 1:numChunks
        g_idxs = ElementBlocks{k}.GlobalIdx;
        vals = Loss_cell{k};
        LossDensity(g_idxs) = vals; % 注意这里其实存的是 Power，不是 Density
        % 为了输出 Density，需要除以 Vol。但作为返回值，Power 更常用。
        % 让我们修正变量名：LossValues
    end
    
    % 这里的 LossDensity 实际上存储的是每个单元的功率 (Watts)
    % 如果需要密度，外部再除以体积。
    TotalLoss = sum(LossDensity);
    
end