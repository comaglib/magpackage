function [TotalLoss, LossDensity] = compute_ohmic_loss(Model, A_curr, V_curr, Param3, Param4)
% COMPUTE_OHMIC_LOSS 计算欧姆损耗 (支持 瞬态/时谐, 修正版)
%
% 用法 1 (瞬态): 
%   [P, dens] = compute_ohmic_loss(Model, A_curr, A_prev, V_curr, dt)
%   计算瞬时功率 P(t)
%
% 用法 2 (时谐):
%   [P, dens] = compute_ohmic_loss(Model, A_phasor, [], V_phasor, freq)
%   计算平均功率 P_avg (假设 A, V 为峰值相量)
%
% 输入:
%   Model   - 模型结构体
%   A_curr  - 当前 A (N_edges x 1)
%   V_curr  - 当前 V (N_active_nodes x 1, 压缩格式)
%   Param3  - (瞬态) A_prev  | (时谐) [] 空占位
%   Param4  - (瞬态) dt      | (时谐) Frequency (Hz)
% COMPUTE_OHMIC_LOSS 计算欧姆损耗 ()
%
% 修复: 改进模式识别逻辑，防止零解(实数)时误判为瞬态模式
% 1. 确保 dt 和 omega 在所有分支路径下都被定义，防止 parfor 广播错误。
% 2. 优化了模式识别逻辑。
    
    
    % 初始化默认值以满足 parfor 静态分析要求
    dt = 1.0;     % 默认占位
    omega = 0.0;  % 默认占位
    
    % 1. 解析输入模式
    if isempty(Param3)
        % --- AC Mode ---
        mode = 'AC';
        freq = Param4;
        omega = 2 * pi * freq;
        TimeFactor = 0.5; % Peak -> RMS Power
        
        % 即使在 AC 模式，保持 dt 有值，防止 parfor 报错
        A_prev = []; 
    else
        % --- Transient Mode ---
        mode = 'Transient';
        A_prev = Param3;
        dt = Param4;
        TimeFactor = 1.0;
    end

    Mesh = Model.Mesh;
    P = Mesh.P; T = Mesh.T; T2E = Mesh.T2E; Signs = double(Model.Mesh.T2E_Sign);
    numElems = size(T, 2);
    
    % 2. 恢复 V 的自由度映射
    DofData = create_dof_map(Model);
    GlobalMap = DofData.Node2DoF;
    V_offset = DofData.NumEdges;
    LocalVMap = zeros(size(GlobalMap));
    mask = GlobalMap > 0;
    LocalVMap(mask) = GlobalMap(mask) - V_offset;
    
    % 3. 识别导电单元
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
    
    % 4. 预分块
    chunkSize = 5000;
    numChunks = ceil(numTarget / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numTarget);
        g_idxs = target_elems(s_idx:e_idx);
        
        ElementBlocks{k}.GlobalIdx = g_idxs;
        ElementBlocks{k}.T = T(:, g_idxs);
        ElementBlocks{k}.T2E = T2E(:, g_idxs);
        ElementBlocks{k}.Signs = Signs(:, g_idxs);
        ElementBlocks{k}.Sigma = Sigma_vec(g_idxs);
        ElementBlocks{k}.Count = length(g_idxs);
    end
    
    % 广播常量
    C_P = parallel.pool.Constant(P);
    C_A_curr = parallel.pool.Constant(A_curr);
    C_V_curr = parallel.pool.Constant(V_curr);
    C_VMap = parallel.pool.Constant(LocalVMap);
    
    % 瞬态历史 A (AC模式为空)
    C_A_prev = parallel.pool.Constant(A_prev);
    
    Loss_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T; loc_T2E = Block.T2E; loc_Signs = Block.Signs;
        loc_Sigma = Block.Sigma;
        
        loc_P_all = C_P.Value;
        loc_A_now = C_A_curr.Value;
        loc_V_now = C_V_curr.Value;
        loc_VMap  = C_VMap.Value;
        loc_A_old = C_A_prev.Value;
        
        res_loss = zeros(count, 1);
        
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            
            % 几何导数
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            J_mat = [v21, v31, v41];
            detJ = det(J_mat);
            invJ_T = inv(J_mat)';
            Vol = abs(detJ) / 6.0;
            
            G_ref = [-1 1 0 0; -1 0 1 0; -1 0 0 1];
            G_phy = invJ_T * G_ref; 
            
            % 1. Grad V
            v_indices = loc_VMap(nodes);
            if any(v_indices == 0)
                V_vals = zeros(4, 1); 
            else
                V_vals = loc_V_now(v_indices);
            end
            GradV = G_phy * V_vals;
            
            % 2. dA/dt 或 jw*A
            edge_idx = loc_T2E(:, e);
            signs = loc_Signs(:, e);
            
            a_edges = loc_A_now(edge_idx) .* signs;
            
            dA_dt_vec = zeros(3, 1);
            dA_edges_dt = zeros(6, 1);
            
            if strcmp(mode, 'AC')
                % AC: j * omega * A
                dA_edges_dt = 1i * omega * a_edges;
            else
                % Transient: (A - A_old) / dt
                % 这里 dt 已经在外部定义，parfor 可以安全访问
                if ~isempty(loc_A_old)
                    a_edges_old = loc_A_old(edge_idx) .* signs;
                    dA_edges_dt = (a_edges - a_edges_old) / dt;
                end
            end
            
            % 插值到中心
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            for i = 1:6
                na = edge_pairs(i, 1); nb = edge_pairs(i, 2);
                N_val = 0.25 * (G_phy(:, nb) - G_phy(:, na));
                dA_dt_vec = dA_dt_vec + dA_edges_dt(i) * N_val;
            end
            
            % E = -dA/dt - grad V
            E_vec = -dA_dt_vec - GradV;
            
            % P = sigma * |E|^2 * Vol
            E_sq = sum(E_vec .* conj(E_vec));
            
            res_loss(e) = TimeFactor * loc_Sigma(e) * real(E_sq) * Vol;
        end
        
        Loss_cell{k} = res_loss;
    end
    
    LossDensity = zeros(numElems, 1);
    for k = 1:numChunks
        LossDensity(ElementBlocks{k}.GlobalIdx) = Loss_cell{k};
    end
    TotalLoss = sum(LossDensity);
end