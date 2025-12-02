function [I, J, V] = assemble_curl_curl_kernel(PackedData, Config)
% ASSEMBLE_CURL_CURL_KERNEL 并行计算 Curl-Curl 刚度矩阵 (高性能版)
% 
% 优化详情:
% 1. [Memory] 所有大数组使用 parallel.pool.Constant 封装，杜绝隐式广播。
% 2. [Speed] 移除 meshgrid/repmat，改用预计算索引模板 (Index Templates)。
% 3. [Speed] 积分常数 (q_w, curl_ref) 放入 Constant，避免重复传输。

    if nargin < 2, Config = FemConfig.Default(); end

    % --- 1. 准备常量数据 (Constant Preparation) ---
    % 将所有只读大数组送入 Constant Pool
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Nu    = parallel.pool.Constant(PackedData.Nu);
    
    % --- 2. 准备积分与基函数数据 ---
    order = Config.DefaultQuadratureOrder;
    [q_pts, q_w_raw] = get_quadrature_data('tet', order); 
    [~, curl_ref_raw] = nedelec_tet_p1(q_pts); % [6 x 3 x Nq]
    
    % 将小数组也放入 Constant，防止每次 Task 都传输
    C_Qw      = parallel.pool.Constant(q_w_raw);
    C_CurlRef = parallel.pool.Constant(curl_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    % --- 3. 预计算索引模板 (Index Templates) ---
    % 目的：避免在循环中调用 meshgrid 或 repmat
    % Ke 是 6x6 矩阵，展平后有 36 个元素 (Column-Major)
    % I (Row Idx): [1,2,3,4,5,6, 1,2,3...]' -> 重复的列
    % J (Col Idx): [1,1,1,1,1,1, 2,2,2...]' -> 重复的元素
    
    base_idx = (1:6)';
    % template_I: [1;2;3;4;5;6; 1;2;3;4;5;6; ...]
    template_I = repmat(base_idx, 6, 1); 
    % template_J: [1;1;1;1;1;1; 2;2;2;2;2;2; ...]
    template_J = repelem(base_idx, 6, 1);
    
    % 放入 Constant 以便 Worker 使用
    C_Tpl_I = parallel.pool.Constant(template_I);
    C_Tpl_J = parallel.pool.Constant(template_J);
    
    % --- 4. 并行循环 (Parfor) ---
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        
        % 获取当前 Chunk 的索引范围
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        % --- 获取 Constant 引用 (零拷贝) ---
        % 注意：这里我们在 Worker 内部进行切片 (Slice inside Worker)
        % 虽然此时所有 Worker 都能访问整个 Constant，但只有需要的部分会被读入 Cache
        % 对于极大规模数据，这比让主线程切片发送更稳健
        
        loc_T_all = C_T.Value;
        loc_Dofs_all = C_Dofs.Value;
        loc_Signs_all = C_Signs.Value;
        loc_Nu_all = C_Nu.Value;
        loc_P_all = C_P.Value;
        
        % 本地化的小常量
        loc_qw = C_Qw.Value;
        loc_curl_ref = C_CurlRef.Value;
        tpl_I = C_Tpl_I.Value;
        tpl_J = C_Tpl_J.Value;
        
        % 预分配输出块
        i_block = zeros(36 * num_local, 1);
        j_block = zeros(36 * num_local, 1);
        v_block = zeros(36 * num_local, 1);
        
        ptr = 0;
        
        % 内部循环
        for ii = 1:num_local
            e = range_indices(ii); % 全局单元索引
            
            % 1. 获取几何
            node_ids = loc_T_all(:, e);
            p_elem = loc_P_all(:, node_ids);
            
            % 手动内联 compute_jacobian_tet 的核心部分以减少函数调用开销 (可选，视 JIT 而定)
            % 这里保持调用函数，因为 compute_jacobian_tet 已经很简单
            [J_mat, detJ] = compute_jacobian_tet(p_elem);
            
            invDetJ = 1.0 / detJ;
            nu_e = loc_Nu_all(e);
            
            % 2. 积分计算 Ke
            Ke = zeros(6, 6);
            
            for q = 1:n_q
                % 权重包含 detJ (积分变换)
                w_val = loc_qw(q) * abs(detJ);
                factor = w_val * nu_e;
                
                % 获取参考旋度 (Nq x 3) -> (6 x 3)
                % loc_curl_ref 是 [6 x 3 x Nq]
                C_ref = loc_curl_ref(:, :, q); % Matlab 会自动降维成 6x3
                
                % Piola 变换: curl_phy = (1/detJ) * J * curl_ref
                % 注意公式: u = J^{-T} u_ref => curl(u) = (1/detJ) * J * curl(u_ref)
                % 代码中: C_phy = (C_ref * J_mat') * invDetJ
                % C_ref: [6x3], J_mat': [3x3] -> [6x3]
                C_phy = (C_ref * J_mat') * invDetJ;
                
                % 累加刚度
                Ke = Ke + factor * (C_phy * C_phy');
            end
            
            % 3. 符号修正 (Int8 -> Double 转换)
            s_vec = double(loc_Signs_all(:, e)); 
            
            % 利用广播进行符号修正 (s * s') 是 6x6 矩阵
            % 比如 s=[1; -1], s*s'=[1 -1; -1 1]
            % 这一步比显式循环快得多
            Ke = Ke .* (s_vec * s_vec');
            
            % 4. 填充索引 (使用模板)
            dofs = loc_Dofs_all(:, e); % [6x1]
            
            % 直接通过模板索引生成 36x1 的 I 和 J，避免 meshgrid
            I_local = dofs(tpl_I);
            J_local = dofs(tpl_J);
            
            % 存入 Block
            range = ptr+1 : ptr+36;
            i_block(range) = I_local;
            j_block(range) = J_local;
            v_block(range) = Ke(:); % 列优先展平
            
            ptr = ptr + 36;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    % --- 5. 聚合 (Reduce) ---
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
end