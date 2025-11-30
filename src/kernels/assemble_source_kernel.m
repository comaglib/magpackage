function F = assemble_source_kernel(PackedData, Config)
% ASSEMBLE_SOURCE_KERNEL 组装电流源向量 (MATLAB)
% 
% 输入:
%   PackedData.Js: [3 x NumElems] 的电流密度矩阵
%
% 更新:
%   1. [Fix] 增加 I > 0 过滤，防止 sparse 报错 (当存在非活跃单元时 DoF=0)。

    if nargin < 2, Config = FemConfig.Default(); end

    % --- Constants ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Js    = parallel.pool.Constant(PackedData.Js);
    
    % Quadrature
    order = Config.DefaultQuadratureOrder;
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [val_ref_raw, ~] = nedelec_tet_p1(q_pts); 
    
    C_Qw = parallel.pool.Constant(q_w_raw);
    C_ValRef = parallel.pool.Constant(val_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range = idx_start:idx_end;
        num_local = length(range);
        
        loc_P = C_P.Value; loc_T = C_T.Value; 
        loc_Dofs = C_Dofs.Value; loc_Signs = C_Signs.Value;
        loc_Js = C_Js.Value;
        loc_qw = C_Qw.Value; loc_val_ref = C_ValRef.Value;
        
        i_block = zeros(6 * num_local, 1);
        v_block = zeros(6 * num_local, 1);
        ptr = 0;
        
        for ii = 1:num_local
            e = range(ii);
            
            % 获取该单元的 J [3x1]
            J_vec = loc_Js(:, e);
            
            % 快速跳过无源或无效单元
            if norm(J_vec) < 1e-12
                % 即使跳过，我们也需要处理索引。
                % 简单策略：如果不计算，就填 0 值，但在 sparse 前我们会过滤掉索引为0的项。
                % 如果 dofs 本身是 0 (非活跃单元)，这里填入 i_block 也是 0。
                % 如果 dofs 有效但 J=0，填入 i_block 是有效值但 v_block 是 0。
                % 为了效率，直接跳过计算，但为了对齐，我们可以在最后统一过滤。
                
                dofs = loc_Dofs(:, e);
                idxs = ptr+1 : ptr+6;
                i_block(idxs) = dofs; 
                v_block(idxs) = 0;
                ptr = ptr + 6;
                continue;
            end
            
            nodes = loc_T(:, e);
            p_elem = loc_P(:, nodes);
            [J_mat, detJ] = compute_jacobian_tet(p_elem);
            invJt = inv(J_mat)';
            
            Fe = zeros(6, 1);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                N_ref = reshape(loc_val_ref(:, :, q), 6, 3);
                N_phy = N_ref * invJt'; % [6x3]
                
                dot_val = N_phy * J_vec; % [6x1]
                
                Fe = Fe + w * dot_val;
            end
            
            s = double(loc_Signs(:, e));
            Fe = Fe .* s;
            
            dofs = loc_Dofs(:, e);
            idxs = ptr+1 : ptr+6;
            i_block(idxs) = dofs; % 这里可能含有 0 (非活跃单元)
            v_block(idxs) = Fe;
            ptr = ptr + 6;
        end
        
        I_cell{k} = i_block(1:ptr);
        V_cell{k} = v_block(1:ptr);
    end
    
    I = vertcat(I_cell{:});
    V = vertcat(V_cell{:});
    
    % [CRITICAL FIX] 过滤掉索引为 0 的项
    % 1. 过滤非活跃单元 (Dof=0)
    % 2. 过滤 0 值贡献 (优化稀疏度)
    valid = (I > 0) & (abs(V) > 1e-15);
    
    I = I(valid);
    V = V(valid);
    
    if isempty(I)
        % 返回空矩阵，Assembler 会调整大小
        F = sparse([], [], [], 0, 1);
    else
        % sparse 会自动累加重复索引的值
        F = sparse(I, 1, V); 
    end
end