function F = assemble_source_kernel(PackedData, SourceDataMap, Config)
% ASSEMBLE_SOURCE_KERNEL 并行计算电流源载荷向量 (修正版)
% F_i = Integral( N_i . J ) dV
%
% 修复:
% 1. [Critical] 修正了 Piola 变换矩阵，由 iJt 改为 iJt' (J^{-1})，
%    确保行向量右乘后等效于列向量的 J^{-T} 变换。

    if nargin < 3, Config = FemConfig.Default(); end

    % --- 1. 准备常量数据 ---
    C_P      = parallel.pool.Constant(PackedData.P);
    C_T      = parallel.pool.Constant(PackedData.T);
    C_Dofs   = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs  = parallel.pool.Constant(PackedData.Signs);
    C_Tags   = parallel.pool.Constant(PackedData.RegionTags);
    C_SrcMap = parallel.pool.Constant(SourceDataMap);
    
    % --- 2. 准备积分数据 ---
    order = max(Config.DefaultQuadratureOrder, 2); 
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [val_ref_raw, ~] = nedelec_tet_p1(q_pts); 
    
    C_Qw     = parallel.pool.Constant(q_w_raw);
    C_ValRef = parallel.pool.Constant(val_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    F_idx_cell = cell(numChunks, 1);
    F_val_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T_all = C_T.Value;
        loc_Dofs_all = C_Dofs.Value;
        loc_Signs_all = C_Signs.Value;
        loc_Tags_all = C_Tags.Value;
        loc_P_all = C_P.Value;
        loc_SrcMap = C_SrcMap.Value;
        
        loc_qw = C_Qw.Value;
        loc_val_ref = C_ValRef.Value;
        
        f_idx_block = zeros(6 * num_local, 1);
        f_val_block = zeros(6 * num_local, 1);
        ptr = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            tag = loc_Tags_all(e);
            
            if tag > size(loc_SrcMap, 1) || tag <= 0, continue; end
            J_vec = loc_SrcMap(tag, :); 
            if all(J_vec == 0), continue; end
            
            node_ids = loc_T_all(:, e);
            p_elem = loc_P_all(:, node_ids);
            
            [~, detJ, iJt] = compute_jacobian_tet(p_elem);
            
            Fe = zeros(6, 1);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                N_ref = reshape(loc_val_ref(:, :, q), 6, 3);
                
                % [修正点] Piola 变换: N_phy = N_ref * J^{-1}
                % 因为 iJt = J^{-T}，所以我们需要 iJt' = J^{-1}
                % 行向量变换: Row_phy = Row_ref * J^{-1}
                % 等价于列向量: Col_phy = J^{-T} * Col_ref (正确的 H(curl) 变换)
                N_phy = N_ref * iJt';
                
                dot_val = N_phy * J_vec';
                Fe = Fe + w * dot_val;
            end
            
            s = double(loc_Signs_all(:, e));
            Fe = Fe .* s;
            
            dofs = loc_Dofs_all(:, e);
            range = ptr+1 : ptr+6;
            f_idx_block(range) = dofs;
            f_val_block(range) = Fe;
            ptr = ptr + 6;
        end
        
        F_idx_cell{k} = f_idx_block(1:ptr);
        F_val_cell{k} = f_val_block(1:ptr);
    end
    
    F_idx_all = vertcat(F_idx_cell{:});
    F_val_all = vertcat(F_val_cell{:});
    
    % 自动推断总自由度
    numTotalDofs = 0;
    if ~isempty(PackedData.CellDofs)
        numTotalDofs = max(PackedData.CellDofs(:));
    end
    
    if isempty(F_idx_all)
        F = sparse(numTotalDofs, 1);
    else
        F = sparse(F_idx_all, 1, F_val_all, numTotalDofs, 1);
    end
end