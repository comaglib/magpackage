function C = assemble_winding_kernel(PackedData, WindingObj, Config)
% ASSEMBLE_WINDING_KERNEL 并行计算绕组耦合向量 C
% C_i = Integral( N_i . (N/S * dir) ) dV
%
% 输入:
%   PackedData - 含 P, T, CellDofs, Signs, RegionTags
%   WindingObj - Winding 对象
% 输出:
%   C - [NumGlobalDofs x 1] 稀疏向量

    if nargin < 3, Config = FemConfig.Default(); end

    % --- 1. 准备常量 ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Tags  = parallel.pool.Constant(PackedData.RegionTags);
    
    % 绕组参数
    TargetRegion = WindingObj.RegionID;
    J_eff_vec = WindingObj.getEffectiveCurrentDensity(); % [1x3]
    
    % --- 2. 积分设置 ---
    % N 是线性的，J_eff 是常数 -> 积分是被积函数是线性的 -> Order 2 (4点) 保证精确
    order = max(Config.DefaultQuadratureOrder, 2);
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [val_ref_raw, ~] = nedelec_tet_p1(q_pts); 
    
    C_Qw     = parallel.pool.Constant(q_w_raw);
    C_ValRef = parallel.pool.Constant(val_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    % --- 3. 并行循环 ---
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T = C_T.Value; loc_Dofs = C_Dofs.Value; loc_Signs = C_Signs.Value;
        loc_Tags = C_Tags.Value; loc_P = C_P.Value;
        
        loc_qw = C_Qw.Value; loc_val_ref = C_ValRef.Value;
        
        i_block = zeros(6 * num_local, 1);
        v_block = zeros(6 * num_local, 1);
        ptr = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            
            % 仅处理属于该绕组的单元
            if loc_Tags(e) ~= TargetRegion
                continue;
            end
            
            node_ids = loc_T(:, e);
            p_elem = loc_P(:, node_ids);
            
            [~, detJ, iJt] = compute_jacobian_tet(p_elem);
            
            Ce = zeros(6, 1);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                N_ref = reshape(loc_val_ref(:, :, q), 6, 3);
                
                % Piola 变换 (Value): N = J^{-T} * N_ref
                % iJt = J^{-T} => iJt' = J^{-1}
                % 公式: vec = J^{-T} * vec_hat
                % 代码行向量: vec_row = vec_hat_row * J^{-1} = vec_hat_row * iJt'
                N_phy = N_ref * iJt';
                
                % Dot product with J_eff
                dot_val = N_phy * J_eff_vec';
                Ce = Ce + w * dot_val;
            end
            
            s = double(loc_Signs(:, e));
            Ce = Ce .* s;
            
            dofs = loc_Dofs(:, e);
            
            range = ptr+1 : ptr+6;
            i_block(range) = dofs;
            v_block(range) = Ce;
            ptr = ptr + 6;
        end
        
        I_cell{k} = i_block(1:ptr);
        V_cell{k} = v_block(1:ptr);
    end
    
    % --- 4. 聚合 ---
    I_all = vertcat(I_cell{:});
    V_all = vertcat(V_cell{:});
    
    numDofs = max(PackedData.CellDofs(:));
    
    if isempty(I_all)
        C = sparse(numDofs, 1);
    else
        C = sparse(I_all, 1, V_all, numDofs, 1);
    end
end