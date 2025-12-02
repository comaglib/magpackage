function [I, J, V] = assemble_mass_kernel(PackedData, Config)
% ASSEMBLE_MASS_KERNEL 并行计算质量矩阵 M (v3.1 - Weighted Support)
% M_ij = Integral ( Coeff * Ni . Nj ) dV
%
% 输入: PackedData
%   .P, .T, .CellDofs, .Signs
%   .Coeff - (可选) [Ne x 1] 单元系数向量 (例如 Sigma)。默认为 1.0。

    if nargin < 2, Config = FemConfig.Default(); end

    % --- 1. 准备常量 ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    
    % 处理系数 (Coeff)
    if isfield(PackedData, 'Coeff') && ~isempty(PackedData.Coeff)
        HasCoeff = true;
        C_Coeff = parallel.pool.Constant(PackedData.Coeff);
    else
        HasCoeff = false;
        C_Coeff = [];
    end
    
    % --- 2. 积分设置 ---
    % N.N 是二次的，使用 Order 2 (4点积分)
    order = max(Config.DefaultQuadratureOrder, 2);
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [val_ref_raw, ~] = nedelec_tet_p1(q_pts); % [6 x 3 x Nq]
    
    C_Qw     = parallel.pool.Constant(q_w_raw);
    C_ValRef = parallel.pool.Constant(val_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    % --- 3. 索引模板 ---
    base_idx = (1:6)';
    template_I = repmat(base_idx, 6, 1);
    template_J = repelem(base_idx, 6, 1);
    C_Tpl_I = parallel.pool.Constant(template_I);
    C_Tpl_J = parallel.pool.Constant(template_J);
    
    % --- 4. 并行循环 ---
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T_all = C_T.Value;
        loc_Dofs_all = C_Dofs.Value;
        loc_Signs_all = C_Signs.Value;
        loc_P_all = C_P.Value;
        
        loc_qw = C_Qw.Value;
        loc_val_ref = C_ValRef.Value;
        tpl_I = C_Tpl_I.Value;
        tpl_J = C_Tpl_J.Value;
        
        % 获取本地系数
        if HasCoeff
            loc_Coeff_all = C_Coeff.Value;
            loc_c = loc_Coeff_all(range_indices);
        else
            loc_c = ones(num_local, 1); % 默认为 1
        end
        
        i_block = zeros(36 * num_local, 1);
        j_block = zeros(36 * num_local, 1);
        v_block = zeros(36 * num_local, 1);
        ptr = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            val_c = loc_c(ii);
            
            % 如果系数为 0 (例如绝缘体)，则跳过计算，节省时间
            if val_c == 0
                ptr = ptr + 36; % 仍然占位 (值为0)，或者可以选择完全跳过不存
                % 这里选择存0，保持结构一致性，sparse会自动剔除
                continue; 
            end
            
            node_ids = loc_T_all(:, e);
            p_elem = loc_P_all(:, node_ids);
            
            [~, detJ, iJt] = compute_jacobian_tet(p_elem);
            
            Me = zeros(6, 6);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                N_ref = reshape(loc_val_ref(:, :, q), 6, 3);
                
                % Piola 变换: N_phy = N_ref * J^{-1} = N_ref * iJt'
                N_phy = N_ref * iJt'; 
                
                % 累加: w * c * (N_phy * N_phy')
                Me = Me + (w * val_c) * (N_phy * N_phy');
            end
            
            s = double(loc_Signs_all(:, e));
            Me = Me .* (s * s');
            
            dofs = loc_Dofs_all(:, e);
            I_local = dofs(tpl_I);
            J_local = dofs(tpl_J);
            
            range = ptr+1 : ptr+36;
            i_block(range) = I_local;
            j_block(range) = J_local;
            v_block(range) = Me(:);
            ptr = ptr + 36;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
end