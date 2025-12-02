function C = assemble_winding_kernel(PackedData, WindingObj, Config)
% ASSEMBLE_WINDING_KERNEL (v2.0 - Direction Field Support)
% 
% 更新:
%   支持从 WindingObj.DirectionField 读取每个单元的独立方向向量。

    if nargin < 3, Config = FemConfig.Default(); end

    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Tags  = parallel.pool.Constant(PackedData.RegionTags);
    
    TargetRegion = WindingObj.RegionID;
    
    % 检查是否有场数据
    HasField = ~isempty(WindingObj.DirectionField);
    if HasField
        C_DirField = parallel.pool.Constant(WindingObj.DirectionField);
        ConstDir = [0,0,0]; % 占位
    else
        C_DirField = [];
        ConstDir = WindingObj.Direction;
    end
    
    % 预计算系数 N/S
    if WindingObj.CrossSectionArea <= 0, error('Invalid Area'); end
    CurrentDensityScale = WindingObj.Turns / WindingObj.CrossSectionArea;
    
    % 积分设置
    order = max(Config.DefaultQuadratureOrder, 2);
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
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T = C_T.Value; loc_Dofs = C_Dofs.Value; loc_Signs = C_Signs.Value;
        loc_Tags = C_Tags.Value; loc_P = C_P.Value;
        loc_qw = C_Qw.Value; loc_val_ref = C_ValRef.Value;
        
        if HasField
            loc_DirAll = C_DirField.Value;
            % 注意：DirectionField 是 [3 x NumElems]
            loc_Dirs = loc_DirAll(:, range_indices);
        end
        
        i_block = zeros(6 * num_local, 1);
        v_block = zeros(6 * num_local, 1);
        ptr = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            
            if loc_Tags(e) ~= TargetRegion
                continue;
            end
            
            % 获取方向向量 t
            if HasField
                t_vec = loc_Dirs(:, ii)'; % [1x3]
            else
                t_vec = ConstDir;
            end
            
            % 有效电流密度 J = (N/S) * t
            J_eff = t_vec * CurrentDensityScale;
            
            node_ids = loc_T(:, e);
            p_elem = loc_P(:, node_ids);
            [~, detJ, iJt] = compute_jacobian_tet(p_elem);
            
            Ce = zeros(6, 1);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                N_ref = reshape(loc_val_ref(:, :, q), 6, 3);
                N_phy = N_ref * iJt';
                
                dot_val = N_phy * J_eff';
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
    
    I_all = vertcat(I_cell{:});
    V_all = vertcat(V_cell{:});
    numDofs = max(PackedData.CellDofs(:));
    
    if isempty(I_all)
        C = sparse(numDofs, 1);
    else
        C = sparse(I_all, 1, V_all, numDofs, 1);
    end
end