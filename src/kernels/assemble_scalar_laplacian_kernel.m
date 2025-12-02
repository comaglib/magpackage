function [I, J, V] = assemble_scalar_laplacian_kernel(PackedData, Config)
% ASSEMBLE_SCALAR_LAPLACIAN_KERNEL 并行组装标量拉普拉斯矩阵
% K_ij = Integral ( Coeff * grad(Ni) . grad(Nj) ) dV
%
% 输入: PackedData (含 P, T, CellDofs)
%       .Coeff - [Ne x 1] 单元系数 (如电导率 sigma 或介电常数 epsilon)

    if nargin < 2, Config = FemConfig.Default(); end

    % --- 1. Constants ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    
    if isfield(PackedData, 'Coeff') && ~isempty(PackedData.Coeff)
        C_Coeff = parallel.pool.Constant(PackedData.Coeff);
        HasCoeff = true;
    else
        HasCoeff = false;
        C_Coeff = [];
    end
    
    % --- 2. Quadrature & Basis ---
    % grad(N) . grad(N) 是常数 (P1单元)，1点积分足够，但为了统一用标准规则
    order = 1; 
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    % 需要 lagrange_tet_p1 (Step 2 已创建)
    [~, grad_ref_raw] = lagrange_tet_p1(q_pts); 
    
    C_Qw = parallel.pool.Constant(q_w_raw);
    C_GradRef = parallel.pool.Constant(grad_ref_raw);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    % --- 3. Templates ---
    % Lagrange P1 有 4 个节点
    base_idx = (1:4)';
    tpl_I_sub = repmat(base_idx, 4, 1);
    tpl_J_sub = repelem(base_idx, 4, 1);
    C_Tpl_I = parallel.pool.Constant(tpl_I_sub);
    C_Tpl_J = parallel.pool.Constant(tpl_J_sub);
    
    % --- 4. Parfor ---
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
        
        loc_T = C_T.Value; 
        loc_Dofs = C_Dofs.Value;
        loc_P = C_P.Value;
        
        loc_qw = C_Qw.Value; 
        loc_grad_ref = C_GradRef.Value;
        
        if HasCoeff
            loc_Coeff_all = C_Coeff.Value;
            loc_c = loc_Coeff_all(range_indices);
        else
            loc_c = ones(num_local, 1);
        end
        
        % 预分配: 每个单元 4*4 = 16 个非零元
        alloc_size = 16 * num_local;
        i_block = zeros(alloc_size, 1);
        j_block = zeros(alloc_size, 1);
        v_block = zeros(alloc_size, 1);
        
        ptr = 0;
        tpl_I = C_Tpl_I.Value;
        tpl_J = C_Tpl_J.Value;
        
        for ii = 1:num_local
            e = range_indices(ii);
            val_c = loc_c(ii);
            
            if val_c == 0, ptr = ptr + 16; continue; end
            
            node_ids = loc_T(:, e);
            p_elem = loc_P(:, node_ids);
            
            % 计算几何 Jacobian
            % 注意: compute_jacobian_tet 返回 iJt = inv(J)'
            [~, detJ, iJt] = compute_jacobian_tet(p_elem);
            
            Ke = zeros(4, 4);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                % grad_ref: [4 x 3] (at q point)
                % loc_grad_ref: [4 x 3 x Nq]
                G_ref = reshape(loc_grad_ref(:, :, q), 4, 3);
                
                % 物理梯度: grad_phy = grad_ref * J^{-1} = grad_ref * iJt'
                G_phy = G_ref * iJt';
                
                % Integral: w * c * (grad . grad)
                % G_phy * G_phy' -> [4 x 4]
                Ke = Ke + (w * val_c) * (G_phy * G_phy');
            end
            
            dofs = loc_Dofs(:, e); % [4 x 1]
            
            range = ptr+1 : ptr+16;
            i_block(range) = dofs(tpl_I);
            j_block(range) = dofs(tpl_J);
            v_block(range) = Ke(:);
            ptr = ptr + 16;
        end
        
        I_cell{k} = i_block(1:ptr);
        J_cell{k} = j_block(1:ptr);
        V_cell{k} = v_block(1:ptr);
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
end