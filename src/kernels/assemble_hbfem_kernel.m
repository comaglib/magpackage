function [I, J, V, Res] = assemble_hbfem_kernel(PackedData, SolHarmonics, AFT_Obj, MatLibData, Config, RequestJacobian)
% ASSEMBLE_HBFEM_KERNEL 并行组装 HBFEM 系统的刚度矩阵与残差 (Optimized)
%
% 功能:
%   1. [Critical] 修正了 nu_dc 的计算方式。
%      之前错误地使用 nu_freq(1) (对应 h=1) 作为 DC 分量。
%      当谐波列表不含 0 时，这会导致使用基波磁阻率(接近0)作为刚度，导致奇异。
%      现在直接使用 mean(nu_time) 计算真实的 DC 磁阻率。
%   2. 保留了 Jacobian 缩放 (AC=2.0) 和实数投影修复。
%   3. 增加 RequestJacobian 标志。若为 false，跳过刚度矩阵计算 (只算 Res)。

    if nargin < 5, Config = FemConfig.Default(); end
    if nargin < 6, RequestJacobian = true; end % 默认计算雅可比

    % --- 1. 准备常量 ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Tags  = parallel.pool.Constant(PackedData.RegionTags);
    
    C_Sol   = parallel.pool.Constant(SolHarmonics); 
    C_Dmat  = parallel.pool.Constant(AFT_Obj.D_matrix); 
    C_Pmat  = parallel.pool.Constant(AFT_Obj.P_matrix); 
    C_MatLib = parallel.pool.Constant(MatLibData);
    
    % Jacobian Scaling (AC=2.0, DC=1.0)
    H_vec = AFT_Obj.Harmonics;
    Scalings = ones(length(H_vec), 1);
    Scalings(H_vec > 0) = 2.0; 
    C_Scalings = parallel.pool.Constant(Scalings);
    
    % --- 2. 积分设置 ---
    order = max(Config.DefaultQuadratureOrder, 2); 
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [~, curl_ref_raw] = nedelec_tet_p1(q_pts); 
    
    C_Qw = parallel.pool.Constant(q_w_raw);
    C_CurlRef = parallel.pool.Constant(curl_ref_raw);
    
    numElems = size(PackedData.T, 2);
    numHarm = AFT_Obj.NumHarmonics;
    n_q = length(q_w_raw);
    
    % --- 3. 索引模板 ---
    base_idx = (1:6)';
    tpl_I_sub = repmat(base_idx, 6, 1);
    tpl_J_sub = repelem(base_idx, 6, 1);
    C_Tpl_I = parallel.pool.Constant(tpl_I_sub);
    C_Tpl_J = parallel.pool.Constant(tpl_J_sub);
    
    % --- 4. 并行循环 ---
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    Res_idx_cell = cell(numChunks, 1);
    Res_val_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T = C_T.Value; loc_Dofs = C_Dofs.Value; loc_Signs = C_Signs.Value;
        loc_P = C_P.Value; loc_Tags = C_Tags.Value; loc_Sol = C_Sol.Value;
        loc_Dmat = C_Dmat.Value; loc_Pmat = C_Pmat.Value;
        loc_MatLib = C_MatLib.Value;
        
        loc_qw = C_Qw.Value; loc_curl_ref = C_CurlRef.Value;
        tpl_I = C_Tpl_I.Value; tpl_J = C_Tpl_J.Value;
        loc_Scalings = C_Scalings.Value;
        
        % 预分配 (根据 RequestJacobian 决定是否分配 K 相关内存)
        if RequestJacobian
            alloc_size = 36 * num_local;
            i_block = zeros(alloc_size, 1);
            j_block = zeros(alloc_size, 1);
            v_block = zeros(alloc_size, numHarm); 
        else
            i_block = []; j_block = []; v_block = [];
        end
        
        alloc_R = 6 * num_local;
        r_idx_block = zeros(alloc_R, 1); 
        r_val_block = zeros(alloc_R, numHarm);
        
        ptr_K = 0;
        ptr_R = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            
            node_ids = loc_T(:, e);
            p_elem = loc_P(:, node_ids);
            dofs = loc_Dofs(:, e);
            s = double(loc_Signs(:, e));
            
            A_elem_global = loc_Sol(dofs, :); 
            A_elem_local  = s .* A_elem_global; 
            
            A_time = real(A_elem_local * loc_Dmat.'); 
            
            [J_mat, detJ] = compute_jacobian_tet(p_elem);
            invDetJ = 1.0 / detJ;
            
            tag = loc_Tags(e);
            matData = loc_MatLib(tag);
            
            if RequestJacobian
                Ke_harm_stack = zeros(6, 6, numHarm);
            else
                Ke_harm_stack = []; % 占位
            end
            Re_harm_stack = zeros(6, numHarm);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                C_ref = reshape(loc_curl_ref(:, :, q), 6, 3);
                C_phy = (C_ref * J_mat') * invDetJ; 
                
                B_time = C_phy' * A_time; 
                B_sq_time = sum(B_time.^2, 1); 
                B_sq_time = real(B_sq_time);
                
                [nu_time, ~] = MaterialLib.evaluate(B_sq_time, matData);
                
                H_time = nu_time .* B_time;
                
                % Freq Transform (Residual)
                H_freq = H_time * loc_Pmat.'; 
                Re_harm_stack = Re_harm_stack + w * (C_phy * H_freq);
                
                % Jacobian Approximation (Only if requested)
                if RequestJacobian
                    nu_dc = mean(nu_time); 
                    
                    Ke_base = w * nu_dc * (C_phy * C_phy'); 
                    
                    for h = 1:numHarm
                        factor = loc_Scalings(h); 
                        Ke_harm_stack(:, :, h) = Ke_harm_stack(:, :, h) + (Ke_base * factor);
                    end
                end
            end
            
            % --- Store Residual ---
            Re_harm_corr = s .* Re_harm_stack;
            range_vec = ptr_R+1 : ptr_R+6;
            r_idx_block(range_vec) = dofs;
            r_val_block(range_vec, :) = Re_harm_corr;
            ptr_R = ptr_R + 6;
            
            % --- Store Jacobian ---
            if RequestJacobian
                range_K = ptr_K+1 : ptr_K+36;
                i_block(range_K) = dofs(tpl_I);
                j_block(range_K) = dofs(tpl_J);
                
                S_mat = s * s';
                vals_corrected = zeros(36, numHarm);
                for h=1:numHarm
                    tmp = Ke_harm_stack(:,:,h) .* S_mat;
                    vals_corrected(:,h) = tmp(:);
                end
                
                v_block(range_K, :) = vals_corrected;
                ptr_K = ptr_K + 36;
            end
        end
        
        if RequestJacobian
            I_cell{k} = i_block(1:ptr_K);
            J_cell{k} = j_block(1:ptr_K);
            V_cell{k} = v_block(1:ptr_K, :);
        end
        
        Res_idx_cell{k} = r_idx_block(1:ptr_R);
        Res_val_cell{k} = r_val_block(1:ptr_R, :);
    end
    
    % --- 5. 聚合 ---
    if RequestJacobian
        I_base = vertcat(I_cell{:});
        J_base = vertcat(J_cell{:});
        V_all  = vertcat(V_cell{:}); 
        
        numDofs = max(PackedData.CellDofs(:));
        totalNNZ = length(I_base) * numHarm;
        I_final = zeros(totalNNZ, 1);
        J_final = zeros(totalNNZ, 1);
        V_final = zeros(totalNNZ, 1);
        
        ptr = 0;
        len = length(I_base);
        for h = 1:numHarm
            offset = (h-1) * numDofs;
            range = ptr+1 : ptr+len;
            I_final(range) = I_base + offset;
            J_final(range) = J_base + offset;
            V_final(range) = V_all(:, h);
            ptr = ptr + len;
        end
        I = I_final; J = J_final; V = V_final;
    else
        I = []; J = []; V = [];
    end
    
    R_idx = vertcat(Res_idx_cell{:}); 
    R_val = vertcat(Res_val_cell{:}); 
    
    numDofs = max(PackedData.CellDofs(:));
    n_entries = length(R_idx);
    Rows = repmat(R_idx, numHarm, 1);
    Cols = repelem(1:numHarm, n_entries, 1); 
    Cols = Cols(:); 
    Vals = R_val(:);
    Res = sparse(Rows, Cols, Vals, numDofs, numHarm);
end