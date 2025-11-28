function [I, J, V, Residual] = assemble_jacobian_kernel(PackedData, SolutionA, MatLibData, Config, RequestJacobian)
% ASSEMBLE_JACOBIAN_KERNEL 并行组装切线刚度矩阵 J 和残差向量 R (Opt v2)
%
% 输入:
%   ...
%   RequestJacobian - (可选, 默认true) 若为 false，则跳过 I,J,V 计算，仅返回 Residual
%
% 输出:
%   I, J, V    - 雅可比矩阵三元组 (若 RequestJacobian=false，则为空)
%   Residual   - 全局残差向量

    if nargin < 4, Config = FemConfig.Default(); end
    if nargin < 5, RequestJacobian = true; end % 默认计算雅可比

    % --- 1. 准备常量 ---
    C_P     = parallel.pool.Constant(PackedData.P);
    C_T     = parallel.pool.Constant(PackedData.T);
    C_Dofs  = parallel.pool.Constant(PackedData.CellDofs);
    C_Signs = parallel.pool.Constant(PackedData.Signs);
    C_Tags  = parallel.pool.Constant(PackedData.RegionTags);
    C_SolA  = parallel.pool.Constant(SolutionA);
    C_MatLib = parallel.pool.Constant(MatLibData);
    
    order = max(Config.DefaultQuadratureOrder, 2); 
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    [~, curl_ref_raw] = nedelec_tet_p1(q_pts); 
    
    C_Qw = parallel.pool.Constant(q_w_raw);
    C_CurlRef = parallel.pool.Constant(curl_ref_raw);
    
    base_idx = (1:6)';
    template_I = repmat(base_idx, 6, 1);
    template_J = repelem(base_idx, 6, 1);
    C_Tpl_I = parallel.pool.Constant(template_I);
    C_Tpl_J = parallel.pool.Constant(template_J);
    
    numElems = size(PackedData.T, 2);
    n_q = length(q_w_raw);
    
    % --- 2. 并行循环 ---
    chunkSize = Config.AssemblyChunkSize;
    numChunks = ceil(numElems / chunkSize);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    R_idx_cell = cell(numChunks, 1);
    R_val_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        idx_start = (k-1)*chunkSize + 1;
        idx_end = min(k*chunkSize, numElems);
        range_indices = idx_start:idx_end;
        num_local = length(range_indices);
        
        loc_T_all = C_T.Value;
        loc_Dofs_all = C_Dofs.Value;
        loc_Signs_all = C_Signs.Value;
        loc_P_all = C_P.Value;
        loc_Tags_all = C_Tags.Value;
        loc_SolA = C_SolA.Value;
        loc_MatLib = C_MatLib.Value;
        
        loc_qw = C_Qw.Value;
        loc_curl_ref = C_CurlRef.Value;
        tpl_I = C_Tpl_I.Value;
        tpl_J = C_Tpl_J.Value;
        
        % 仅在需要时分配雅可比内存
        if RequestJacobian
            i_block = zeros(36 * num_local, 1);
            j_block = zeros(36 * num_local, 1);
            v_block = zeros(36 * num_local, 1);
        else
            i_block = []; j_block = []; v_block = [];
        end
        
        r_idx_block = zeros(6 * num_local, 1);
        r_val_block = zeros(6 * num_local, 1);
        
        ptr_K = 0;
        ptr_R = 0;
        Identity = eye(3);
        
        for ii = 1:num_local
            e = range_indices(ii);
            
            node_ids = loc_T_all(:, e);
            p_elem = loc_P_all(:, node_ids);
            dofs = loc_Dofs_all(:, e);
            s = double(loc_Signs_all(:, e));
            
            A_elem_dofs = loc_SolA(dofs);
            A_local = s .* A_elem_dofs; 
            
            [J_mat, detJ] = compute_jacobian_tet(p_elem);
            invDetJ = 1.0 / detJ;
            
            tag = loc_Tags_all(e);
            matData = loc_MatLib(tag);
            
            Ke_tan = zeros(6, 6); % 即使不存，计算过程可能需要
            Re = zeros(6, 1);
            
            for q = 1:n_q
                w = loc_qw(q) * abs(detJ);
                
                C_ref = reshape(loc_curl_ref(:, :, q), 6, 3);
                C_phy = (C_ref * J_mat') * invDetJ; 
                
                B_vec = C_phy' * A_local;
                B_sq = dot(B_vec, B_vec);
                
                [nu, dnu_db2] = MaterialLib.evaluate(B_sq, matData);
                
                % 计算残差 (Internal Force)
                H_vec = nu * B_vec;
                Re = Re + w * (C_phy * H_vec);
                
                % 仅当需要时计算雅可比
                if RequestJacobian
                    M_diff = nu * Identity + (2 * dnu_db2) * (B_vec * B_vec');
                    Temp = C_phy * M_diff;
                    Ke_tan = Ke_tan + w * (Temp * C_phy');
                end
            end
            
            % Residual (始终计算)
            Re = Re .* s;
            range_R = ptr_R+1 : ptr_R+6;
            r_idx_block(range_R) = dofs;
            r_val_block(range_R) = Re;
            ptr_R = ptr_R + 6;
            
            % Jacobian (可选)
            if RequestJacobian
                Ke_tan = Ke_tan .* (s * s');
                range_K = ptr_K+1 : ptr_K+36;
                i_block(range_K) = dofs(tpl_I);
                j_block(range_K) = dofs(tpl_J);
                v_block(range_K) = Ke_tan(:);
                ptr_K = ptr_K + 36;
            end
        end
        
        if RequestJacobian
            I_cell{k} = i_block(1:ptr_K);
            J_cell{k} = j_block(1:ptr_K);
            V_cell{k} = v_block(1:ptr_K);
        end
        
        R_idx_cell{k} = r_idx_block(1:ptr_R);
        R_val_cell{k} = r_val_block(1:ptr_R);
    end
    
    if RequestJacobian
        I = vertcat(I_cell{:});
        J = vertcat(J_cell{:});
        V = vertcat(V_cell{:});
    else
        I = []; J = []; V = [];
    end
    
    R_idx = vertcat(R_idx_cell{:});
    R_val = vertcat(R_val_cell{:});
    
    numDofs = length(SolutionA);
    if isempty(R_idx)
        Residual = sparse(numDofs, 1);
    else
        Residual = sparse(R_idx, 1, R_val, numDofs, 1);
    end
end