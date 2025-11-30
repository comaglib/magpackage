function [I, J, V] = assemble_coupling_kernel(MeshP, MeshT, Dofs_Edge, Signs_Edge, Dofs_Node, Coeffs, Config)
% ASSEMBLE_COUPLING_KERNEL 计算 Nedelec-Lagrange 耦合矩阵 (高性能优化版)
%
% 优化日志:
% 1. [Speed] 移除 inv(J_mat) 和 det(J_mat) 函数调用。
%    改用显式代数公式计算 3x3 逆矩阵，消除 BLAS/LAPACK 调用开销。
%    对于 3x3 矩阵，显式计算比 inv() 快 10-20 倍。
% 2. [Memory] 使用 parallel.pool.Constant 管理大数组。
% 3. [Vectorization] 移除循环内的 cross() 函数，使用直接索引乘法。
%
% 公式: C_ij = \int sigma * (N_i . grad(phi_j)) dV

    if nargin < 7, Config = FemConfig.Default(); end

    % --- 1. 准备常量数据 (Constant Pool) ---
    C_P     = parallel.pool.Constant(MeshP);
    C_T     = parallel.pool.Constant(MeshT);
    C_DofsE = parallel.pool.Constant(Dofs_Edge);
    C_Signs = parallel.pool.Constant(Signs_Edge);
    C_DofsN = parallel.pool.Constant(Dofs_Node);
    C_Coeff = parallel.pool.Constant(Coeffs);
    
    % --- 2. 准备积分数据 ---
    order = max(2, Config.DefaultQuadratureOrder);
    [q_pts, q_w_raw] = get_quadrature_data('tet', order);
    n_q = length(q_w_raw);
    
    % 获取参考基函数值
    [ref_ned_val, ~] = nedelec_tet_p1(q_pts); % [6 x 3 x Nq]
    [~, ref_lag_grad] = lagrange_tet_p1(q_pts); % [4 x 3 x Nq]
    
    C_Qw      = parallel.pool.Constant(q_w_raw);
    C_NedVal  = parallel.pool.Constant(ref_ned_val);
    C_LagGrad = parallel.pool.Constant(ref_lag_grad);
    
    numElems = size(MeshT, 2);
    
    % --- 3. 预计算索引模板 (Index Templates) ---
    % 局部矩阵大小: [6 (Edges) x 4 (Nodes)] = 24 元素
    base_rows = (1:6)';
    base_cols = (1:4);
    template_I = repmat(base_rows, 4, 1);     % [1;2;3;4;5;6; 1...]
    template_J = repelem(base_cols', 6, 1);   % [1;1;1;1;1;1; 2...]
    
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
        
        % 提取本地数据块 (减少 Value 访问次数)
        loc_T = C_T.Value;
        loc_P = C_P.Value;
        loc_DofsE = C_DofsE.Value;
        loc_Signs = C_Signs.Value;
        loc_DofsN = C_DofsN.Value;
        loc_Coeff = C_Coeff.Value;
        
        loc_qw = C_Qw.Value;
        loc_ned = C_NedVal.Value;
        loc_lag = C_LagGrad.Value;
        
        tpl_I = C_Tpl_I.Value;
        tpl_J = C_Tpl_J.Value;
        
        i_block = zeros(24 * num_local, 1);
        j_block = zeros(24 * num_local, 1);
        v_block = zeros(24 * num_local, 1);
        
        ptr = 0;
        
        for ii = 1:num_local
            e = range_indices(ii);
            sigma = loc_Coeff(e);
            
            % 快速跳过绝缘体 (优化稀疏度)
            if abs(sigma) < 1e-14
                 % 获取 DoFs 以填 0 占位 (或完全跳过如果 Assembler 支持)
                 % 这里为了代码简单，填 0
                 dofs_E = loc_DofsE(:, e);
                 dofs_N = loc_DofsN(:, e);
                 range = ptr+1 : ptr+24;
                 i_block(range) = dofs_E(tpl_I);
                 j_block(range) = dofs_N(tpl_J);
                 v_block(range) = 0;
                 ptr = ptr + 24;
                 continue;
            end
            
            % 1. 显式获取坐标
            n_ids = loc_T(:, e);
            p1 = loc_P(:, n_ids(1));
            p2 = loc_P(:, n_ids(2));
            p3 = loc_P(:, n_ids(3));
            p4 = loc_P(:, n_ids(4));
            
            % 2. 显式计算雅可比矩阵列向量 (Edge Vectors)
            x21 = p2(1)-p1(1); y21 = p2(2)-p1(2); z21 = p2(3)-p1(3);
            x31 = p3(1)-p1(1); y31 = p3(2)-p1(2); z31 = p3(3)-p1(3);
            x41 = p4(1)-p1(1); y41 = p4(2)-p1(2); z41 = p4(3)-p1(3);
            
            % 3. 显式计算行列式 (Scalar Triple Product)
            % detJ = (v21 x v31) . v41
            % 提前计算交叉项以复用
            cr_x_1 = y21*z31 - z21*y31; % v21 x v31 (X)
            cr_y_1 = z21*x31 - x21*z31; % v21 x v31 (Y)
            cr_z_1 = x21*y31 - y21*x31; % v21 x v31 (Z)
            
            detJ = cr_x_1*x41 + cr_y_1*y41 + cr_z_1*z41;
            invDetJ = 1.0 / detJ;
            absDetJ = abs(detJ);
            
            % 4. 显式计算逆矩阵 J^{-1} (Cramer's Rule / Cross Products)
            % J = [v21, v31, v41]
            % invJ = [ (v31 x v41)^T ; (v41 x v21)^T ; (v21 x v31)^T ] / detJ
            
            % Row 1: (v31 x v41) / detJ
            rx1 = (y31*z41 - z31*y41) * invDetJ;
            ry1 = (z31*x41 - x31*z41) * invDetJ;
            rz1 = (x31*y41 - y31*x41) * invDetJ;
            
            % Row 2: (v41 x v21) / detJ
            rx2 = (y41*z21 - z41*y21) * invDetJ;
            ry2 = (z41*x21 - x41*z21) * invDetJ;
            rz2 = (x41*y21 - y41*x21) * invDetJ;
            
            % Row 3: (v21 x v31) / detJ (利用已计算的 cr_1)
            rx3 = cr_x_1 * invDetJ;
            ry3 = cr_y_1 * invDetJ;
            rz3 = cr_z_1 * invDetJ;
            
            % 构造 invJ (3x3)
            invJ = [rx1, ry1, rz1; 
                    rx2, ry2, rz2; 
                    rx3, ry3, rz3];
            
            % 5. 积分循环
            Ce = zeros(6, 4);
            
            for q = 1:n_q
                w = loc_qw(q) * absDetJ;
                factor = w * sigma;
                
                % 获取参考值 (Matlab 会自动将 3D 切片降维为 2D)
                N_ref = loc_ned(:, :, q); % [6x3]
                G_ref = loc_lag(:, :, q); % [4x3]
                
                % Piola 变换: u = u_ref * invJ (行向量规则)
                % 这里用矩阵乘法，因为 6x3 * 3x3 很快，且比显式展开更易读
                N_phy = N_ref * invJ; 
                G_phy = G_ref * invJ;
                
                % 累加: [6x4]
                % 这里是瓶颈所在，确保矩阵乘法高效
                Ce = Ce + factor * (N_phy * G_phy');
            end
            
            % 6. 符号修正 (Broadcast)
            s_vec = double(loc_Signs(:, e));
            Ce = Ce .* s_vec; 
            
            % 7. 填充
            dofs_E = loc_DofsE(:, e);
            dofs_N = loc_DofsN(:, e);
            
            range = ptr+1 : ptr+24;
            i_block(range) = dofs_E(tpl_I); 
            j_block(range) = dofs_N(tpl_J); 
            v_block(range) = Ce(:);         
            
            ptr = ptr + 24;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
end