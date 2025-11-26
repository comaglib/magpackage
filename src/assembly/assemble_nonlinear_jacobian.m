function [K_tan, Res] = assemble_nonlinear_jacobian(Model, x, b_ext)
% ASSEMBLE_NONLINEAR_JACOBIAN 组装非线性切线刚度矩阵和残差向量 (HPC优化版)
%
% 输入:
%   Model  - 模型结构体
%   x      - 当前全局解向量 (N_edges x 1)
%   b_ext  - 外部载荷向量 (N_edges x 1) (即 RHS)
%
% 输出:
%   K_tan  - 切线刚度矩阵 (Jacobian)
%   Res    - 残差向量 (Internal Force - External Force)
%
% 优化:
%   1. 使用 Cell 预分块消除 T, T2E, Signs, ElemA 的广播。
%   2. 使用 parallel.pool.Constant 封装 MatLib 和 P。

    % fprintf('  [Assembly] Assembling Jacobian and Residual...\n');
    
    Mesh = Model.Mesh;
    
    % 1. 全局数据准备
    % ----------------------------------------------------
    % 将全局解 x 映射到每个单元的局部向量 (含符号修正)
    % ElemA: 6 x N_elems
    ElemA = prepare_element_A(Model, x);
    
    % 提取拓扑与材料映射
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    MatMap_raw = Model.Materials.ActiveMap;
    
    % [优化点] 封装广播变量
    C_P = parallel.pool.Constant(Mesh.P);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    
    numElems = size(T_raw, 2);
    numEdges = size(Mesh.Edges, 2);
    
    % 2. 数据预分块 (Pre-chunking)
    % ----------------------------------------------------
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    
    ElementBlocks = cell(numChunks, 1);
    
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, numElems);
        indices = idx_s:idx_e;
        
        % 仅切片存储必要数据
        ElementBlocks{k}.T = T_raw(:, indices);
        ElementBlocks{k}.T2E = T2E_raw(:, indices);
        ElementBlocks{k}.Signs = Signs_raw(:, indices);
        ElementBlocks{k}.MatIDs = MatMap_raw(indices);
        ElementBlocks{k}.A_loc = ElemA(:, indices); 
        ElementBlocks{k}.Count = length(indices);
    end
    
    % 结果缓存
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);     % For Jacobian values
    F_int_cell = cell(numChunks, 1); % For Internal Force values
    idx_F_cell = cell(numChunks, 1); % For Force indices
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        % 本地数据提取
        loc_T = Block.T;
        loc_T2E = Block.T2E;
        loc_Signs = Block.Signs;
        loc_MatIDs = Block.MatIDs;
        loc_A = Block.A_loc;
        
        % 获取常量
        loc_P_all = C_P.Value;
        loc_MatLib = C_MatLib.Value; % [修复] 在 Worker 端获取本地副本
        
        % 预分配
        i_block = zeros(36 * count, 1);
        j_block = zeros(36 * count, 1);
        v_block = zeros(36 * count, 1);
        
        f_val_block = zeros(6 * count, 1);
        f_idx_block = zeros(6 * count, 1);
        
        ptr_mat = 0;
        ptr_vec = 0;
        
        % -------------------------------------------------
        % 单元循环
        % -------------------------------------------------
        for e = 1:count
            nodes = loc_T(:, e);
            p_elem = loc_P_all(:, nodes);
            
            a_vec = loc_A(:, e);
            
            % 计算 B 场
            B_vec = calc_element_B(p_elem, a_vec); 
            B_sq = sum(B_vec.^2);
            
            % [修复] 使用本地材料库副本
            mat_id = loc_MatIDs(e);
            mat_info = loc_MatLib(mat_id);
            
            [nu, dnu] = eval_material_nu(B_sq, mat_info);
            
            % --- 1. 计算残差 (Internal Force) ---
            % F_int = K_secant(nu) * a
            Ke_sec = Ke_curl_curl(p_elem, nu);
            f_local = Ke_sec * a_vec;
            
            s = loc_Signs(:, e);
            f_corrected = f_local .* s;
            
            range_vec = ptr_vec+1 : ptr_vec+6;
            f_idx_block(range_vec) = loc_T2E(:, e);
            f_val_block(range_vec) = f_corrected;
            ptr_vec = ptr_vec + 6;
            
            % --- 2. 计算雅可比 (Tangent Stiffness) ---
            if dnu == 0
                Ke_tan = Ke_sec;
            else
                bx = B_vec(1); by = B_vec(2); bz = B_vec(3);
                BBt = [bx*bx, bx*by, bx*bz;
                       bx*by, by*by, by*bz;
                       bx*bz, by*bz, bz*bz];
                   
                % 微分磁阻率张量
                % nu_diff = nu*I + 2*dnu * B*B'
                Nu_tensor = nu * eye(3) + (2 * dnu) * BBt;
                
                nu_vec = [Nu_tensor(1,1), Nu_tensor(1,2), Nu_tensor(1,3), ...
                          Nu_tensor(2,2), Nu_tensor(2,3), Nu_tensor(3,3)];
                
                Ke_tan = Ke_aniso(p_elem, nu_vec);
            end
            
            % 符号修正 K_corrected = K .* (s * s')
            Ke_tan_corr = Ke_tan .* (s * s');
            
            row_indices = loc_T2E(:, e);
            R_idx = repmat(row_indices, 1, 6);
            C_idx = repmat(row_indices', 6, 1);
            
            range_mat = ptr_mat+1 : ptr_mat+36;
            i_block(range_mat) = R_idx;
            j_block(range_mat) = C_idx;
            v_block(range_mat) = Ke_tan_corr;
            ptr_mat = ptr_mat + 36;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
        
        F_int_cell{k} = f_val_block;
        idx_F_cell{k} = f_idx_block;
    end
    
    % 3. 全局组装 (Reduce)
    % ----------------------------------------------------
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    K_tan = sparse(I, J, V, numEdges, numEdges);
    
    idx_F = vertcat(idx_F_cell{:});
    val_F = vertcat(F_int_cell{:});
    F_int = sparse(idx_F, 1, val_F, numEdges, 1);
    F_int = full(F_int);
    
    % 4. 计算残差 Res = F_int - F_ext
    Res = F_int - b_ext;
end