function [K_tan, Res] = assemble_nonlinear_jacobian(Model, x, b_ext)
% ASSEMBLE_NONLINEAR_JACOBIAN 组装非线性切线刚度矩阵和残差向量
%
% 输入:
%   Model  - 模型结构体
%   x      - 当前全局解向量 (N_edges x 1)
%   b_ext  - 外部载荷向量 (N_edges x 1) (即 RHS)
%
% 输出:
%   K_tan  - 切线刚度矩阵 (Jacobian)
%   Res    - 残差向量 (Internal Force - External Force)

    % fprintf('  [Assembly] Assembling Jacobian and Residual...\n');
    % t_start = tic;

    Mesh = Model.Mesh;
    
    % 1. 全局数据准备
    % ----------------------------------------------------
    % 将全局解 x 映射到每个单元的局部向量 (含符号修正)
    % ElemA: 6 x N_elems
    ElemA = prepare_element_A(Model, x);
    
    % 提取拓扑与材料映射
    T = Mesh.T;
    T2E = Mesh.T2E;
    Signs = double(Mesh.T2E_Sign);
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    % P 需要广播，封装
    C_P = parallel.pool.Constant(Mesh.P);
    
    numElems = size(T, 2);
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
        
        ElementBlocks{k}.T = T(:, indices);
        ElementBlocks{k}.T2E = T2E(:, indices);
        ElementBlocks{k}.Signs = Signs(:, indices);
        ElementBlocks{k}.MatIDs = MatMap(indices);
        ElementBlocks{k}.A_loc = ElemA(:, indices); % 关键: 只有切片数据传入
        ElementBlocks{k}.Count = length(indices);
    end
    
    % 结果缓存
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);     % For Jacobian values
    F_int_cell = cell(numChunks, 1); % For Internal Force values
    idx_F_cell = cell(numChunks, 1); % For Force indices
    
    % 广播材料库 (通常很小，直接拷入 Worker 也没问题)
    % 但为了严谨，我们假设 Lib 可以在 parfor 中只读访问
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        % 本地数据提取
        loc_T = Block.T;
        loc_T2E = Block.T2E;
        loc_Signs = Block.Signs;
        loc_MatIDs = Block.MatIDs;
        loc_A = Block.A_loc;
        
        loc_P_all = C_P.Value;
        
        % 预分配
        % Jacobian: 36 entries per elem
        % Residual: 6 entries per elem
        
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
            % 几何信息
            nodes = loc_T(:, e);
            p_elem = loc_P_all(:, nodes);
            
            % 当前单元解向量 (已修正符号)
            a_vec = loc_A(:, e);
            
            % 计算 B 场 (B = Curl N * a)
            % 使用 Kernel 生成的轻量函数
            B_vec = calc_element_B(p_elem, a_vec); % 3x1
            
            B_sq = sum(B_vec.^2);
            
            % 材料属性查询
            mat_id = loc_MatIDs(e);
            mat_info = MatLib(mat_id);
            
            [nu, dnu] = eval_material_nu(B_sq, mat_info);
            
            % ---------------------------------------------
            % 1. 计算残差贡献 (Internal Force)
            % F_int = K_secant(nu) * a
            % ---------------------------------------------
            
            % 调用各向同性内核计算 K_sec
            Ke_sec = Ke_curl_curl(p_elem, nu);
            
            f_local = Ke_sec * a_vec;
            
            % 符号修正 (仅用于组装回全局)
            s = loc_Signs(:, e);
            f_corrected = f_local .* s;
            
            % 存入 Buffer
            range_vec = ptr_vec+1 : ptr_vec+6;
            f_idx_block(range_vec) = loc_T2E(:, e);
            f_val_block(range_vec) = f_corrected;
            ptr_vec = ptr_vec + 6;
            
            % ---------------------------------------------
            % 2. 计算雅可比贡献 (Tangent Stiffness)
            % J = K_aniso(nu_tensor)
            % ---------------------------------------------
            
            if dnu == 0
                % 线性情况: Jacobian = Secant Stiffness
                Ke_tan = Ke_sec;
            else
                % 非线性情况: 构造微分磁阻率张量
                % nu_diff = nu*I + 2*dnu * B*B'
                % 注意: B*B' 是 3x3 矩阵
                
                % 手动构建 B*B' 以加速
                bx = B_vec(1); by = B_vec(2); bz = B_vec(3);
                BBt = [bx*bx, bx*by, bx*bz;
                       bx*by, by*by, by*bz;
                       bx*bz, by*bz, bz*bz];
                   
                Nu_tensor = nu * eye(3) + (2 * dnu) * BBt;
                
                % 提取唯一的对称分量输入给 Kernel
                % [xx, xy, xz, yy, yz, zz]
                nu_vec = [Nu_tensor(1,1), Nu_tensor(1,2), Nu_tensor(1,3), ...
                          Nu_tensor(2,2), Nu_tensor(2,3), Nu_tensor(3,3)];
                
                % 调用各向异性内核
                Ke_tan = Ke_aniso(p_elem, nu_vec);
            end
            
            % 符号修正 (行和列)
            % K_corrected = K .* (s * s')
            Ke_tan_corr = Ke_tan .* (s * s');
            
            % 存入 Buffer (Matrix)
            row_indices = loc_T2E(:, e);
            
            % 快速展开索引
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
    % Jacobian Matrix
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    K_tan = sparse(I, J, V, numEdges, numEdges);
    
    % Internal Force Vector
    idx_F = vertcat(idx_F_cell{:});
    val_F = vertcat(F_int_cell{:});
    F_int = sparse(idx_F, 1, val_F, numEdges, 1); % sparse 自动累加
    % 转为全向量 (如果维度不大) 或者保持 sparse
    F_int = full(F_int);
    
    % 4. 计算残差 Res = F_int - F_ext
    % ----------------------------------------------------
    Res = F_int - b_ext;
    
    % t_end = toc(t_start);
    % fprintf('  [Assembly] Done in %.4f s. Res Norm: %.4e\n', t_end, norm(Res));
end