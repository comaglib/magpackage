function K = assemble_magnetic_stiffness(Model, Nu_Override)
% ASSEMBLE_MAGNETIC_STIFFNESS 并行组装磁场刚度矩阵 (HPC重构版)
%
% 输入:
%   Model       - 模型结构体
%   Nu_Override - (可选) N_elems x 1 向量。如果提供，将忽略 Model.Materials
%                 直接使用此向量作为每个单元的磁阻率 nu。
%                 (用于非线性迭代更新)

    % fprintf('正在组装刚度矩阵 (并行模式)...\n');
    % t_start = tic;

    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    numEdges = size(Mesh.Edges, 2);
    
    % 1. 确定磁阻率向量 Nu_vec
    if nargin > 1 && ~isempty(Nu_Override)
        if length(Nu_Override) ~= numElems
            error('Nu_Override 维度 (%d) 与单元数 (%d) 不匹配', length(Nu_Override), numElems);
        end
        Nu_vec = Nu_Override;
        % fprintf('  - 使用外部传入的磁阻率向量\n');
    else
        % 从材料库读取 (标准流程)
        MatMap = Model.Materials.ActiveMap;
        MatLib = Model.Materials.Lib;
        Nu_vec = zeros(numElems, 1);
        
        % 预计算材料库的 nu
        % 假设 MatLib 很小，循环无所谓
        LibNu = zeros(length(MatLib), 1);
        mu0 = 4*pi*1e-7;
        for i = 1:length(MatLib)
            LibNu(i) = 1.0 / (MatLib(i).Mu_r * mu0);
        end
        
        % 向量化赋值
        Nu_vec = LibNu(MatMap);
    end

    % 2. 数据预分块 (Pre-chunking)
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    
    for k = 1:numChunks
        idx_s = (k-1)*chunkSize + 1;
        idx_e = min(k*chunkSize, numElems);
        indices = idx_s:idx_e;
        
        ElementBlocks{k}.T = T_raw(:, indices);
        ElementBlocks{k}.T2E = T2E_raw(:, indices);
        ElementBlocks{k}.Signs = Signs_raw(:, indices);
        ElementBlocks{k}.Nu = Nu_vec(indices); % 直接存入切片
        ElementBlocks{k}.Count = length(indices);
    end
    
    C_P = parallel.pool.Constant(Mesh.P);
    
    I_cell = cell(numChunks, 1);
    J_cell = cell(numChunks, 1);
    V_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T; loc_T2E = Block.T2E; loc_Signs = Block.Signs;
        loc_Nu = Block.Nu;
        loc_P_all = C_P.Value;
        
        i_block = zeros(36 * count, 1);
        j_block = zeros(36 * count, 1);
        v_block = zeros(36 * count, 1);
        ptr = 0;
        
        for e = 1:count
            nodes = loc_T(:, e);
            p_elem = loc_P_all(:, nodes);
            nu_val = loc_Nu(e);
            
            % Kernel
            Ke = Ke_curl_curl(p_elem, nu_val);
            
            s = loc_Signs(:, e);
            Ke_corr = Ke .* (s * s');
            
            rows = loc_T2E(:, e);
            R_idx = repmat(rows, 1, 6);
            C_idx = repmat(rows', 6, 1);
            
            range = ptr+1 : ptr+36;
            i_block(range) = R_idx;
            j_block(range) = C_idx;
            v_block(range) = Ke_corr;
            ptr = ptr + 36;
        end
        
        I_cell{k} = i_block;
        J_cell{k} = j_block;
        V_cell{k} = v_block;
    end
    
    I = vertcat(I_cell{:});
    J = vertcat(J_cell{:});
    V = vertcat(V_cell{:});
    
    K = sparse(I, J, V, numEdges, numEdges);
    
    % t_end = toc(t_start);
    % fprintf('刚度矩阵组装完成。耗时: %.4f 秒, 非零元: %d\n', t_end, nnz(K));
end