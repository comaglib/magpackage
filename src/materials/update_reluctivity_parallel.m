function [Nu_new, max_B] = update_reluctivity_parallel(Model, A_phasor, Coil)
% UPDATE_RELUCTIVITY_PARALLEL 并行更新所有单元的磁阻率
%
% 功能:
%   根据当前的磁矢量位 A 和线圈电流，计算每个单元的 B 场模值，
%   并查询 B-H 曲线（通过 eval_material_nu）得到新的磁阻率 nu。
%
% 输入:
%   Model    - 模型结构体 (含 Mesh, Materials)
%   A_phasor - 当前 A 场 (N_edges x 1)
%   Coil     - 线圈对象 (含当前电流 I)
%
% 输出:
%   Nu_new   - 更新后的磁阻率向量 (N_elems x 1)
%   max_B    - 全场最大磁通密度 (Scalar)

    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    
    % 1. 数据预分块 (Pre-chunking)
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    MatMap_raw = Model.Materials.ActiveMap;
    
    for k = 1:numChunks
        idx = (k-1)*chunkSize + 1 : min(k*chunkSize, numElems);
        ElementBlocks{k}.T = T_raw(:, idx);
        ElementBlocks{k}.T2E = T2E_raw(:, idx);
        ElementBlocks{k}.Signs = Signs_raw(:, idx);
        ElementBlocks{k}.MatIDs = MatMap_raw(idx);
        ElementBlocks{k}.Count = length(idx);
    end
    
    % 2. 广播常量
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(A_phasor);
    C_Coil = parallel.pool.Constant(Coil);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    
    Nu_cell = cell(numChunks, 1);
    B_cell = cell(numChunks, 1);
    
    % 3. 并行计算
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        cnt = Block.Count;
        
        loc_P = C_P.Value;
        loc_A = C_A.Value;
        loc_Coil = C_Coil.Value;
        loc_MatLib = C_MatLib.Value;
        
        nu_local = zeros(cnt, 1);
        b_max_local = 0;
        
        % 批量计算单元重心
        Centers = zeros(3, cnt);
        for e = 1:cnt
            nodes = Block.T(:, e);
            Centers(:, e) = mean(loc_P(:, nodes), 2);
        end
        
        % 计算源场 Bs (串行版函数，避免嵌套并行)
        Bs_all = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e = 1:cnt
            nodes = Block.T(:, e);
            pts = loc_P(:, nodes);
            edges = Block.T2E(:, e);
            s = Block.Signs(:, e);
            
            % 计算反应场 Br
            a_vec = loc_A(edges) .* s;
            Br = calc_element_B(pts, a_vec);
            
            % 总场 B
            B_tot = Br + Bs_all(:, e);
            
            % B_mag (有效值或峰值，取决于定义)
            B_mag = norm([abs(B_tot(1)), abs(B_tot(2)), abs(B_tot(3))]);
            
            if B_mag > b_max_local
                b_max_local = B_mag;
            end
            
            % 查表更新 Nu
            mat_id = Block.MatIDs(e);
            mat_info = loc_MatLib(mat_id);
            
            [nu, ~] = eval_material_nu(B_mag^2, mat_info);
            nu_local(e) = nu;
        end
        
        Nu_cell{k} = nu_local;
        B_cell{k} = b_max_local;
    end
    
    % 4. 结果合并
    Nu_new = vertcat(Nu_cell{:});
    max_B = max([B_cell{:}]);
end