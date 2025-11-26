function TotalEnergy = compute_magnetic_energy(Model, A_curr, Coil)
% COMPUTE_MAGNETIC_ENERGY 计算全域磁场储能 (HPC修正版)
% 
% W = sum (0.5 * nu * |B_total|^2 * Vol)
% B_total = Curl Ar + Bs
% 
% 修正: 使用 compute_biot_savart_B_serial 避免 parfor 嵌套错误

    % fprintf('  [Post] Calculating Magnetic Energy...\n');
    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    
    T_raw = Mesh.T;
    T2E_raw = Mesh.T2E;
    Signs_raw = double(Mesh.T2E_Sign);
    MatMap_raw = Model.Materials.ActiveMap;
    
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numElems);
        indices = s_idx:e_idx;
        
        ElementBlocks{k}.T = T_raw(:, indices);
        ElementBlocks{k}.T2E = T2E_raw(:, indices);
        ElementBlocks{k}.Signs = Signs_raw(:, indices);
        ElementBlocks{k}.MatIDs = MatMap_raw(indices);
        ElementBlocks{k}.Count = length(indices);
    end
    
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(A_curr);
    C_Coil = parallel.pool.Constant(Coil);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib); 
    
    Energy_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T;
        loc_T2E = Block.T2E;
        loc_Signs = Block.Signs;
        loc_MatIDs = Block.MatIDs;
        
        loc_P_all = C_P.Value;
        loc_A = C_A.Value;
        loc_Coil = C_Coil.Value;
        loc_MatLib = C_MatLib.Value;
        
        local_energy = 0;
        
        Centers = zeros(3, count);
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            Centers(:, e) = mean(pts, 2);
        end
        
        % [核心修正] 调用串行版
        Bs_all = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            
            edges = loc_T2E(:, e);
            s = loc_Signs(:, e);
            a_vec = loc_A(edges) .* s;
            
            Br = calc_element_B(pts, a_vec);
            
            B_tot = Br + Bs_all(:, e);
            B_sq = sum(B_tot.^2);
            
            mat_id = loc_MatIDs(e);
            mat_info = loc_MatLib(mat_id);
            [nu, ~] = eval_material_nu(B_sq, mat_info);
            
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            cp = [v21(2)*v31(3)-v21(3)*v31(2);
                  v21(3)*v31(1)-v21(1)*v31(3);
                  v21(1)*v31(2)-v21(2)*v31(1)];
            Vol = abs(sum(cp .* v41)) / 6.0;
            
            local_energy = local_energy + 0.5 * nu * B_sq * Vol;
        end
        
        Energy_cell{k} = local_energy;
    end
    
    TotalEnergy = sum([Energy_cell{:}]);
end