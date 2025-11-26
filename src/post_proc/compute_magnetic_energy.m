function TotalEnergy = compute_magnetic_energy(Model, A_curr, Coil)
% COMPUTE_MAGNETIC_ENERGY 计算全域磁场储能 (支持复数/实数)
%
% 公式:
%   AC (Peak): W_avg = 1/4 * int( B . H* ) dV
%   Transient: W     = 1/2 * int( B . H ) dV
%
% 注意: 假设材料局部线性 (H = nu * B)

    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    
    % 判断模式
    is_complex = ~isreal(A_curr) || ~isreal(Coil.I);
    if is_complex
        EnergyFactor = 0.25;
        % fprintf('  [Post] Calculating AC Magnetic Energy (Average)...\n');
    else
        EnergyFactor = 0.5;
        % fprintf('  [Post] Calculating Transient Magnetic Energy...\n');
    end
    
    % 预分块
    T_raw = Mesh.T; T2E_raw = Mesh.T2E; Signs_raw = double(Mesh.T2E_Sign);
    MatMap_raw = Model.Materials.ActiveMap;
    
    chunkSize = 5000;
    numChunks = ceil(numElems / chunkSize);
    ElementBlocks = cell(numChunks, 1);
    
    for k = 1:numChunks
        s_idx = (k-1)*chunkSize + 1;
        e_idx = min(k*chunkSize, numElems);
        idx = s_idx:e_idx;
        ElementBlocks{k}.T = T_raw(:, idx);
        ElementBlocks{k}.T2E = T2E_raw(:, idx);
        ElementBlocks{k}.Signs = Signs_raw(:, idx);
        ElementBlocks{k}.MatIDs = MatMap_raw(idx);
        ElementBlocks{k}.Count = length(idx);
    end
    
    C_P = parallel.pool.Constant(Mesh.P);
    C_A = parallel.pool.Constant(A_curr);
    C_Coil = parallel.pool.Constant(Coil);
    C_MatLib = parallel.pool.Constant(Model.Materials.Lib);
    
    Energy_cell = cell(numChunks, 1);
    
    parfor k = 1:numChunks
        Block = ElementBlocks{k};
        count = Block.Count;
        
        loc_T = Block.T; loc_T2E = Block.T2E; loc_Signs = Block.Signs;
        loc_MatIDs = Block.MatIDs;
        
        loc_P_all = C_P.Value;
        loc_A = C_A.Value;
        loc_Coil = C_Coil.Value;
        loc_MatLib = C_MatLib.Value;
        
        local_energy = 0;
        
        % 重心计算
        Centers = zeros(3, count);
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            Centers(:, e) = mean(pts, 2);
        end
        
        % 使用串行版 Bs 计算器
        Bs_all = compute_biot_savart_B_serial(loc_Coil, Centers);
        
        for e = 1:count
            nodes = loc_T(:, e);
            pts = loc_P_all(:, nodes);
            
            edges = loc_T2E(:, e);
            s = loc_Signs(:, e);
            a_vec = loc_A(edges) .* s;
            
            Br = calc_element_B(pts, a_vec);
            B_tot = Br + Bs_all(:, e);
            
            % |B|^2 = B . conj(B)
            B_sq_val = sum(B_tot .* conj(B_tot));
            B_sq_real = real(B_sq_val); % 能量必须是实数
            
            mat_id = loc_MatIDs(e);
            mat_info = loc_MatLib(mat_id);
            
            % 这里仅处理线性 nu (对于非线性，后处理通常只取近似)
            % 使用 eval_material_nu 兼容两者
            [nu, ~] = eval_material_nu(B_sq_real, mat_info);
            
            v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
            cp = cross(v21, v31);
            Vol = abs(dot(cp, v41)) / 6.0;
            
            local_energy = local_energy + EnergyFactor * nu * B_sq_real * Vol;
        end
        Energy_cell{k} = local_energy;
    end
    
    TotalEnergy = sum([Energy_cell{:}]);
end