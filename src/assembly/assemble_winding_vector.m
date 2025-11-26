function W = assemble_winding_vector(Model, Coil)
% ASSEMBLE_WINDING_VECTOR 组装绕组耦合向量
% 
% 计算磁链 Psi = W' * A
% 原理: 对线圈路径上的基函数 N 进行线积分。
%
% 输入: Model, Coil (含 P1, P2)
% 输出: W (N_edges x 1 稀疏向量)

    fprintf('正在组装绕组耦合向量 (Winding Vector)...\n');
    t_start = tic;

    Mesh = Model.Mesh;
    P = Mesh.P;
    T = Mesh.T;
    
    % 构建 triangulation 对象用于快速点定位
    TR = triangulation(T', P');
    
    % 自动计算积分步长
    vols = Mesh.Volumes;
    avg_len = mean(power(abs(vols), 1/3));
    step_size = avg_len * 0.2; 
    
    idx_list = [];
    val_list = [];
    
    numSegs = size(Coil.P1, 2);
    total_integrated_len = 0;
    
    % 循环线段
    for i = 1:numSegs
        p1 = Coil.P1(:, i);
        p2 = Coil.P2(:, i);
        L = norm(p2 - p1);
        
        if L < 1e-12, continue; end
        
        num_sub = ceil(L / step_size);
        t = linspace(0, 1, num_sub+1);
        
        % 中点积分
        mid_pts = p1 + (p2 - p1) * (t(1:end-1) + t(2:end))/2;
        dl = (p2 - p1) / num_sub;
        
        % 批量查找单元
        [elem_ids, ~] = pointLocation(TR, mid_pts');
        
        valid_mask = ~isnan(elem_ids);
        if ~any(valid_mask), continue; end
        
        valid_elems = elem_ids(valid_mask);
        valid_pts = mid_pts(:, valid_mask);
        num_valid = length(valid_elems);
        
        total_integrated_len = total_integrated_len + norm(dl) * num_valid;
        
        for k = 1:num_valid
            e = valid_elems(k);
            pt = valid_pts(:, k);
            
            nodes = T(:, e); 
            el_pts = P(:, nodes);
            
            % 几何导数
            v21=el_pts(:,2)-el_pts(:,1); 
            v31=el_pts(:,3)-el_pts(:,1); 
            v41=el_pts(:,4)-el_pts(:,1);
            J = [v21, v31, v41];
            invJ_T = inv(J)';
            G_phy = invJ_T * [-1 1 0 0; -1 0 1 0; -1 0 0 1];
            
            % 局部坐标
            uvw = J \ (pt - el_pts(:,1));
            L_vec = [1-sum(uvw); uvw];
            
            global_edge_idx = Mesh.T2E(:, e);
            signs = double(Mesh.T2E_Sign(:, e));
            edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            
            for edge_local = 1:6
                na = edge_pairs(edge_local, 1);
                nb = edge_pairs(edge_local, 2);
                
                % N = La * gLb - Lb * gLa
                N_vec = L_vec(na) * G_phy(:, nb) - L_vec(nb) * G_phy(:, na);
                
                % 投影: N dot dl
                val = dot(N_vec, dl) * signs(edge_local);
                
                if abs(val) > 1e-15
                    idx_list(end+1) = global_edge_idx(edge_local); %#ok<AGROW>
                    val_list(end+1) = val; %#ok<AGROW>
                end
            end
        end
    end
    
    numEdges = size(Mesh.Edges, 2);
    W = sparse(idx_list, 1, val_list, numEdges, 1);
    
    t_end = toc(t_start);
    coil_len = sum(sqrt(sum((Coil.P2-Coil.P1).^2,1)));
    fprintf('  - Winding assembled. NNZ=%d. IntegLen=%.4f m (Total Coil=%.4f m)\n', ...
        nnz(W), total_integrated_len, coil_len);
end