classdef PostProcessor < handle
    % POSTPROCESSOR 有限元后处理模块 (v3.1 - Field Probe Support)
    % 
    % 更新:
    %   1. 新增 probeSolution 方法，用于直接插值向量场 A。
    %      这对于计算涡流密度 J = -j*w*sigma*A 至关重要。
    
    properties
        Assembler
        Mesh
        DofHandler
        Triangulation 
    end
    
    methods
        function obj = PostProcessor(assembler)
            obj.Assembler = assembler;
            obj.Mesh = assembler.Mesh;
            obj.DofHandler = assembler.DofHandler;
            
            if ~isempty(obj.Mesh.T)
                T_double = double(obj.Mesh.T');
                P_double = double(obj.Mesh.P');
                obj.Triangulation = triangulation(T_double, P_double);
            end
        end
        
        function centers = computeCentroids(obj)
            P = obj.Mesh.P; T = obj.Mesh.T;
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
        end
        
        function B_elems = computeElementB(obj, A_sol, spaceName)
            if nargin < 3
                keys = obj.DofHandler.DofMaps.keys; targetKey = keys{1}; 
                dofMap = obj.DofHandler.DofMaps(targetKey);
            else
                dofMap = obj.DofHandler.DofMaps(spaceName);
            end
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            numCols = size(A_sol, 2);
            B_elems = zeros(3, numElems, numCols);
            
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign));
            C_Sol = parallel.pool.Constant(A_sol);
            C_DofMap = parallel.pool.Constant(dofMap);
            
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center); 
            C_ref = reshape(curl_ref, 6, 3);
            
            parfor e = 1:numElems
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Signs = C_Signs.Value;
                loc_Sol = C_Sol.Value; loc_DofMap = C_DofMap.Value;
                
                nodes = loc_T(:, e); p_elem = loc_P(:, nodes); 
                [J_mat, detJ] = compute_jacobian_tet(p_elem); invDetJ = 1.0 / detJ;
                C_phy = (C_ref * J_mat') * invDetJ;
                
                dofs = loc_DofMap(:, e); s = loc_Signs(:, e);
                A_local_all = loc_Sol(dofs, :) .* s;
                B_elems(:, e, :) = C_phy' * A_local_all;
            end
        end
        
        function [B_val, elem_idx] = probeB(obj, A_sol, point)
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                warning('Point is outside the mesh.');
                B_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            keys = obj.DofHandler.DofMaps.keys; dofMap = obj.DofHandler.DofMaps(keys{1});
            nodes = obj.Mesh.T(:, elem_idx); p_elem = obj.Mesh.P(:, nodes);
            [J_mat, detJ] = compute_jacobian_tet(p_elem); invDetJ = 1.0 / detJ;
            q_center = [0.25; 0.25; 0.25]; [~, curl_ref] = nedelec_tet_p1(q_center);
            C_ref = reshape(curl_ref, 6, 3); C_phy = (C_ref * J_mat') * invDetJ;
            
            dofs = dofMap(:, elem_idx); s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            A_local = A_sol(dofs, :) .* s;
            B_val = C_phy' * A_local;
        end
        
        function [A_val, elem_idx] = probeSolution(obj, A_sol, point)
            % PROBESOLUTION 探测磁矢位 A 的值
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                A_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            keys = obj.DofHandler.DofMaps.keys; dofMap = obj.DofHandler.DofMaps(keys{1});
            nodes = obj.Mesh.T(:, elem_idx); p_elem = obj.Mesh.P(:, nodes);
            
            % 计算局部坐标和基函数
            [J_mat, ~] = compute_jacobian_tet(p_elem);
            x_ref = J_mat \ (point(:) - p_elem(:,1));
            [val_ref, ~] = nedelec_tet_p1(x_ref); 
            
            % Piola 变换 (Value): N = J^{-T} * N_hat
            invJt = inv(J_mat)';
            N_phy = val_ref * invJt'; 
            
            dofs = dofMap(:, elem_idx); s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            A_local = A_sol(dofs, :) .* s; 
            A_val = N_phy' * A_local;
        end
        
        function [J_val, elem_idx] = probeEddyCurrent(obj, A_sol, point, omega, sigmaMap)
            % PROBEEDDYCURRENT 封装电流密度提取功能
            % J = -j * w * sigma * A
            
            % 1. 获取 A
            [A_val, elem_idx] = obj.probeSolution(A_sol, point);
            
            if isnan(elem_idx)
                J_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            % 2. 获取 Sigma
            % sigmaMap 是 [MaxTag x 1] 的数组
            tag = obj.Mesh.RegionTags(elem_idx);
            if tag > length(sigmaMap) || tag < 1
                sigma = 0;
            else
                sigma = sigmaMap(tag);
            end
            
            % 3. 计算 J
            J_val = -1j * omega * sigma * A_val;
        end
        
        function mag = computeMagnitude(obj, VectorField)
            mag_sq = sum(abs(VectorField).^2, 1);
            mag = sqrt(mag_sq);
            [~, Ne, K] = size(VectorField);
            mag = reshape(mag, Ne, K);
        end
        
        function Energy = computeMagneticEnergy(obj, A_sol, MatLibData)
            all_linear = true;
            for i = 1:length(MatLibData)
                if ~strcmp(MatLibData(i).Type, 'Linear')
                    all_linear = false; break;
                end
            end
            keys = obj.DofHandler.DofMaps.keys;
            space = FunctionSpace('Nedelec', 1); 
            
            if all_linear
                NuMap = zeros(length(MatLibData), 1);
                for i = 1:length(MatLibData)
                    NuMap(i) = MatLibData(i).Nu_Linear;
                end
                K = obj.Assembler.assembleStiffness(space, NuMap);
                Energy = 0.5 * A_sol' * K * A_sol;
            else
                packedData = obj.Assembler.preparePackedData(space);
                packedData.RegionTags = obj.Mesh.RegionTags;
                Energy = obj.integrate_energy_kernel(packedData, A_sol, MatLibData);
            end
        end
        
        function total_W = integrate_energy_kernel(obj, PackedData, A_sol, MatLibData)
            C_P = parallel.pool.Constant(PackedData.P);
            C_T = parallel.pool.Constant(PackedData.T);
            C_Dofs = parallel.pool.Constant(PackedData.CellDofs);
            C_Signs = parallel.pool.Constant(PackedData.Signs);
            C_Tags = parallel.pool.Constant(PackedData.RegionTags);
            C_Sol = parallel.pool.Constant(A_sol);
            C_MatLib = parallel.pool.Constant(MatLibData);
            
            [q_pts, q_w] = get_quadrature_data('tet', 2);
            [~, curl_ref] = nedelec_tet_p1(q_pts);
            C_Qw = parallel.pool.Constant(q_w);
            C_CurlRef = parallel.pool.Constant(curl_ref);
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            n_q = length(q_w);
            
            W_local = zeros(numElems, 1);
            
            parfor e = 1:numElems
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Sol = C_Sol.Value;
                loc_MatLib = C_MatLib.Value; loc_Tags = C_Tags.Value;
                loc_Signs = C_Signs.Value; loc_Dofs = C_Dofs.Value;
                loc_qw = C_Qw.Value; loc_curl_ref = C_CurlRef.Value;
                
                nodes = loc_T(:, e); p_elem = loc_P(:, nodes);
                dofs = loc_Dofs(:, e); s = loc_Signs(:, e);
                tag = loc_Tags(e); matData = loc_MatLib(tag);
                A_local = loc_Sol(dofs) .* double(s);
                
                [J_mat, detJ] = compute_jacobian_tet(p_elem); invDetJ = 1.0 / detJ;
                
                w_sum = 0;
                for q = 1:n_q
                    weight = loc_qw(q) * abs(detJ);
                    C_ref = reshape(loc_curl_ref(:, :, q), 6, 3);
                    C_phy = (C_ref * J_mat') * invDetJ;
                    B_vec = C_phy' * A_local;
                    B_mag = norm(B_vec);
                    B_sq = B_mag^2;
                    
                    w_density = 0;
                    if strcmp(matData.Type, 'Linear')
                        w_density = 0.5 * matData.Nu_Linear * B_sq;
                    else
                        b_samples = [0, 0.5*B_mag, B_mag];
                        b_sq_samples = b_samples.^2;
                        nu_vals = zeros(1,3);
                        for k=1:3
                             [nu_k, ~] = MaterialLib.evaluate(b_sq_samples(k), matData);
                             nu_vals(k) = nu_k;
                        end
                        h_vals = nu_vals .* b_samples;
                        w_density = (B_mag / 6.0) * (0 + 4*h_vals(2) + h_vals(3));
                    end
                    w_sum = w_sum + weight * w_density;
                end
                W_local(e) = w_sum;
            end
            total_W = sum(W_local);
        end
    end
end