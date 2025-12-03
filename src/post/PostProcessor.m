classdef PostProcessor < handle
    % POSTPROCESSOR 有限元后处理模块 (v3.5 - Matrix Source Support)
    % 
    % 功能:
    %   1. 场量探测 (Probe): 支持 B 场、A 场、以及包含电位梯度的总电流密度 J。
    %   2. 能量计算: 磁场能量积分 (支持线性和非线性材料)。
    %   3. 节点平滑: 将单元值映射到节点以便绘图。
    %
    % 更新日志:
    %   v3.5: probeTotalCurrent 支持传入全场源电流矩阵 (3 x Ne)，
    %         可直接处理空间变化的电流分布（如自动生成的线圈电流场）。
    
    properties
        Assembler
        Mesh
        DofHandler
        Triangulation %用于快速查找点所在的单元
    end
    
    methods
        function obj = PostProcessor(assembler)
            obj.Assembler = assembler;
            obj.Mesh = assembler.Mesh;
            obj.DofHandler = assembler.DofHandler;
            
            % 初始化三角剖分对象用于快速点定位 (pointLocation)
            %这对于 probe 操作是必须的
            if ~isempty(obj.Mesh.T)
                T_double = double(obj.Mesh.T');
                P_double = double(obj.Mesh.P');
                obj.Triangulation = triangulation(T_double, P_double);
            end
        end
        
        function centers = computeCentroids(obj)
            % 计算所有单元的几何中心
            P = obj.Mesh.P; T = obj.Mesh.T;
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
        end
        
        % =================================================================
        %  核心场量计算 (批量处理)
        % =================================================================
        
        function B_elems = computeElementB(obj, A_sol, spaceName)
            % COMPUTEELEMENTB 计算每个单元中心的磁通密度 B = curl(A)
            % 输入: 
            %   A_sol - 磁矢位解向量
            %   spaceName - (可选) 有限元空间名称，默认取第一个
            
            if nargin < 3
                keys = obj.DofHandler.DofMaps.keys; 
                spaceName = keys{1}; 
            end
            
            dofMap = obj.DofHandler.DofMaps(spaceName);
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            numCols = size(A_sol, 2); % 支持多列解 (如多谐波或多时间步)
            B_elems = zeros(3, numElems, numCols);
            
            % 使用 Parallel Pool Constant 优化大数据传输
            % 避免在 parfor 中重复传递大数组
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign));
            C_Sol = parallel.pool.Constant(A_sol);
            C_DofMap = parallel.pool.Constant(dofMap);
            
            % 预计算参考单元的 Curl 值 (常数, 任意点即可)
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center); 
            C_ref = reshape(curl_ref, 6, 3);
            
            parfor e = 1:numElems
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Signs = C_Signs.Value;
                loc_Sol = C_Sol.Value; loc_DofMap = C_DofMap.Value;
                
                % 提取单元节点坐标
                nodes = loc_T(:, e); 
                p_elem = loc_P(:, nodes); 
                
                % 计算雅可比矩阵
                [J_mat, detJ] = compute_jacobian_tet(p_elem); 
                invDetJ = 1.0 / detJ;
                
                % Piola 变换: curl_phy = (1/detJ) * J * curl_ref
                C_phy = (C_ref * J_mat') * invDetJ;
                
                % 提取局部自由度系数
                dofs = loc_DofMap(:, e); 
                s = loc_Signs(:, e);
                A_local_all = loc_Sol(dofs, :) .* s;
                
                % 计算 B = Sum( Ni * ai )
                B_elems(:, e, :) = C_phy' * A_local_all;
            end
        end
        
        function Val_nodes = mapElementsToNodes(obj, Val_elems)
            % MAPELEMENTSTONODES 将单元中心值平均到节点上 (用于平滑云图显示)
            % 算法: 简单的体积平均或算术平均 (此处为算术平均)
            
            numNodes = obj.Mesh.NumNodes;
            numElems = obj.Mesh.NumElements;
            [dims, ~] = size(Val_elems);
            
            Val_nodes = zeros(dims, numNodes);
            T = obj.Mesh.T; % [4 x Ne]
            
            for d = 1:dims
                v = Val_elems(d, :); 
                
                % 将单元值广播到其4个节点
                node_idxs = T(:);
                vals_expanded = repmat(v, 4, 1);
                vals_expanded = vals_expanded(:);
                
                % 累加到节点
                sum_d = accumarray(node_idxs, vals_expanded, [numNodes, 1]);
                Val_nodes(d, :) = sum_d';
            end
            
            % 计算每个节点被多少个单元共享
            ones_expanded = ones(4 * numElems, 1);
            node_counts = accumarray(T(:), ones_expanded, [numNodes, 1]);
            
            % 平均
            node_counts(node_counts == 0) = 1; % 防止除零
            Val_nodes = Val_nodes ./ repmat(node_counts', dims, 1);
        end

        % =================================================================
        %  单点探测工具 (Probe)
        % =================================================================
        
        function [J_total, elem_idx] = probeTotalCurrent(obj, point, A_sol, V_sol, omega, sigmaMap, sourceMap)
            % PROBETOTALCURRENT 计算点处的总电流密度
            % 公式: J_total = J_source + J_eddy
            %              = J_source + sigma * (-j*omega*A - grad(V))
            
            % 1. 查找点所在的单元
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                J_total = [0;0;0];
                return;
            end
            
            tag = obj.Mesh.RegionTags(elem_idx);
            
            % 2. 获取源电流 (J_source)
            % [Update] 支持全场矩阵 (3 x Ne) 或 传统的 Map
            if isnumeric(sourceMap) && size(sourceMap, 2) == obj.Mesh.NumElements
                % 如果是全场矩阵，直接根据单元索引取值 (适用于复杂线圈)
                J_s = sourceMap(:, elem_idx);
            else
                % 否则按区域 Tag 查找 (适用于简单块状电流)
                J_s = obj.getPropertyValue(sourceMap, tag, [0,0,0]);
            end
            J_s = J_s(:);
            
            % 3. 获取导电率
            sigma = obj.getPropertyValue(sigmaMap, tag, 0);
            
            % 4. 计算涡流 J_eddy
            if abs(sigma) > 1e-12
                % 仅在导电区域计算电场项，避免在空气中查找 Lagrange 自由度导致越界
                [A_val, ~] = obj.internalProbeNedelec(A_sol, elem_idx, point, 'Nedelec_P1');
                
                if ~isempty(V_sol)
                    GradV_val = obj.internalProbeLagrangeGrad(V_sol, elem_idx, 'Lagrange_P1');
                else
                    GradV_val = [0;0;0];
                end
                
                E_eddy = -1j * omega * A_val - GradV_val;
                J_total = J_s + sigma * E_eddy;
            else
                % 非导电区域 (如空气、线圈源区域)，没有涡流
                J_total = J_s;
            end
        end
        
        function [B_val, elem_idx] = probeB(obj, A_sol, point)
            % PROBEB 计算点处的磁通密度 B = curl(A)
            
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                B_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            % 复用内部逻辑
            dofMap = obj.DofHandler.DofMaps('Nedelec_P1');
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            
            [J_mat, detJ] = compute_jacobian_tet(p_elem); 
            invDetJ = 1.0 / detJ;
            
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center);
            C_ref = reshape(curl_ref, 6, 3);
            
            % Piola Transform
            C_phy = (C_ref * J_mat') * invDetJ;
            
            dofs = dofMap(:, elem_idx); 
            s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            A_local = A_sol(dofs, :) .* s;
            
            B_val = C_phy' * A_local;
        end
        
        function [A_val, elem_idx] = probeSolution(obj, A_sol, point)
            % PROBESOLUTION 计算点处的磁矢位 A (Nedelec 插值)
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                A_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            [A_val, ~] = obj.internalProbeNedelec(A_sol, elem_idx, point, 'Nedelec_P1');
        end
        
        function result = probeLine(obj, A_sol, startPt, endPt, numPts, type, varargin)
            % PROBELINE 沿线段采样
            % type: 'B', 'A', 'J'
            pts = [linspace(startPt(1), endPt(1), numPts)', ...
                   linspace(startPt(2), endPt(2), numPts)', ...
                   linspace(startPt(3), endPt(3), numPts)'];
            
            result = zeros(numPts, 3); 
            
            for k = 1:numPts
                pt = pts(k, :);
                if strcmpi(type, 'B')
                    val = obj.probeB(A_sol, pt);
                elseif strcmpi(type, 'A')
                    val = obj.probeSolution(A_sol, pt);
                elseif strcmpi(type, 'J')
                    % 对于 J，varargin 需包含 V_sol, omega, sigmaMap, sourceMap
                    val = obj.probeTotalCurrent(pt, A_sol, varargin{:});
                else
                    error('Unknown probe type');
                end
                result(k, :) = val.'; 
            end
        end

        function mag = computeMagnitude(obj, VectorField)
            mag_sq = sum(abs(VectorField).^2, 1);
            mag = sqrt(mag_sq);
            [~, Ne, K] = size(VectorField);
            mag = reshape(mag, Ne, K);
        end
        
        % =================================================================
        %  能量计算
        % =================================================================
        
        function Energy = computeMagneticEnergy(obj, A_sol, MatLibData)
            % COMPUTEMAGNETICENERGY 计算全场磁能
            % 如果全是线性材料，使用快速矩阵乘法 E = 0.5 * A' * K * A
            % 如果包含非线性材料，使用逐单元数值积分
            
            all_linear = true;
            for i = 1:length(MatLibData)
                if ~strcmp(MatLibData(i).Type, 'Linear')
                    all_linear = false; break;
                end
            end
            space = FunctionSpace('Nedelec', 1); 
            
            if all_linear
                % 线性情况优化
                NuMap = zeros(length(MatLibData), 1);
                for i = 1:length(MatLibData)
                    NuMap(i) = MatLibData(i).Nu_Linear;
                end
                K = obj.Assembler.assembleStiffness(space, NuMap);
                Energy = 0.5 * real(A_sol' * K * A_sol);
            else
                % 非线性情况
                packedData = obj.Assembler.preparePackedData(space);
                packedData.RegionTags = obj.Mesh.RegionTags;
                Energy = obj.integrate_energy_kernel(packedData, A_sol, MatLibData);
            end
        end
    end
    
    % =================================================================
    %  私有辅助函数
    % =================================================================
    methods (Access = private)
        function [val, J_mat] = internalProbeNedelec(obj, solVec, elem_idx, point, spaceName)
            % 封装 Nedelec 插值逻辑
            dofMap = obj.DofHandler.DofMaps(spaceName);
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            
            [J_mat, ~] = compute_jacobian_tet(p_elem);
            
            % 物理坐标 -> 参考坐标
            x_ref = J_mat \ (point(:) - p_elem(:,1));
            
            [val_ref, ~] = nedelec_tet_p1(x_ref); 
            
            % H(curl) Piola 变换: u = J^{-T} * u_ref
            invJt = inv(J_mat)'; 
            N_phy = val_ref * invJt'; 
            
            dofs = dofMap(:, elem_idx); 
            s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            
            local_coeffs = solVec(dofs) .* s; 
            val = N_phy' * local_coeffs;
        end
        
        function grad_val = internalProbeLagrangeGrad(obj, solVec, elem_idx, spaceName)
            % 封装 Lagrange P1 梯度计算逻辑
            dofMap = obj.DofHandler.DofMaps(spaceName);
            
            % 处理全局索引到局部向量索引的转换
            global_dofs = dofMap(:, elem_idx);
            if obj.DofHandler.SpaceOffsets.isKey(spaceName)
                offset = obj.DofHandler.SpaceOffsets(spaceName);
            else
                offset = 0;
            end
            local_dofs = global_dofs - offset;
            
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            [J_mat, ~] = compute_jacobian_tet(p_elem);
            
            % P1 梯度参考值 (常数)
            q_ref = [0.25; 0.25; 0.25];
            [~, grad_ref] = lagrange_tet_p1(q_ref); % [4x3]
            
            % H(grad) Piola 变换: grad = J^{-T} * grad_ref
            invJt = inv(J_mat)';
            Grad_phy = grad_ref * invJt'; % [4x3]
            
            local_coeffs = solVec(local_dofs);
            grad_val = Grad_phy' * local_coeffs; % [3x1]
        end
        
        function val = getPropertyValue(~, mapOrArray, tag, defaultVal)
            % 统一处理 Map 或 Array 形式的属性查询
            if isa(mapOrArray, 'containers.Map')
                if mapOrArray.isKey(tag)
                    val = mapOrArray(tag);
                else
                    val = defaultVal;
                end
            elseif isnumeric(mapOrArray)
                if isscalar(mapOrArray)
                    val = mapOrArray; % 标量
                elseif numel(mapOrArray) == 3 && isvector(mapOrArray)
                     % 可能是 [Jx, Jy, Jz] 向量作为标量传入?
                     val = mapOrArray;
                elseif tag <= length(mapOrArray)
                    val = mapOrArray(tag);
                else
                    val = defaultVal;
                end
            else
                val = defaultVal;
            end
        end
        
        function total_W = integrate_energy_kernel(obj, PackedData, A_sol, MatLibData)
            % INTEGRATE_ENERGY_KERNEL 逐单元积分计算磁场能量
            % 并行计算加速
            
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
                
                [J_mat, detJ] = compute_jacobian_tet(p_elem); 
                invDetJ = 1.0 / detJ;
                
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
                        % 简化的共能/能量密度积分 (梯形公式近似)
                        % W = \int H dB
                        b_samples = [0, 0.5*B_mag, B_mag];
                        b_sq_samples = b_samples.^2;
                        nu_vals = zeros(1,3);
                        for k=1:3
                             [nu_k, ~] = MaterialLib.evaluate(b_sq_samples(k), matData);
                             nu_vals(k) = nu_k;
                        end
                        h_vals = nu_vals .* b_samples;
                        % Simpson / Trapezoidal rule
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