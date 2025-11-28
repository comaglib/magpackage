classdef PostProcessor < handle
    % POSTPROCESSOR 有限元后处理模块 (v2.1 - Final Robust Fix)
    % 
    % 修复与功能整合:
    %   1. [Critical] 修复 parfor 循环范围错误：显式将 numElems 转换为 double 整数。
    %   2. [Optimization] 使用 parallel.pool.Constant 封装大数组，避免广播开销。
    %   3. [Feature] 保留 computeCentroids 以支持 Visualizer。
    %   4. [Math] computeMagnitude 支持复数场 (abs.^2)。
    
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
                % Triangulation requires double inputs
                T_double = double(obj.Mesh.T');
                P_double = double(obj.Mesh.P');
                obj.Triangulation = triangulation(T_double, P_double);
            end
        end
        
        function centers = computeCentroids(obj)
            % COMPUTECENTROIDS 计算所有单元的几何重心 (用于可视化)
            % 输出: centers [3 x NumElements]
            P = obj.Mesh.P;
            T = obj.Mesh.T;
            % 向量化计算
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
        end
        
        function B_elems = computeElementB(obj, A_sol, spaceName)
            % COMPUTEELEMENTB 计算所有单元重心的磁通密度 B
            % 支持并行计算和多列解向量(HBFEM)
            
            if nargin < 3
                keys = obj.DofHandler.DofMaps.keys;
                targetKey = keys{1}; 
                for k = 1:length(keys)
                    if contains(keys{k}, 'Nedelec'), targetKey = keys{k}; break; end
                end
                dofMap = obj.DofHandler.DofMaps(targetKey);
            else
                dofMap = obj.DofHandler.DofMaps(spaceName);
            end
            
            % [Critical Fix] 显式强制转换为标量 double 整数
            % 解决 "parfor 循环变量范围 'e' 的端点计算结果必须为整数" 错误
            numElems = fix(full(double(obj.Mesh.NumElements)));
            
            numCols = size(A_sol, 2);
            B_elems = zeros(3, numElems, numCols);
            
            % [Optimization] 使用 Constant 避免广播
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign));
            C_Sol = parallel.pool.Constant(A_sol);
            C_DofMap = parallel.pool.Constant(dofMap);
            
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center); 
            C_ref = reshape(curl_ref, 6, 3);
            
            parfor e = 1:numElems
                % 获取 Worker 本地数据副本
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Signs = C_Signs.Value;
                loc_Sol = C_Sol.Value; loc_DofMap = C_DofMap.Value;
                
                nodes = loc_T(:, e);
                p_elem = loc_P(:, nodes); 
                
                [J_mat, detJ] = compute_jacobian_tet(p_elem);
                invDetJ = 1.0 / detJ;
                C_phy = (C_ref * J_mat') * invDetJ;
                
                dofs = loc_DofMap(:, e);
                s = loc_Signs(:, e);
                A_local_all = loc_Sol(dofs, :) .* s;
                
                B_val = C_phy' * A_local_all;
                B_elems(:, e, :) = B_val;
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
            
            keys = obj.DofHandler.DofMaps.keys;
            dofMap = obj.DofHandler.DofMaps(keys{1});
            nodes = obj.Mesh.T(:, elem_idx);
            p_elem = obj.Mesh.P(:, nodes);
            [J_mat, detJ] = compute_jacobian_tet(p_elem);
            invDetJ = 1.0 / detJ;
            
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center);
            C_ref = reshape(curl_ref, 6, 3);
            C_phy = (C_ref * J_mat') * invDetJ;
            
            dofs = dofMap(:, elem_idx);
            s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            A_local = A_sol(dofs, :) .* s;
            B_val = C_phy' * A_local;
        end
        
        function mag = computeMagnitude(obj, VectorField)
            % [Fix] 支持复数场 (abs.^2)
            mag_sq = sum(abs(VectorField).^2, 1);
            mag = sqrt(mag_sq);
            [~, Ne, K] = size(VectorField);
            mag = reshape(mag, Ne, K);
        end
    end
end