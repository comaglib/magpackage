classdef PostProcessor < handle
    % POSTPROCESSOR 有限元后处理模块 (v1.4 - Broadcast Optimized)
    % 
    % 修复:
    %   1. [Performance] 使用 parallel.pool.Constant 封装 P, T, A_sol 等大数组，
    %      彻底解决 parfor 中的广播变量警告，显著降低内存开销。
    %   2. [Robustness] 修正了 computeMagnitude 中的维度重塑逻辑。
    %   3. [Safety] 保持了 numElems 的整数强制转换。
    
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
        
        function B_elems = computeElementB(obj, A_sol, spaceName)
            % COMPUTEELEMENTB 计算所有单元重心的磁通密度 B
            % 优化: 解决了 P 和 A_sol 的广播变量问题
            
            if nargin < 3
                keys = obj.DofHandler.DofMaps.keys;
                targetKey = '';
                for k = 1:length(keys)
                    if contains(keys{k}, 'Nedelec')
                        targetKey = keys{k}; break;
                    end
                end
                if isempty(targetKey), error('No Nedelec space found.'); end
                dofMap = obj.DofHandler.DofMaps(targetKey);
            else
                dofMap = obj.DofHandler.DofMaps(spaceName);
            end
            
            % [Fix] 显式转换为标量 Double 整数
            numElems = fix(full(double(obj.Mesh.NumElements)));
            
            numCols = size(A_sol, 2);
            B_elems = zeros(3, numElems, numCols);
            
            % [Optimization] 将大数组放入 Constant Pool，避免广播
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign));
            C_Sol = parallel.pool.Constant(A_sol);
            C_DofMap = parallel.pool.Constant(dofMap);
            
            % 预计算基函数 (P1 Nedelec Curl 在单元内为常数)
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center); 
            C_ref = reshape(curl_ref, 6, 3);
            
            parfor e = 1:numElems
                % 获取本地引用 (Worker 端零拷贝访问)
                loc_P = C_P.Value;
                loc_T = C_T.Value;
                loc_Signs = C_Signs.Value;
                loc_Sol = C_Sol.Value;
                loc_DofMap = C_DofMap.Value;
                
                % 1. 几何计算 (此时 loc_P 的索引访问不再触发广播警告)
                nodes = loc_T(:, e);
                p_elem = loc_P(:, nodes); 
                
                [J_mat, detJ] = compute_jacobian_tet(p_elem);
                invDetJ = 1.0 / detJ;
                
                % 2. Piola 变换: Curl_phy = (1/detJ) * J * Curl_ref
                C_phy = (C_ref * J_mat') * invDetJ;
                
                % 3. 获取系数
                dofs = loc_DofMap(:, e);
                s = loc_Signs(:, e);
                
                % A_local_all: [6 x K]
                A_local_all = loc_Sol(dofs, :) .* s;
                
                % 4. 计算 B = Curl * A
                B_val = C_phy' * A_local_all;
                
                % 5. 存入切片
                B_elems(:, e, :) = B_val;
            end
        end
        
        function [B_val, elem_idx] = probeB(obj, A_sol, point)
            % PROBEB 在指定空间点探测 B 场
            if isempty(obj.Triangulation)
                error('Triangulation not initialized.');
            end
            
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
            % 计算向量场模值 B_mag [Ne x K]
            % VectorField 输入维度: [3 x Ne x K]
            
            mag_sq = sum(VectorField.^2, 1);
            mag = sqrt(mag_sq);
            
            % [Fix] 健壮的 reshape 逻辑
            [~, Ne, K] = size(VectorField);
            mag = reshape(mag, Ne, K);
        end
    end
end