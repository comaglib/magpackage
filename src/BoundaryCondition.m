classdef BoundaryCondition < handle
    % BOUNDARYCONDITION 边界条件处理工具 (v3.1 - EntityMap Support)
    %
    % 更新:
    %   1. findOuterBoundaryDofs 现在利用 DofHandler 的 EntityMap 进行查询。
    %      这对于 Reduced Domain (部分区域) 自由度分配至关重要。
    
    methods (Static)
        function is_bnd_dof = findOuterBoundaryDofs(mesh, dofHandler, space)
            % FINDOUTERBOUNDARYDOFS 自动识别外边界上的自由度
            
            fprintf('[BC] Detecting boundary entities for %s...\n', space.toString());
            
            % 1. 识别边界"面" (Boundary Faces)
            T = mesh.T; 
            F1 = T([1 2 3], :); F2 = T([1 2 4], :);
            F3 = T([1 3 4], :); F4 = T([2 3 4], :);
            
            AllFaces = sort([F1, F2, F3, F4], 1);
            [U, ~, IC] = unique(AllFaces', 'rows');
            counts = accumarray(IC, 1);
            BndFaces = U(counts == 1, :); 
            
            if isempty(BndFaces)
                warning('No boundary faces detected!');
                is_bnd_dof = false(dofHandler.NumGlobalDofs, 1);
                return;
            end
            
            % 2. 获取该 Space 的实体映射表 (Entity ID -> Global DoF)
            spaceKey = space.toString();
            if ~dofHandler.EntityMaps.isKey(spaceKey)
                 error('DofHandler EntityMap missing for %s. Did you call distributeDofs?', spaceKey);
            end
            entityMap = dofHandler.EntityMaps(spaceKey);
            
            % 3. 提取边界实体并查找其全局 DoF
            bnd_global_dofs = [];
            
            switch lower(space.Type)
                case 'nedelec'
                    % --- Nedelec: 边界棱 ---
                    E1 = BndFaces(:, [1 2]); E2 = BndFaces(:, [1 3]); E3 = BndFaces(:, [2 3]);
                    BndEdgesRaw = unique([E1; E2; E3], 'rows');
                    
                    % 全局 Edge 列表
                    GlobalEdges = mesh.Edges'; 
                    
                    % 找到边界 Edge 的 ID (行号)
                    [is_on_bnd, ~] = ismember(GlobalEdges, BndEdgesRaw, 'rows');
                    bnd_edge_ids = find(is_on_bnd);
                    
                    % 查表获取全局 DoF
                    if ~isempty(bnd_edge_ids)
                        bnd_global_dofs = entityMap(bnd_edge_ids);
                    end
                    
                case 'lagrange'
                    % --- Lagrange: 边界节点 ---
                    BndNodes = unique(BndFaces(:));
                    
                    % 查表获取全局 DoF
                    if ~isempty(BndNodes)
                        bnd_global_dofs = entityMap(BndNodes);
                    end
                    
                otherwise
                    error('Unsupported space type.');
            end
            
            % 4. 过滤无效 DoF (0) 并构建掩码
            % 在 Reduced Domain 模式下，某些边界节点可能未分配自由度 (值为0)
            valid_dofs = bnd_global_dofs(bnd_global_dofs > 0);
            
            is_bnd_dof = false(dofHandler.NumGlobalDofs, 1);
            is_bnd_dof(valid_dofs) = true;
            
            fprintf('[BC] Identified %d boundary DoFs for %s.\n', sum(is_bnd_dof), spaceKey);
        end
        
        function [K_sys, F_sys] = applyDirichlet(K, F, is_fixed)
            % APPLYDIRICHLET (标准实现)
            num_dofs = size(K, 1);
            
            if ~islogical(is_fixed)
                idx = is_fixed;
                is_fixed = false(num_dofs, 1);
                valid_mask = (idx >= 1 & idx <= num_dofs);
                is_fixed(idx(valid_mask)) = true;
            else
                len_mask = length(is_fixed);
                if len_mask < num_dofs
                    padding = false(num_dofs - len_mask, 1);
                    is_fixed = [is_fixed(:); padding];
                elseif len_mask > num_dofs
                     error('Mask too long.');
                end
            end
            
            is_free = ~is_fixed;
            M_free = spdiags(double(is_free), 0, num_dofs, num_dofs);
            K_sys = M_free * K * M_free;
            M_fixed = spdiags(double(is_fixed), 0, num_dofs, num_dofs);
            K_sys = K_sys + M_fixed;
            
            F_sys = F;
            if isvector(F)
                F_sys(is_fixed) = 0;
            else
                F_sys(is_fixed, :) = 0;
            end
        end
    end
end