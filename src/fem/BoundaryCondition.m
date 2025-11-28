classdef BoundaryCondition < handle
    % BOUNDARYCONDITION 边界条件处理工具 (修正版 v2)
    %
    % 修复:
    % 1. [Critical] 边界识别算法从 "Edge-based" 改为 "Face-based"。
    %    旧算法错误地假设边界棱只属于1个单元，导致漏掉绝大多数边界，引发矩阵奇异。
    %    新算法正确识别属于边界面的所有棱。
    
    methods (Static)
        function is_bnd_dof = findOuterBoundaryDofs(mesh, dofHandler, space)
            % FINDOUTERBOUNDARYDOFS 识别外边界上的自由度 (PEC)
            % 逻辑:
            % 1. 找到只出现一次的面 (Boundary Faces)。
            % 2. 提取这些面上的所有棱 (Boundary Edges)。
            % 3. 将这些棱映射回全局 Edge ID。
            
            fprintf('[BC] Detecting boundary edges using Face Topology...\n');
            
            T = mesh.T; % [4 x Ne]
            num_elems = size(T, 2);
            
            % 1. 构建所有面 (每个四面体4个面)
            % Faces: [n1, n2, n3]
            F1 = T([1 2 3], :);
            F2 = T([1 2 4], :);
            F3 = T([1 3 4], :);
            F4 = T([2 3 4], :);
            
            % 拼接所有面 [3 x 4*Ne]
            AllFaces = [F1, F2, F3, F4];
            
            % 2. 对面内的节点排序 (保证唯一性)
            AllFaces = sort(AllFaces, 1);
            
            % 3. 统计面出现的次数
            % 转置为 [N_faces x 3] 以便使用 unique('rows')
            [U, ~, IC] = unique(AllFaces', 'rows');
            counts = accumarray(IC, 1);
            
            % 出现次数为 1 的面即为边界/表面
            bnd_face_mask = (counts == 1);
            BndFaces = U(bnd_face_mask, :); % [N_bnd_faces x 3]
            
            fprintf('     - Found %d boundary faces.\n', size(BndFaces, 1));
            
            if isempty(BndFaces)
                warning('No boundary faces detected! Check mesh connectivity.');
                is_bnd_dof = false(dofHandler.NumGlobalDofs, 1);
                return;
            end
            
            % 4. 提取边界上的棱
            % 每个三角形面有3条边: (1-2), (1-3), (2-3)
            % BndFaces 已经是排序过的 [n1 < n2 < n3]
            E1 = BndFaces(:, [1 2]);
            E2 = BndFaces(:, [1 3]);
            E3 = BndFaces(:, [2 3]);
            
            % 拼接并去重
            BndEdgesRaw = [E1; E2; E3];
            BndEdgesRaw = unique(BndEdgesRaw, 'rows'); % [N_bnd_edges x 2]
            
            % 5. 映射回全局 Edge ID (DoF)
            % mesh.Edges 是 [2 x N_global_edges] (已排序)
            % 我们需要找到 mesh.Edges 中哪些列出现在 BndEdgesRaw 中
            
            GlobalEdges = mesh.Edges'; % [N_global x 2]
            
            % 使用 ismember 进行行匹配
            is_bnd_dof = ismember(GlobalEdges, BndEdgesRaw, 'rows');
            
            % 6. 对齐输出维度
            % 如果由 DofHandler 管理的 DoF 比 Edge 多 (例如未来引入 Node DoF)，需要补零
            if length(is_bnd_dof) < dofHandler.NumGlobalDofs
                padding = false(dofHandler.NumGlobalDofs - length(is_bnd_dof), 1);
                is_bnd_dof = [is_bnd_dof; padding];
            end
            
            fprintf('[BC] Identified %d boundary edges (DoFs).\n', sum(is_bnd_dof));
        end
        
        function [K_sys, F_sys] = applyDirichlet(K, F, is_fixed)
            % APPLYDIRICHLET 施加齐次 Dirichlet 边界条件 (Masking Method)
            
            num_dofs = length(F);
            
            if ~islogical(is_fixed)
                idx = is_fixed;
                is_fixed = false(num_dofs, 1);
                is_fixed(idx) = true;
            end
            
            % 构建自由度掩码矩阵 M
            is_free = ~is_fixed;
            M_free = spdiags(double(is_free), 0, num_dofs, num_dofs);
            
            % 矩阵投影: K_sys = M*K*M + I_fixed
            K_sys = M_free * K * M_free;
            
            M_fixed = spdiags(double(is_fixed), 0, num_dofs, num_dofs);
            K_sys = K_sys + M_fixed;
            
            % 处理 RHS
            F_sys = F;
            F_sys(is_fixed) = 0;
            
            fprintf('[BC] Applied Dirichlet BCs to %d DoFs.\n', sum(is_fixed));
        end
    end
end