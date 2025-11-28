classdef BoundaryCondition < handle
    % BOUNDARYCONDITION 边界条件处理工具 (v2.1 - MRHS Support)
    %
    % 更新:
    %   1. applyDirichlet 现在支持矩阵形式的 F [Ndofs x Nrhs]。
    %      这对参数提取和并行激励非常重要。
    
    methods (Static)
        function is_bnd_dof = findOuterBoundaryDofs(mesh, dofHandler, space)
            % FINDOUTERBOUNDARYDOFS 识别外边界上的自由度 (PEC)
            % 逻辑: 基于面拓扑 (Face Topology)
            
            fprintf('[BC] Detecting boundary edges using Face Topology...\n');
            
            T = mesh.T; 
            % 1. 构建所有面
            F1 = T([1 2 3], :);
            F2 = T([1 2 4], :);
            F3 = T([1 3 4], :);
            F4 = T([2 3 4], :);
            
            AllFaces = [F1, F2, F3, F4];
            AllFaces = sort(AllFaces, 1);
            
            % 2. 统计面出现的次数
            [U, ~, IC] = unique(AllFaces', 'rows');
            counts = accumarray(IC, 1);
            
            % 3. 提取边界上的棱
            BndFaces = U(counts == 1, :); 
            
            if isempty(BndFaces)
                warning('No boundary faces detected! Check mesh connectivity.');
                is_bnd_dof = false(dofHandler.NumGlobalDofs, 1);
                return;
            end
            
            E1 = BndFaces(:, [1 2]);
            E2 = BndFaces(:, [1 3]);
            E3 = BndFaces(:, [2 3]);
            
            BndEdgesRaw = unique([E1; E2; E3], 'rows');
            
            % 4. 映射回全局 Edge ID
            GlobalEdges = mesh.Edges';
            is_bnd_dof = ismember(GlobalEdges, BndEdgesRaw, 'rows');
            
            % 5. 对齐维度
            if length(is_bnd_dof) < dofHandler.NumGlobalDofs
                padding = false(dofHandler.NumGlobalDofs - length(is_bnd_dof), 1);
                is_bnd_dof = [is_bnd_dof; padding];
            end
            
            fprintf('[BC] Identified %d boundary edges (DoFs).\n', sum(is_bnd_dof));
        end
        
        function [K_sys, F_sys] = applyDirichlet(K, F, is_fixed)
            % APPLYDIRICHLET 施加 Dirichlet 边界条件
            % 
            % 输入:
            %   K - 刚度矩阵
            %   F - 载荷向量 或 载荷矩阵 [Ndofs x Nrhs]
            %   is_fixed - 固定自由度的逻辑掩码
            
            num_dofs = size(K, 1);
            
            if ~islogical(is_fixed)
                idx = is_fixed;
                is_fixed = false(num_dofs, 1);
                is_fixed(idx) = true;
            end
            
            % 1. 处理刚度矩阵 K (置1法: K_free = M*K*M + I_fixed)
            is_free = ~is_fixed;
            M_free = spdiags(double(is_free), 0, num_dofs, num_dofs);
            K_sys = M_free * K * M_free;
            
            M_fixed = spdiags(double(is_fixed), 0, num_dofs, num_dofs);
            K_sys = K_sys + M_fixed;
            
            % 2. 处理载荷 F
            % 如果 F 是矩阵，需要将所有列中对应固定行的值置零
            F_sys = F;
            
            if isvector(F)
                % 单向量情况
                F_sys(is_fixed) = 0;
            else
                % [Updated] 多向量情况 (MRHS)
                F_sys(is_fixed, :) = 0;
            end
            
            % fprintf('[BC] Applied Dirichlet BCs to %d DoFs.\n', sum(is_fixed));
        end
    end
end