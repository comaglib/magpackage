classdef DofHandler < handle
    % DOFHANDLER 自由度管理器
    % 负责将 FunctionSpace 映射到 Mesh 上的具体索引
    
    properties
        Mesh          % 关联的 Mesh 对象
        NumGlobalDofs % 全局自由度总数
        DofMaps       % 存储每个 Space 的 CellDofs 映射表
    end
    
    methods
        function obj = DofHandler(mesh)
            obj.Mesh = mesh;
            obj.NumGlobalDofs = 0;
            obj.DofMaps = containers.Map();
        end
        
        function cellDofs = distributeDofs(obj, space)
            % DISTRIBUTEDOFS 为指定的 FunctionSpace 分配自由度
            % 返回: [N_dof_per_elem x N_elems] 的索引矩阵
            
            fprintf('[DofHandler] Distributing DoFs for %s...\n', space.toString());
            
            switch lower(space.Type)
                case 'lagrange'
                    if space.Order == 1
                        % Lagrange P1: DoF 直接对应 Mesh 的节点
                        % CellDofs 就是 Mesh.T (Local Node IDs -> Global Node IDs)
                        cellDofs = obj.Mesh.T;
                        
                    else
                        error('Only Lagrange P1 is supported in this version.');
                    end
                    
                case 'nedelec'
                    if space.Order == 1
                        % Nedelec P1: DoF 对应 Mesh 的棱 (Edges)
                        
                        % 1. 确保 Mesh 已经生成了棱信息
                        if isempty(obj.Mesh.Edges)
                            obj.Mesh.generateEdges();
                        end
                        
                        % 2. CellDofs 就是 Mesh.T2E
                        cellDofs = obj.Mesh.T2E;
                        
                    else
                        error('Only Nedelec P1 is supported in this version.');
                    end
                    
                otherwise
                    error('Unknown FunctionSpace type: %s', space.Type);
            end
            
            % 存储映射供后续查询
            key = space.toString();
            obj.DofMaps(key) = cellDofs;
            
            % [核心修正] 自动更新全局自由度总数
            % 这确保了 Assembler 在组装稀疏矩阵时知道矩阵的大小
            current_max_dof = max(cellDofs(:));
            if current_max_dof > obj.NumGlobalDofs
                obj.NumGlobalDofs = current_max_dof;
            end
            
            fprintf('             -> Generated map size: [%d x %d], Max DoF: %d\n', ...
                size(cellDofs, 1), size(cellDofs, 2), obj.NumGlobalDofs);
        end
        
        function packedData = packForKernel(obj, space)
            % PACKFORKERNEL 准备传入 parallel kernel 的纯数据结构
            % 符合 DOD 设计原则
            
            key = space.toString();
            if ~obj.DofMaps.isKey(key)
                error('DoFs for %s have not been distributed yet.', key);
            end
            
            packedData.CellDofs = obj.DofMaps(key);
            
            % 如果是 Nedelec，还需要额外的方向修正信息 (Sign)
            if strcmpi(space.Type, 'nedelec')
                packedData.Signs = obj.Mesh.T2E_Sign;
            end
        end
    end
end