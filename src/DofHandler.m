classdef DofHandler < handle
    % DOFHANDLER 自由度管理器 (v3.2 - Safe Indexing)
    
    properties
        Mesh
        NumGlobalDofs
        DofMaps
        SpaceOffsets
        SpaceLocalSizes
        EntityMaps
    end
    
    methods
        function obj = DofHandler(mesh)
            obj.Mesh = mesh;
            obj.NumGlobalDofs = 0;
            obj.DofMaps = containers.Map();
            obj.SpaceOffsets = containers.Map();
            obj.SpaceLocalSizes = containers.Map();
            obj.EntityMaps = containers.Map();
        end
        
        function [cellDofs, offset] = distributeDofs(obj, space, targetRegionTags)
            if nargin < 3, targetRegionTags = []; end
            
            key = space.toString();
            if obj.DofMaps.isKey(key)
                cellDofs = obj.DofMaps(key);
                offset = obj.SpaceOffsets(key);
                return;
            end
            
            % 1. 确定活跃单元
            numElems = obj.Mesh.NumElements;
            if numElems == 0, error('Mesh has no elements!'); end
            
            if isempty(targetRegionTags)
                isActiveElem = true(1, numElems);
            else
                elementTags = obj.Mesh.RegionTags;
                if isempty(elementTags)
                    warning('Mesh has no RegionTags. Assuming Full Domain.');
                    isActiveElem = true(1, numElems);
                else
                    isActiveElem = ismember(elementTags, targetRegionTags);
                end
            end
            
            % 2. 分配逻辑
            currentOffset = obj.NumGlobalDofs;
            
            switch lower(space.Type)
                case 'lagrange' 
                    if space.Order ~= 1, error('Only Lagrange P1 supported.'); end
                    numEntities = obj.Mesh.NumNodes;
                    Topology = obj.Mesh.T; 
                    
                    activeElemsIdx = find(isActiveElem);
                    if isempty(activeElemsIdx)
                        activeNodes = [];
                    else
                        activeNodes = unique(Topology(:, activeElemsIdx));
                    end
                    
                    numActiveDoFs = length(activeNodes);
                    entityMap = zeros(numEntities, 1);
                    if numActiveDoFs > 0
                        entityMap(activeNodes) = currentOffset + (1:numActiveDoFs)';
                    end
                    cellDofs = entityMap(Topology);
                    
                case 'nedelec'
                    if space.Order ~= 1, error('Only Nedelec P1 supported.'); end
                    
                    if isempty(obj.Mesh.Edges) || isempty(obj.Mesh.T2E)
                        fprintf('   [DofHandler] Generating Edges...\n');
                        obj.Mesh.generateEdges();
                    end
                    
                    numEntities = size(obj.Mesh.Edges, 2);
                    Topology = obj.Mesh.T2E; 
                    
                    % 检查 T2E 是否有效 (不应全为 0)
                    if max(Topology(:)) == 0
                        warning('DofHandler:InvalidTopology', 'Mesh.T2E contains all zeros. Edges might not be generated correctly.');
                    end
                    
                    activeElemsIdx = find(isActiveElem);
                    if isempty(activeElemsIdx)
                        activeEdges = [];
                    else
                        activeEdges = unique(Topology(:, activeElemsIdx));
                        
                        % [FIX] 安全移除 0 索引 (非边)
                        if ~isempty(activeEdges) && activeEdges(1) == 0
                            activeEdges(1) = [];
                        end
                    end
                    
                    numActiveDoFs = length(activeEdges);
                    entityMap = zeros(numEntities, 1);
                    if numActiveDoFs > 0
                        entityMap(activeEdges) = currentOffset + (1:numActiveDoFs)';
                    end
                    cellDofs = entityMap(Topology);
                    
                otherwise
                    error('Unknown space: %s', space.Type);
            end
            
            % 3. 存储
            obj.DofMaps(key) = cellDofs;
            obj.SpaceOffsets(key) = currentOffset;
            obj.SpaceLocalSizes(key) = numActiveDoFs;
            obj.EntityMaps(key) = entityMap;
            
            obj.NumGlobalDofs = obj.NumGlobalDofs + numActiveDoFs;
            offset = currentOffset;
            
            fprintf('[DofHandler] %s: Active DoFs %d / %d. Range [%d - %d]\n', ...
                key, numActiveDoFs, numEntities, currentOffset+1, obj.NumGlobalDofs);
        end
        
        function packedData = packForKernel(obj, space)
            key = space.toString();
            packedData.CellDofs = obj.DofMaps(key);
            if strcmpi(space.Type, 'nedelec')
                packedData.Signs = obj.Mesh.T2E_Sign;
            end
        end
        
        function globalIndices = getGlobalIndices(obj, space, localIndices)
            key = space.toString();
            if ~obj.EntityMaps.isKey(key), error('Not distributed'); end
            entityMap = obj.EntityMaps(key);
            globalIndices = entityMap(localIndices);
        end
    end
end