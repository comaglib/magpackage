function MeshRaw = read_msh(filename)
% READ_MSH 读取 Gmsh .msh 文件 (鲁棒增强版)
% 
% 修复:
% 1. 增强了对空行和文件指针位置的处理，防止读取错位。
% 2. 增加了对节点和单元数据块的维度检查。

    if ~isfile(filename)
        error('文件未找到: %s', filename);
    end
    
    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));
    
    MeshRaw.P = []; 
    MeshRaw.T = []; 
    MeshRaw.RegionTags = [];
    MeshRaw.Faces = []; 
    MeshRaw.FaceTags = [];
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        
        if strcmp(line, '$Nodes')
            % 读取块头: numEntityBlocks numNodes minNodeTag maxNodeTag
            info = sscanf_safe(fid);
            if isempty(info), continue; end
            numBlocks = info(1);
            numNodes = info(2);
            
            P_all = zeros(3, numNodes);
            cnt = 0;
            
            for k = 1:numBlocks
                % 块信息: dim tag parametric numNodes
                b_info = sscanf_safe(fid);
                nInBlock = b_info(4);
                
                % 读取数据: ID X Y Z (共4列)
                % 使用 %f 读取所有数据，包括 ID
                data = fscanf(fid, '%f', 4 * nInBlock);
                
                if length(data) ~= 4 * nInBlock
                    error('读取节点数据失败: 预期 %d 个数据，实际读取 %d 个。', 4 * nInBlock, length(data));
                end
                
                data = reshape(data, 4, nInBlock);
                
                % 提取坐标 (2:4 行)
                P_all(:, cnt+1:cnt+nInBlock) = data(2:4, :);
                cnt = cnt + nInBlock;
            end
            MeshRaw.P = P_all;
            
        elseif strcmp(line, '$Elements')
            % 读取块头: numEntityBlocks numElements minTag maxTag
            info = sscanf_safe(fid);
            if isempty(info), continue; end
            numBlocks = info(1);
            
            T_list = {}; 
            Tag_list = {}; 
            F_list = {}; 
            FTag_list = {};
            
            for k = 1:numBlocks
                % 块信息: dim tag type numElements
                b_info = sscanf_safe(fid);
                entTag = b_info(2);
                elemType = b_info(3);
                nInBlock = b_info(4);
                
                if elemType == 2 % Triangle (3-node) -> 4 cols (ID n1 n2 n3)
                    data = fscanf(fid, '%f', 4 * nInBlock); % 使用 %f 防止溢出，稍后转整
                    if length(data) ~= 4 * nInBlock
                         error('读取三角形单元失败');
                    end
                    data = reshape(data, 4, nInBlock);
                    F_list{end+1} = data(2:4, :);
                    FTag_list{end+1} = repmat(entTag, 1, nInBlock);
                    
                elseif elemType == 4 % Tetrahedron (4-node) -> 5 cols (ID n1 n2 n3 n4)
                    data = fscanf(fid, '%f', 5 * nInBlock);
                    if length(data) ~= 5 * nInBlock
                         error('读取四面体单元失败');
                    end
                    data = reshape(data, 5, nInBlock);
                    T_list{end+1} = data(2:5, :);
                    Tag_list{end+1} = repmat(entTag, 1, nInBlock);
                    
                else
                    % 跳过未知类型 (Line=1, Point=15)
                    % 需要根据类型确定列数
                    n_skip = 0;
                    if elemType == 1, n_skip = 3; % ID n1 n2
                    elseif elemType == 15, n_skip = 2; % ID n1
                    end
                    
                    if n_skip > 0
                        fscanf(fid, '%f', n_skip * nInBlock);
                    else
                         % 如果遇到完全未知的类型，可能会导致读取错位
                         warning('遇到未知单元类型 %d，尝试跳过但可能失败。', elemType);
                    end
                end
            end
            
            if ~isempty(T_list)
                MeshRaw.T = cat(2, T_list{:});
                MeshRaw.RegionTags = cat(2, Tag_list{:});
            end
            if ~isempty(F_list)
                MeshRaw.Faces = cat(2, F_list{:});
                MeshRaw.FaceTags = cat(2, FTag_list{:});
            end
        end
    end
    
    % 内部辅助函数: 安全读取一行数字，跳过空行
    function val = sscanf_safe(fid)
        val = [];
        while isempty(val) && ~feof(fid)
            l = strtrim(fgetl(fid));
            if isempty(l), continue; end
            val = sscanf(l, '%f');
        end
    end
end