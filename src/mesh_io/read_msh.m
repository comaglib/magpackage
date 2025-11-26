function MeshRaw = read_msh(filename)
% READ_MSH 读取 Gmsh .msh 文件 (支持 Format 4.1) - 高性能版
% 优化: 将 str2num 替换为 sscanf，显著提升头部解析速度
% 修正: 保持了对换行符的健壮处理

    if ~isfile(filename)
        error('文件未找到: %s', filename);
    end
    
    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));
    
    MeshRaw.P = [];
    MeshRaw.T = [];           % Tet4 单元
    MeshRaw.RegionTags = [];  % Tet4 物理标签
    MeshRaw.Faces = [];       % Tri3 单元 (边界)
    MeshRaw.FaceTags = [];    % Tri3 物理标签 (边界)
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        
        if strcmp(line, '$Nodes')
            read_nodes();
        elseif strcmp(line, '$Elements')
            read_elements();
        end
    end
    
    % -------------------------------------------------------
    % 内部函数: 读取节点
    % -------------------------------------------------------
    function read_nodes()
        % 使用 sscanf 替代 str2num 解析 "numBlocks numNodes ..."
        header = sscanf(fgetl(fid), '%f'); 
        if isempty(header), return; end
        
        numBlocks = header(1);
        numNodes = header(2);
        
        P_all = zeros(3, numNodes);
        cnt = 0;
        
        for i = 1:numBlocks
            % 解析块信息: entityDim entityTag parametric numNodesInBlock
            info = sscanf(fgetl(fid), '%f');
            if isempty(info), info = sscanf(fgetl(fid), '%f'); end % 跳过可能的空行
            
            nInBlock = info(4);
            
            % 读取并跳过 Node Tags
            textscan(fid, '%d', nInBlock); 
            
            % 读取坐标
            data = fscanf(fid, '%f', 3 * nInBlock);
            P_all(:, cnt+1 : cnt+nInBlock) = reshape(data, 3, nInBlock);
            
            % 吃掉行末的换行符
            fgetl(fid); 
            
            cnt = cnt + nInBlock;
        end
        MeshRaw.P = P_all;
        fprintf('  - 读取节点数: %d\n', numNodes);
    end

    % -------------------------------------------------------
    % 内部函数: 读取单元
    % -------------------------------------------------------
    function read_elements()
        header = sscanf(fgetl(fid), '%f');
        if isempty(header), return; end
        numBlocks = header(1);
        
        T_list = {}; Tag_list = {}; 
        F_list = {}; FTag_list = {}; 
        
        totalTets = 0;
        totalTris = 0;
        
        for i = 1:numBlocks
            % 读取块信息头
            lineStr = fgetl(fid);
            info = sscanf(lineStr, '%f');
            
            % 如果读到了空行，再读一次
            if isempty(info)
                info = sscanf(fgetl(fid), '%f'); 
            end
            
            entTag = info(2); 
            elemType = info(3);
            nInBlock = info(4);
            
            if elemType == 4 % 4-node Tetrahedron
                data = fscanf(fid, '%d', 5 * nInBlock);
                data = reshape(data, 5, nInBlock);
                T_list{end+1} = data(2:5, :); %#ok<AGROW>
                Tag_list{end+1} = repmat(entTag, 1, nInBlock); %#ok<AGROW>
                totalTets = totalTets + nInBlock;
                
                fgetl(fid); % 吃掉换行符
                
            elseif elemType == 2 % 3-node Triangle
                data = fscanf(fid, '%d', 4 * nInBlock);
                data = reshape(data, 4, nInBlock);
                F_list{end+1} = data(2:4, :); %#ok<AGROW>
                FTag_list{end+1} = repmat(entTag, 1, nInBlock); %#ok<AGROW>
                totalTris = totalTris + nInBlock;
                
                fgetl(fid); % 吃掉换行符
                
            else
                % 跳过不支持的单元类型
                nodesPerElem = 0;
                if elemType == 1, nodesPerElem = 2;
                elseif elemType == 15, nodesPerElem = 1;
                end
                
                if nodesPerElem > 0
                    fscanf(fid, '%d', (1 + nodesPerElem) * nInBlock);
                    fgetl(fid); % 吃掉换行符
                end
            end
        end
        
        if ~isempty(T_list)
            MeshRaw.T = [T_list{:}];
            MeshRaw.RegionTags = [Tag_list{:}];
        end
        if ~isempty(F_list)
            MeshRaw.Faces = [F_list{:}];
            MeshRaw.FaceTags = [FTag_list{:}];
        end
        fprintf('  - 读取单元: Tet=%d, Tri=%d\n', totalTets, totalTris);
    end
end