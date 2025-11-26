function MeshRaw = read_msh(filename)
% READ_MSH 读取 Gmsh .msh 文件 (项目内部专用格式修正版)
% 
% 修正:
% 1. read_nodes: 修改为读取 "ID X Y Z" 格式 (4列数据)，适配 create_delaunay_mesh_robust 的输出。
% 2. 保持了 read_elements 对交替格式的支持。
% 3. 这是一个高性能版本，使用 sscanf/fscanf 加速。

    if ~isfile(filename)
        error('文件未找到: %s', filename);
    end
    
    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));
    
    MeshRaw.P = [];
    MeshRaw.T = [];           % Tet4
    MeshRaw.RegionTags = [];  % Tet4 Tags
    MeshRaw.Faces = [];       % Tri3
    MeshRaw.FaceTags = [];    % Tri3 Tags
    
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
    % 内部函数: 读取节点 (修正版)
    % -------------------------------------------------------
    function read_nodes()
        header = sscanf(fgetl(fid), '%f'); 
        if isempty(header), return; end
        
        numBlocks = header(1);
        numNodes = header(2);
        
        P_all = zeros(3, numNodes);
        cnt = 0;
        
        for i = 1:numBlocks
            % 块头: dim tag parametric numNodesInBlock
            info = sscanf(fgetl(fid), '%f');
            if isempty(info), info = sscanf(fgetl(fid), '%f'); end
            
            nInBlock = info(4);
            
            % [修正]: 假设文件格式为每行 "ID X Y Z" (共4个数)
            % 之前的代码假设 Tags 和 Coords 是分开的 (Gmsh 4.1 标准)，这里改为适应生成器
            
            % 读取 4 * N 个浮点数
            data = fscanf(fid, '%f', 4 * nInBlock);
            
            % 重塑为 4 x N 矩阵
            % Row 1: ID
            % Row 2: X, Row 3: Y, Row 4: Z
            data = reshape(data, 4, nInBlock);
            
            % 提取坐标 (2:4 行)
            P_all(:, cnt+1 : cnt+nInBlock) = data(2:4, :);
            
            % 吃掉换行符
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
        
        totalTets = 0; totalTris = 0;
        
        for i = 1:numBlocks
            lineStr = fgetl(fid);
            info = sscanf(lineStr, '%f');
            if isempty(info), info = sscanf(fgetl(fid), '%f'); end
            
            entTag = info(2); 
            elemType = info(3);
            nInBlock = info(4);
            
            if elemType == 4 % Tet4 (ID + 4 Nodes = 5 cols)
                data = fscanf(fid, '%d', 5 * nInBlock);
                data = reshape(data, 5, nInBlock);
                
                T_list{end+1} = data(2:5, :); 
                Tag_list{end+1} = repmat(entTag, 1, nInBlock); 
                totalTets = totalTets + nInBlock;
                fgetl(fid);
                
            elseif elemType == 2 % Tri3 (ID + 3 Nodes = 4 cols)
                data = fscanf(fid, '%d', 4 * nInBlock);
                data = reshape(data, 4, nInBlock);
                
                F_list{end+1} = data(2:4, :); 
                FTag_list{end+1} = repmat(entTag, 1, nInBlock); 
                totalTris = totalTris + nInBlock;
                fgetl(fid);
                
            else
                % Skip unknown (estimate line length)
                nodesPerElem = 0;
                if elemType == 1, nodesPerElem = 2;
                elseif elemType == 15, nodesPerElem = 1;
                end
                
                if nodesPerElem > 0
                    fscanf(fid, '%d', (1 + nodesPerElem) * nInBlock);
                    fgetl(fid); 
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