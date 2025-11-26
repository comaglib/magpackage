function MeshRaw = read_mphtxt(filename, unit_tag)
% READ_MPHTXT 读取 COMSOL Multiphysics 网格文件 (.mphtxt) (含单位转换)
%
% 输入:
%   filename - 文件路径
%   unit_tag - (可选) 长度单位字符串: 'm'(默认), 'mm', 'cm', 'um', 'nm', 'in', 'ft'
%              代码会自动将其转换为标准国际单位制 (米)。
%
% 示例:
%   Mesh = read_mphtxt('motor.mphtxt', 'mm');

    if nargin < 2, unit_tag = 'm'; end

    % 1. 确定单位缩放因子
    scale_factor = 1.0;
    switch lower(unit_tag)
        case 'm',  scale_factor = 1.0;
        case 'cm', scale_factor = 1e-2;
        case 'mm', scale_factor = 1e-3;
        case 'um', scale_factor = 1e-6;
        case 'nm', scale_factor = 1e-9;
        case 'in', scale_factor = 0.0254;
        case 'ft', scale_factor = 0.3048;
        otherwise
            warning('未知的单位标签: "%s"。默认按米 (scale=1.0) 处理。', unit_tag);
    end
    
    if scale_factor ~= 1.0
        fprintf('  [Unit] 将坐标从 "%s" 转换为 "m" (Scale: %.2e)\n', unit_tag, scale_factor);
    end

    if ~isfile(filename)
        error('文件未找到: %s', filename);
    end
    
    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));
    
    fprintf('正在读取 COMSOL 网格: %s ...\n', filename);
    
    MeshRaw.P = [];
    MeshRaw.T = [];
    MeshRaw.RegionTags = [];
    MeshRaw.Faces = [];
    MeshRaw.FaceTags = [];
    
    num_nodes_expected = 0;
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        
        % --- 节点数量 ---
        if contains(line, '# number of mesh vertices')
            C = textscan(line, '%d');
            if ~isempty(C) && ~isempty(C{1})
                num_nodes_expected = C{1};
                fprintf('  - 检测到节点数量: %d\n', num_nodes_expected);
            end
        end
        
        % --- 节点坐标 ---
        if contains(line, '# Mesh vertex coordinates')
            if num_nodes_expected == 0
                warning('未找到节点数量定义，尝试继续读取。');
            end
            fprintf('  - 读取节点坐标...\n');
            
            % 假设 3D
            data = fscanf(fid, '%f', 3 * num_nodes_expected);
            
            if length(data) ~= 3 * num_nodes_expected
                % 尝试回退并报错 (或支持2D逻辑)
                error('读取到的坐标数据量 (%d) 与预期 (%d) 不符。', length(data), 3*num_nodes_expected);
            end
            
            MeshRaw.P = reshape(data, 3, num_nodes_expected);
            
            % [核心] 应用单位缩放
            if scale_factor ~= 1.0
                MeshRaw.P = MeshRaw.P * scale_factor;
            end
            
            fprintf('    已存储 %d 个节点。\n', size(MeshRaw.P, 2));
        end
        
        % --- 单元块识别 ---
        if contains(line, '# type name')
            if contains(line, ' tet ')
                fprintf('  - 发现四面体块 (tet)...\n');
                [elems, tags] = parse_element_block_v2(fid);
                if ~isempty(elems)
                    MeshRaw.T = [MeshRaw.T, elems];
                    MeshRaw.RegionTags = [MeshRaw.RegionTags, tags];
                end
                
            elseif contains(line, ' tri ')
                fprintf('  - 发现三角形块 (tri)...\n');
                [elems, tags] = parse_element_block_v2(fid);
                if ~isempty(elems)
                    MeshRaw.Faces = [MeshRaw.Faces, elems];
                    MeshRaw.FaceTags = [MeshRaw.FaceTags, tags];
                end
            end
        end
    end
    
    % 子函数: 解析特定单元块 (保持不变)
    function [elems, tags] = parse_element_block_v2(fid)
        elems = []; tags = [];
        num_elems = 0;
        nodes_per_elem = 0;
        
        found_elems = false;
        found_tags = false;
        
        while ~feof(fid)
            pos = ftell(fid);
            line_sub = strtrim(fgetl(fid));
            if isempty(line_sub), continue; end
            
            if contains(line_sub, '# type name')
                fseek(fid, pos, 'bof'); return; 
            end
            
            if contains(line_sub, '# number of') && contains(line_sub, 'per element')
                C = textscan(line_sub, '%d');
                if ~isempty(C) && ~isempty(C{1}), nodes_per_elem = C{1}; end
            end
            
            if contains(line_sub, '# number of elements')
                C = textscan(line_sub, '%d');
                if ~isempty(C) && ~isempty(C{1}), num_elems = C{1}; end
            end
            
            if contains(line_sub, '# Elements') && ~found_elems
                if num_elems == 0 || nodes_per_elem == 0
                    warning('读取拓扑前未找到数量定义。'); return;
                end
                raw_conn = fscanf(fid, '%d', nodes_per_elem * num_elems);
                if length(raw_conn) == nodes_per_elem * num_elems
                    elems = reshape(raw_conn, nodes_per_elem, num_elems) + 1; 
                    found_elems = true;
                    fprintf('    已读取拓扑: %d 个单元\n', num_elems);
                end
            end
            
            if contains(line_sub, '# Geometric entity indices') && ~found_tags
                if num_elems == 0, return; end
                raw_tags = fscanf(fid, '%d', num_elems);
                if length(raw_tags) == num_elems
                    tags = reshape(raw_tags, 1, num_elems);
                    found_tags = true;
                    return; 
                end
            end
        end
    end
end