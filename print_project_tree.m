function print_project_tree(startPath)
% PRINT_PROJECT_TREE 递归打印项目目录结构树
%
% 用法:
%   print_project_tree()            % 打印当前工作目录的结构
%   print_project_tree('src')       % 打印 src 文件夹的结构
%   print_project_tree('D:\MyProj') % 打印指定绝对路径的结构

    if nargin < 1
        startPath = pwd; % 默认为当前目录
    end

    if ~exist(startPath, 'dir')
        error('错误: 找不到目录 "%s"', startPath);
    end

    % 获取根目录名称用于显示
    [~, rootName] = fileparts(startPath);
    if isempty(rootName), rootName = startPath; end
    
    fprintf('%s/\n', rootName);
    
    % 开始递归遍历
    traverse_dir(startPath, '');
end

function traverse_dir(currentPath, prefix)
    % 获取当前目录下的所有文件和文件夹
    items = dir(currentPath);
    
    % --- 过滤列表 ---
    % 1. 移除 '.' 和 '..'
    items = items(~ismember({items.name}, {'.', '..'}));
    
    % 2. 移除以 '.' 开头的隐藏文件/文件夹 (如 .git, .vscode)
    % (兼容旧版 MATLAB 的写法)
    isHidden = cellfun(@(x) x(1)=='.', {items.name});
    items = items(~isHidden);
    
    % 3. (可选) 自定义忽略列表，例如忽略二进制输出目录
    ignoreList = {'bin', 'obj', 'DerivedData', '__pycache__'};
    items = items(~ismember({items.name}, ignoreList));
    
    % --- 排序优化 ---
    % 将文件夹排在文件前面 (可选，取消注释即可)
    % [~, sortIdx] = sort([items.isdir], 'descend');
    % items = items(sortIdx);
    
    count = length(items);
    
    for i = 1:count
        item = items(i);
        isLast = (i == count);
        
        % 确定连接符和下一级的缩进前缀
        if isLast
            connector = '└── ';
            nextPrefix = [prefix '    '];
        else
            connector = '├── ';
            nextPrefix = [prefix '│   '];
        end
        
        % 打印当前项
        if item.isdir
            fprintf('%s%s%s/\n', prefix, connector, item.name);
            % 递归进入子目录
            traverse_dir(fullfile(currentPath, item.name), nextPrefix);
        else
            fprintf('%s%s%s\n', prefix, connector, item.name);
        end
    end
end