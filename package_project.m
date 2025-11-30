function package_project()
% PACKAGE_PROJECT 自动将 MagPackage 项目打包为可分发的 ZIP 文件
%
% 功能:
%   1. 扫描项目目录，排除无关文件 (.git, .asv, .vscode 等)。
%   2. 生成带时间戳的发布包。
    
    % --- 配置 ---
    projectName = 'MagPackage';
    outputZipName = sprintf('%s_%s.zip', projectName, datestr(date));
    
    rootDir = pwd;
    
    fprintf('======================================================\n');
    fprintf('   %s 自动打包工具 \n', projectName);
    fprintf('======================================================\n');
    
    % 1. 扫描文件
    % ------------------------------------------------------
    fprintf('[Scan] 正在扫描项目文件...\n');
    allFiles = getAllFiles(rootDir);
    
    % 2. 过滤文件
    % ------------------------------------------------------
    fprintf('[Filter] 正在过滤无关文件...\n');
    filesToZip = filterFiles(allFiles, rootDir, outputZipName);
    
    if isempty(filesToZip)
        fprintf('[Error] 未找到可打包的文件。\n');
        return;
    end
    
    fprintf('   - 扫描到文件: %d 个\n', length(allFiles));
    fprintf('   - 实际打包:   %d 个\n', length(filesToZip));
    
    % 3. 生成 ZIP
    % ------------------------------------------------------
    fprintf('[Zip] 正在生成压缩包: %s ...\n', outputZipName);
    try
        % 使用相对路径打包，这样解压后结构正确
        zip(outputZipName, filesToZip);
        
        fileInfo = dir(outputZipName);
        fileSizeMB = fileInfo.bytes / 1024 / 1024;
        
        fprintf('\n[Success] 打包完成！\n');
        fprintf('   - 文件名: %s\n', outputZipName);
        fprintf('   - 大小:   %.2f MB\n', fileSizeMB);
        fprintf('   - 位置:   %s\n', fullfile(rootDir, outputZipName));
        
    catch ME
        fprintf('[Error] 打包失败: %s\n', ME.message);
    end
end

% --- 辅助函数 ---

function fileList = getAllFiles(dirName)
    % 递归获取所有文件路径
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    
    % 获取当前层级的文件
    fileList = {dirData(~dirIndex).name}';
    if ~isempty(fileList)
        % 转换为完整路径
        fileList = cellfun(@(f) fullfile(dirName, f), fileList, 'UniformOutput', false);
    end
    
    % 递归遍历子目录
    subDirs = {dirData(dirIndex).name};
    validIndex = ~ismember(subDirs, {'.', '..'});
    
    for i = find(validIndex)
        nextDir = fullfile(dirName, subDirs{i});
        fileList = [fileList; getAllFiles(nextDir)];
    end
end

function filtered = filterFiles(files, rootDir, outputZipName)
    % 过滤掉不需要打包的文件和文件夹
    
    % 定义排除规则 (正则表达式)
    excludePatterns = {
        '[/\\]\.git[/\\]', ...      % .git 文件夹内容
        '[/\\]\.git$', ...          % .git 文件夹本身
        '[/\\]\.vscode[/\\]', ...   % .vscode 文件夹
        '[/\\]\.idea[/\\]', ...     % .idea 文件夹
        '\.DS_Store', ...           % Mac 系统文件
        '\.asv$', ...               % MATLAB 自动保存文件
        '\.autosave$', ...          % Simulink 自动保存
        '~$' , ...                  % 临时备份文件
        'unnamed\.png', ...         % 临时截图
        ['^' strrep(outputZipName, '.', '\.') '$'], ... % 排除正在生成的 zip 本身
        '\.gitignore'               % git文件
    };
    
    keepMask = true(size(files));
    
    % 规范化根目录路径以确保替换正确
    if rootDir(end) ~= filesep, rootDir = [rootDir filesep]; end
    
    for i = 1:length(files)
        filePath = files{i};
        
        % 1. 匹配排除规则
        for p = 1:length(excludePatterns)
            if ~isempty(regexp(filePath, excludePatterns{p}, 'once'))
                keepMask(i) = false;
                break;
            end
        end
        
        % 2. 转换为相对路径 (用于 zip)
        if keepMask(i)
            % 移除根目录前缀，保留相对结构 (如 src/Assembler.m)
            files{i} = strrep(filePath, rootDir, '');
        end
    end
    
    filtered = files(keepMask);
end