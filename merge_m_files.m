%% MATLAB 代码合并工具
% 功能：选择一个文件夹，将其下所有 .m 文件合并到一个新的 .m 文件中
clear; clc;

%% 1. 选择源文件夹
disp('请选择包含 .m 文件的源文件夹...');
srcDir = uigetdir('', '选择源文件夹');

if srcDir == 0
    disp('已取消操作。');
    return;
end

%% 2. 选择保存路径和文件名
disp('请选择合并后文件的保存位置...');
[fileName, pathName] = uiputfile('*.m', '保存合并后的文件', 'All_Functions_Combined.m');

if fileName == 0
    disp('已取消操作。');
    return;
end

outputFilePath = fullfile(pathName, fileName);

%% 3. 获取所有 .m 文件
fileList = dir(fullfile(srcDir, '*.m'));
numFiles = length(fileList);

if numFiles == 0
    error('所选文件夹中没有找到任何 .m 文件！');
end

fprintf('正在处理，共找到 %d 个文件...\n', numFiles);

%% 4. 开始合并
fid_out = fopen(outputFilePath, 'w', 'n', 'UTF-8'); % 使用 UTF-8 打开写入，防止中文乱码

if fid_out == -1
    error('无法创建输出文件，请检查权限。');
end

% 写入文件头信息
% fprintf(fid_out, '%% ==================================================\n');
% fprintf(fid_out, '%% 合并生成时间: %s\n', char(datetime('now')));
% fprintf(fid_out, '%% 源文件夹: %s\n', srcDir);
% fprintf(fid_out, '%% ==================================================\n\n');

try
    for k = 1:numFiles
        thisFileName = fileList(k).name;
        
        % 跳过合并后的文件本身（如果它也被保存在源目录中）
        if strcmp(thisFileName, fileName)
            continue;
        end
        
        fullPath = fullfile(srcDir, thisFileName);
        
        % 读取源文件内容
        fileContent = fileread(fullPath);
        
        % 写入分隔符和文件名作为注释
        % fprintf(fid_out, '%% \n');
        fprintf(fid_out, '%% ==================================================\n');
        fprintf(fid_out, '%% === 文件开始: %s ===\n', thisFileName);
        fprintf(fid_out, '%% ==================================================\n');
        
        % 写入代码内容
        fprintf(fid_out, '%s', fileContent);
        
        % 确保末尾有换行，防止文件粘连
        if ~endsWith(fileContent, newline)
            fprintf(fid_out, '\n');
        end
        fprintf(fid_out, '\n\n'); % 额外加空行
        
        fprintf('已合并: %s\n', thisFileName);
    end
    
    disp('-----------------------------------');
    disp(['成功！所有文件已合并至: ' outputFilePath]);
    
catch ME
    % 错误处理
    fclose(fid_out);
    rethrow(ME);
end

% 关闭文件
fclose(fid_out);