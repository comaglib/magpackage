classdef FileUtils
    % FILEUTILS 文件操作辅助工具类
    % 提供文件查找、路径解析等通用功能
    
    methods (Static)
        function fullPath = locateFile(filename, searchDirs)
            % LOCATEFILE 在路径和指定文件夹中查找文件，并返回绝对路径
            %
            % 输入:
            %   filename   - 目标文件名 (例如 'Team7.mphtxt')
            %   searchDirs - (可选) 额外的搜索目录列表 (Cell array of strings)
            %
            % 输出:
            %   fullPath   - 文件的绝对路径。如果找不到则报错。
            
            if nargin < 2, searchDirs = {}; end
            if ischar(searchDirs), searchDirs = {searchDirs}; end
            
            fullPath = '';
            
            % 1. 尝试使用 MATLAB 内置的 which 查找 (覆盖当前目录和 Path)
            onPath = which(filename);
            if ~isempty(onPath)
                fullPath = onPath;
                return;
            end
            
            % 2. 检查当前目录及指定的相对路径
            % 组合候选路径列表
            candidates = {filename}; % 首先检查当前工作目录
            for i = 1:length(searchDirs)
                candidates{end+1} = fullfile(searchDirs{i}, filename); %#ok<AGROW>
            end
            
            % 遍历查找
            for i = 1:length(candidates)
                target = candidates{i};
                if exist(target, 'file')
                    % 获取文件的绝对路径信息
                    info = dir(target);
                    if ~isempty(info)
                        fullPath = fullfile(info(1).folder, info(1).name);
                        return;
                    end
                end
            end
            
            % 3. 如果仍未找到，抛出错误
            if isempty(fullPath)
                error('FileUtils:FileNotFound', ...
                    'File "%s" not found in MATLAB path or specified search directories.', filename);
            end
        end
    end
end