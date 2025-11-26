function build_distribution()
% BUILD_DISTRIBUTION 项目打包脚本
% 
% 功能:
% 1. 创建 build 目录
% 2. 复制 src, data, docs, examples
% 3. 清理临时文件 (.asv, .mex, etc.)
% 4. 生成版本信息

    fprintf('==============================================\n');
    fprintf('   Building Distribution Package              \n');
    fprintf('==============================================\n');
    
    version = 'v1.0_alpha';
    dist_name = ['magpackage_' version];
    
    % 1. 创建目录
    if exist(dist_name, 'dir')
        rmdir(dist_name, 's');
    end
    mkdir(dist_name);
    
    % 2. 复制文件
    fprintf('  - Copying Source Code...\n');
    copyfile('src', fullfile(dist_name, 'src'));
    
    fprintf('  - Copying Data...\n');
    copyfile('data', fullfile(dist_name, 'data'));
    
    fprintf('  - Copying Documents...\n');
    copyfile('docs', fullfile(dist_name, 'docs'));
    
    fprintf('  - Copying Examples...\n');
    copyfile('examples', fullfile(dist_name, 'examples'));
    
    % 3. 创建启动脚本
    fid = fopen(fullfile(dist_name, 'startup_fem.m'), 'w');
    fprintf(fid, '%% Startup script\n');
    fprintf(fid, 'addpath(genpath(''src''));\n');
    fprintf(fid, 'addpath(genpath(''data''));\n');
    fprintf(fid, 'addpath(genpath(''docs''));\n');
    fprintf(fid, 'addpath(genpath(''examples''));\n');
    fprintf(fid, 'fprintf(''MAGPACKAGE Loaded.\\n'');\n');
    fclose(fid);
    
    % 3. 创建清理脚本
    fid = fopen(fullfile(dist_name, 'cleanup_fem.m'), 'w');
    fprintf(fid, '%% Cleanup script\n');
    fprintf(fid, 'rmpath(genpath(''src''));\n');
    fprintf(fid, 'rmpath(genpath(''data''));\n');
    fprintf(fid, 'rmpath(genpath(''docs''));\n');
    fprintf(fid, 'rmpath(genpath(''examples''));\n');
    fprintf(fid, 'fprintf(''MAGPACKAGE Cleaned.\\n'');\n');
    fclose(fid);
    
    % 5. 清理垃圾文件
    fprintf('  - Cleaning up...\n');
    files = dir(fullfile(dist_name, '**', '*.*'));
    for i = 1:length(files)
        [~, ~, ext] = fileparts(files(i).name);
        if strcmp(ext, '.asv') || strcmp(ext, '.autosave')
            delete(fullfile(files(i).folder, files(i).name));
        end
    end
    
    fprintf('  [Success] Build created in folder: %s\n', dist_name);
end