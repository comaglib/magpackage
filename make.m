function make(varargin)
% MAKE 自动编译 MagPackage 的所有 C++ MEX 内核
% (自动适配 MinGW64 / MSVC / GCC / Clang)
%
% 修复日志:
%   - [Critical] 增加文件解锁逻辑 (clear mexName)。修复 Windows 下 "Permission denied" 错误。
%   - 保持 -R2017b 兼容模式。
%
% 用法:
%   make        - 编译所有文件
%   make clean  - 删除已编译的 mex 文件

    % --- 1. 配置路径 ---
    baseDir = pwd;
    srcDir = fullfile(baseDir, 'src', 'mex');
    outDir = fullfile(pwd, 'bin');
    
    if ~exist(srcDir, 'dir')
        error('Source directory not found: %s', srcDir);
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % --- 2. 处理清理命令 ---
    if nargin > 0 && strcmpi(varargin{1}, 'clean')
        fprintf('Cleaning MEX files...\n');
        delete(fullfile(outDir, '*.mex*'));
        fprintf('Done.\n');
        return;
    end

    % --- 3. 定义要编译的文件列表 ---
    mexFiles = { ...
        % 'assemble_curl_curl_kernel_mex.cpp', ...
        % 'assemble_mass_kernel_mex.cpp', ...
        % 'assemble_source_kernel_mex.cpp', ...
        % 'assemble_coupling_kernel_mex.cpp', ...
        % 'assemble_scalar_laplacian_kernel_mex.cpp', ...
        % 'assemble_jacobian_kernel_mex.cpp', ...
        % 'assemble_winding_kernel_mex.cpp', ...
        'assemble_hbfem_kernel_mex.cpp' ...
    };

    % --- 4. 设置编译器标志 ---
    includeFlag = ['-I' srcDir];
    
    % 使用兼容 API (-R2017b) 以支持 mxGetPi
    apiFlag = '-R2017b';
    
    % 检测编译器类型
    isMinGW = false;
    try
        cc = mex.getCompilerConfigurations('C++', 'Selected');
        if ~isempty(cc)
            if contains(cc.Name, 'MinGW', 'IgnoreCase', true) || ...
               contains(cc.Manufacturer, 'GNU', 'IgnoreCase', true) || ...
               contains(cc.ShortName, 'mingw', 'IgnoreCase', true)
               isMinGW = true;
            end
        end
    catch
    end

    if ispc && ~isMinGW
        % Windows (MSVC)
        cxxFlags = 'CXXFLAGS="$CXXFLAGS /O2 /openmp /std:c++11"';
        ldFlags = 'LDFLAGS="$LDFLAGS /openmp"';
        statusMsg = 'Building for Windows (MSVC)...';
    else
        % MinGW64 / Linux / Mac
        cxxFlags = 'CXXFLAGS="$CXXFLAGS -O3 -fopenmp -std=c++11"';
        ldFlags = 'LDFLAGS="$LDFLAGS -fopenmp"';
        if ispc
            statusMsg = 'Building for Windows (MinGW64)...';
        else
            cxxFlags = [cxxFlags ' -fPIC'];
            statusMsg = 'Building for Linux/Mac (GCC/Clang)...';
        end
    end

    fprintf('==================================================\n');
    fprintf('   MagPackage MEX Builder (Unlock Safe)\n');
    fprintf('==================================================\n');
    fprintf('%s\n', statusMsg);
    fprintf('Output Directory: %s\n', outDir);
    fprintf('API Mode: %s\n', apiFlag);
    fprintf('--------------------------------------------------\n');

    % --- 5. 编译循环 ---
    successCount = 0;
    tTotal = tic;

    for k = 1:length(mexFiles)
        fileName = mexFiles{k};
        filePath = fullfile(srcDir, fileName);
        [~, mexName, ~] = fileparts(fileName);
        
        if ~exist(filePath, 'file')
            fprintf('[SKIP] File not found: %s\n', fileName);
            continue;
        end
        
        fprintf('Compiling: %-40s ... ', fileName);
        tStart = tic;
        
        % [关键修复] 在编译前强制清除内存中的函数，解除 Windows 文件锁
        clear(mexName); 
        
        try
            mex(apiFlag, ...
                '-outdir', outDir, ...
                includeFlag, ...
                cxxFlags, ...
                ldFlags, ...
                filePath);
            
            tElapsed = toc(tStart);
            fprintf('[OK] (%.2fs)\n', tElapsed);
            successCount = successCount + 1;
            
        catch ME
            fprintf('\n[FAILED]\n');
            fprintf('%s\n', ME.message);
        end
    end

    fprintf('--------------------------------------------------\n');
    fprintf('Build Complete. Success: %d / %d. Total Time: %.2fs\n', ...
        successCount, length(mexFiles), toc(tTotal));
end