function make()
    % make.m - Compile ALL MEX files including Nonlinear Kernels
    srcDir = fullfile(pwd, 'src', 'mex');
    outDir = fullfile(pwd, 'bin');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    
    files = { ...
        'assemble_coupling_kernel_mex.cpp', ...
        'assemble_curl_curl_kernel_mex.cpp', ...
        'assemble_mass_kernel_mex.cpp', ...
        'assemble_source_kernel_mex.cpp', ...
        'assemble_winding_kernel_mex.cpp', ...
        'assemble_scalar_laplacian_kernel_mex.cpp', ...
        'assemble_jacobian_kernel_mex.cpp', ... 
        'assemble_hbfem_kernel_mex.cpp' ... 
        };
    
    % 检查头文件是否存在
    if ~exist(fullfile(srcDir, 'MexMaterialUtils.hpp'), 'file')
        error('MexMaterialUtils.hpp is missing in src/mex/');
    end
    
    if ispc
        cxxFlags = 'COMPFLAGS="$COMPFLAGS /openmp /O2"';
        ldFlags  = '';
    else
        cxxFlags = 'CXXFLAGS="$CXXFLAGS -fopenmp -O3 -march=native"';
        ldFlags  = 'LDFLAGS="$LDFLAGS -fopenmp"';
    end
    
    fprintf('Compiling MEX files...\n');
    for i = 1:length(files)
        fprintf('  Building %s ... ', files{i});
        sourceFile = fullfile(srcDir, files{i});
        try
            mex(sourceFile, '-outdir', outDir, cxxFlags, ldFlags);
            fprintf('Done.\n');
        catch ME
            fprintf('Failed.\n%s\n', ME.message);
        end
    end
    addpath(outDir);
    fprintf('Compilation finished.\n');
end