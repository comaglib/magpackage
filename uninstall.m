function uninstall
    baseDir = fileparts(mfilename('fullpath'));
    rmpath(genpath(fullfile(baseDir, 'src')));
    rmpath(genpath(fullfile(baseDir, 'tests')));
    rmpath(genpath(fullfile(baseDir, 'tutorials')));
    rmpath(fullfile(baseDir, 'bin')); 
    savepath; clc;
    fprintf('MagPackage uninstalled successfully.\n');
end